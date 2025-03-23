# cython: cpow=True

import re
import warnings
import numpy as np
from scipy._lib._util import check_random_state, _transition_to_rng
from ._rotation_groups import create_group

cimport numpy as np
cimport cython
from cython.view cimport array
from libc.math cimport sqrt, sin, cos, atan2, acos, hypot, isnan, NAN, pi

np.import_array()

# utilities for empty array initialization
cdef inline double[:] _empty1(int n) noexcept:
    return array(shape=(n,), itemsize=sizeof(double), format=b"d")
cdef inline double[:, :] _empty2(int n1, int n2) noexcept :
    return array(shape=(n1, n2), itemsize=sizeof(double), format=b"d")
cdef inline double[:, :, :] _empty3(int n1, int n2, int n3) noexcept:
    return array(shape=(n1, n2, n3), itemsize=sizeof(double), format=b"d")
cdef inline double[:, :] _zeros2(int n1, int n2) noexcept:
    cdef double[:, :] arr = array(shape=(n1, n2),
        itemsize=sizeof(double), format=b"d")
    arr[:, :] = 0
    return arr

# flat implementations of numpy functions
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline double[:] _cross3(const double[:] a, const double[:] b) noexcept:
    cdef double[:] result = _empty1(3)
    result[0] = a[1]*b[2] - a[2]*b[1]
    result[1] = a[2]*b[0] - a[0]*b[2]
    result[2] = a[0]*b[1] - a[1]*b[0]
    return result

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline double _dot3(const double[:] a, const double[:] b) noexcept nogil:
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline double _norm3(const double[:] elems) noexcept nogil:
    return sqrt(_dot3(elems, elems))

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline double _normalize4(double[:] elems) noexcept nogil:
    cdef double norm = sqrt(_dot3(elems, elems) + elems[3]*elems[3])

    if norm == 0:
        return NAN

    elems[0] /= norm
    elems[1] /= norm
    elems[2] /= norm
    elems[3] /= norm

    return norm

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int _argmax4(const double[:] a) noexcept nogil:
    cdef int imax = 0
    cdef double vmax = a[0]

    for i in range(1, 4):
        if a[i] > vmax:
            imax = i
            vmax = a[i]

    return imax

ctypedef unsigned char uchar

cdef double[3] _ex = [1, 0, 0]
cdef double[3] _ey = [0, 1, 0]
cdef double[3] _ez = [0, 0, 1]

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline const double[:] _elementary_basis_vector(uchar axis) noexcept:
    if axis == b'x': return _ex
    elif axis == b'y': return _ey
    elif axis == b'z': return _ez

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int _elementary_basis_index(uchar axis) noexcept:
    if axis == b'x': return 0
    elif axis == b'y': return 1
    elif axis == b'z': return 2

# Reduce the quaternion double coverage of the rotation group to a unique
# canonical "positive" single cover
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline void _quat_canonical_single(double[:] q) noexcept nogil:
    if ((q[3] < 0)
        or (q[3] == 0 and q[0] < 0)
        or (q[3] == 0 and q[0] == 0 and q[1] < 0)
        or (q[3] == 0 and q[0] == 0 and q[1] == 0 and q[2] < 0)):
        q[0] *= -1.0
        q[1] *= -1.0
        q[2] *= -1.0
        q[3] *= -1.0

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline void _quat_canonical(double[:, :] q) noexcept:
    cdef Py_ssize_t n = q.shape[0]
    for ind in range(n):
        _quat_canonical_single(q[ind])

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline void _get_angles(
    double[:] angles, bint extrinsic, bint symmetric, bint sign,
    double lamb, double a, double b, double c, double d):

    # intrinsic/extrinsic conversion helpers
    cdef int angle_first, angle_third
    if extrinsic:
        angle_first = 0
        angle_third = 2
    else:
        angle_first = 2
        angle_third = 0

    cdef double half_sum, half_diff
    cdef int case

    # Step 2
    # Compute second angle...
    angles[1] = 2 * atan2(hypot(c, d), hypot(a, b))

    # ... and check if equal to is 0 or pi, causing a singularity
    if abs(angles[1]) <= 1e-7:
        case = 1
    elif abs(angles[1] - <double>pi) <= 1e-7:
        case = 2
    else:
        case = 0 # normal case

    # Step 3
    # compute first and third angles, according to case
    half_sum = atan2(b, a)
    half_diff = atan2(d, c)

    if case == 0:  # no singularities
        angles[angle_first] = half_sum - half_diff
        angles[angle_third] = half_sum + half_diff

    else:  # any degenerate case
        angles[2] = 0
        if case == 1:
            angles[0] = 2 * half_sum
        else:
            angles[0] = 2 * half_diff * (-1 if extrinsic else 1)

    # for Tait-Bryan/asymmetric sequences
    if not symmetric:
        angles[angle_third] *= sign
        angles[1] -= lamb

    for idx in range(3):
        if angles[idx] < -pi:
            angles[idx] += 2 * pi
        elif angles[idx] > pi:
            angles[idx] -= 2 * pi

    if case != 0:
        warnings.warn("Gimbal lock detected. Setting third angle to zero "
                      "since it is not possible to uniquely determine "
                      "all angles.", stacklevel=3)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:, :] _compute_davenport_from_quat(
    np.ndarray[double, ndim=2] quat, np.ndarray[double, ndim=1] n1,
    np.ndarray[double, ndim=1] n2, np.ndarray[double, ndim=1] n3,
    bint extrinsic
):
    # The algorithm assumes extrinsic frame transformations. The algorithm
    # in the paper is formulated for rotation quaternions, which are stored
    # directly by Rotation.
    # Adapt the algorithm for our case by reversing both axis sequence and
    # angles for intrinsic rotations when needed

    if not extrinsic:
        n1, n3 = n3, n1

    cdef double[:] n_cross = _cross3(n1, n2)
    cdef double lamb = atan2(_dot3(n3, n_cross), _dot3(n3, n1))

    cdef int correct_set = False
    if lamb < 0:
        # alternative set of angles compatible with as_euler implementation
        n2 = -n2
        lamb = -lamb
        n_cross[0] = -n_cross[0]
        n_cross[1] = -n_cross[1]
        n_cross[2] = -n_cross[2]
        correct_set = True

    cdef double[:] quat_lamb = np.array([
            sin(lamb / 2) * n2[0],
            sin(lamb / 2) * n2[1],
            sin(lamb / 2) * n2[2],
            cos(lamb / 2)]
    )

    cdef Py_ssize_t num_rotations = quat.shape[0]

    # some forward definitions
    cdef double[:, :] angles = _empty2(num_rotations, 3)
    cdef double[:] quat_transformed = _empty1(4)
    cdef double a, b, c, d

    for ind in range(num_rotations):
        _compose_quat_single(quat_lamb, quat[ind], quat_transformed)

        # Step 1
        # Permutate quaternion elements
        a = quat_transformed[3]
        b = _dot3(quat_transformed[:3], n1)
        c = _dot3(quat_transformed[:3], n2)
        d = _dot3(quat_transformed[:3], n_cross)

        _get_angles(angles[ind], extrinsic, False, 1, lamb, a, b, c, d)

        if correct_set:
            angles[ind, 1] = -angles[ind, 1]

    return angles

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline void _compose_quat_single( # calculate p * q into r
    const double[:] p, const double[:] q, double[:] r
) noexcept:
    cdef double[:] cross = _cross3(p[:3], q[:3])

    r[0] = p[3]*q[0] + q[3]*p[0] + cross[0]
    r[1] = p[3]*q[1] + q[3]*p[1] + cross[1]
    r[2] = p[3]*q[2] + q[3]*p[2] + cross[2]
    r[3] = p[3]*q[3] - p[0]*q[0] - p[1]*q[1] - p[2]*q[2]

@cython.boundscheck(False)
@cython.wraparound(False)
def compose_quat(
    const double[:, :] p, const double[:, :] q
) -> double[:, :]:
    return _compose_quat(p, q)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline double[:, :] _compose_quat(
    const double[:, :] p, const double[:, :] q
) noexcept:
    cdef Py_ssize_t n = max(p.shape[0], q.shape[0])
    cdef double[:, :] product = _empty2(n, 4)

    # dealing with broadcasting
    if p.shape[0] == 1:
        for ind in range(n):
            _compose_quat_single(p[0], q[ind], product[ind])
    elif q.shape[0] == 1:
        for ind in range(n):
            _compose_quat_single(p[ind], q[0], product[ind])
    else:
        for ind in range(n):
            _compose_quat_single(p[ind], q[ind], product[ind])

    return product

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline double[:, :] _make_elementary_quat(
    uchar axis, const double[:] angles
) noexcept:
    cdef Py_ssize_t n = angles.shape[0]
    cdef double[:, :] quat = _zeros2(n, 4)

    cdef int axis_ind
    if axis == b'x':   axis_ind = 0
    elif axis == b'y': axis_ind = 1
    elif axis == b'z': axis_ind = 2

    for ind in range(n):
        quat[ind, 3] = cos(angles[ind] / 2)
        quat[ind, axis_ind] = sin(angles[ind] / 2)
    return quat

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:, :] _elementary_quat_compose(
    const uchar[:] seq, const double[:, :] angles, bint intrinsic=False
) noexcept:
    cdef double[:, :] result = _make_elementary_quat(seq[0], angles[:, 0])
    cdef Py_ssize_t seq_len = seq.shape[0]

    for idx in range(1, seq_len):
        if intrinsic:
            result = _compose_quat(
                result,
                _make_elementary_quat(seq[idx], angles[:, idx]))
        else:
            result = _compose_quat(
                _make_elementary_quat(seq[idx], angles[:, idx]),
                result)
    return result

def _format_angles(angles, degrees, num_axes):
    angles = np.asarray(angles, dtype=float)
    if degrees:
        angles = np.deg2rad(angles)

    is_single = False
    # Prepare angles to have shape (num_rot, num_axes)
    if num_axes == 1:
        if angles.ndim == 0:
            # (1, 1)
            angles = angles.reshape((1, 1))
            is_single = True
        elif angles.ndim == 1:
            # (N, 1)
            angles = angles[:, None]
        elif angles.ndim == 2 and angles.shape[-1] != 1:
            raise ValueError("Expected `angles` parameter to have shape "
                             "(N, 1), got {}.".format(angles.shape))
        elif angles.ndim > 2:
            raise ValueError("Expected float, 1D array, or 2D array for "
                             "parameter `angles` corresponding to `seq`, "
                             "got shape {}.".format(angles.shape))
    else:  # 2 or 3 axes
        if angles.ndim not in [1, 2] or angles.shape[-1] != num_axes:
            raise ValueError("Expected `angles` to be at most "
                             "2-dimensional with width equal to number "
                             "of axes specified, got "
                             "{} for shape".format(angles.shape))

        if angles.ndim == 1:
            # (1, num_axes)
            angles = angles[None, :]
            is_single = True

    # By now angles should have shape (num_rot, num_axes)
    # sanity check
    if angles.ndim != 2 or angles.shape[-1] != num_axes:
        raise ValueError("Expected angles to have shape (num_rotations, "
                         "num_axes), got {}.".format(angles.shape))

    return angles, is_single

@cython.boundscheck(False)
@cython.wraparound(False)
def from_quat(quat: cython.double[:, :], normalize: cython.bint = True, copy: cython.bint = True, scalar_first: cython.bint = False):
    if quat.ndim != 2 or quat.shape[1] != 4 or quat.shape[0] == 0:
        raise ValueError(f"Expected `quat` to have shape (N, 4), got {quat.shape}.")

    # TODO: We always assume quat is a 2D array with shape (N, 4). Instead, we should broadcast over
    # all dimensions.

    cdef Py_ssize_t num_rotations = quat.shape[0]

    if scalar_first:
        quat = np.roll(quat, -1, axis=1)
    elif normalize or copy:
        quat = quat.copy()

    if normalize:
        for ind in range(num_rotations):
            if isnan(_normalize4(quat[ind, :])):
                raise ValueError("Found zero norm quaternions in `quat`.")
        # TODO: A broadcasted version of this would be faster for rotations > 100 and work for any
        # number of leading dimensions.
        # quat_norm = np.linalg.norm(quat, axis=-1, keepdims=True)
        # if np.any(quat_norm == 0):
        #     raise ValueError("Found zero norm quaternions in `quat`.")
        # quat /= quat_norm

    return np.asarray(quat, dtype=float)

@cython.boundscheck(False)
@cython.wraparound(False)
def from_euler(seq, angles, degrees=False):
    num_axes = len(seq)
    if num_axes < 1 or num_axes > 3:
        raise ValueError("Expected axis specification to be a non-empty "
                            "string of upto 3 characters, got {}".format(seq))

    intrinsic = (re.match(r'^[XYZ]{1,3}$', seq) is not None)
    extrinsic = (re.match(r'^[xyz]{1,3}$', seq) is not None)
    if not (intrinsic or extrinsic):
        raise ValueError("Expected axes from `seq` to be from ['x', 'y', "
                            "'z'] or ['X', 'Y', 'Z'], got {}".format(seq))

    if any(seq[i] == seq[i+1] for i in range(num_axes - 1)):
        raise ValueError("Expected consecutive axes to be different, "
                            "got {}".format(seq))

    seq = seq.lower()

    angles, is_single = _format_angles(angles, degrees, num_axes)

    quat = _elementary_quat_compose(seq.encode(), angles, intrinsic)

    if is_single:
        return quat[0]
    return quat

@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def from_matrix(matrix) -> double[:, :]:
    is_single = False
    matrix = np.array(matrix, dtype=float)

    if (matrix.ndim not in [2, 3] or
        matrix.shape[len(matrix.shape)-2:] != (3, 3)):
        raise ValueError("Expected `matrix` to have shape (3, 3) or "
                            "(N, 3, 3), got {}".format(matrix.shape))

    # If a single matrix is given, convert it to 3D 1 x 3 x 3 matrix but
    # set is_single to True so that we can return appropriate objects in
    # the `to_...` methods
    if matrix.shape == (3, 3):
        matrix = matrix[np.newaxis, :, :]
        is_single = True

    # Calculate the determinant of the rotation matrix
    # (should be positive for right-handed rotations)
    dets = np.linalg.det(matrix)
    if np.any(dets <= 0):
        ind = np.where(dets <= 0)[0][0]
        raise ValueError("Non-positive determinant (left-handed or null "
                            f"coordinate frame) in rotation matrix {ind}: "
                            f"{matrix[ind]}.")

    # Gramian orthogonality check
    # (should be the identity matrix for orthogonal matrices)
    # Note that we have already ruled out left-handed cases above
    gramians = matrix @ np.transpose(matrix, (0, 2, 1))
    is_orthogonal = np.all(np.isclose(gramians, np.eye(3), atol=1e-12),
                            axis=(1, 2))
    indices_to_orthogonalize = np.where(~is_orthogonal)[0]

    # Orthogonalize the rotation matrices where necessary
    if len(indices_to_orthogonalize) > 0:
        # Exact solution to the orthogonal Procrustes problem using singular
        # value decomposition
        U, _, Vt = np.linalg.svd(matrix[indices_to_orthogonalize, :, :])
        matrix[indices_to_orthogonalize, :, :] = U @ Vt

    # Convert the orthogonal rotation matrices to quaternions using the
    # algorithm described in [3]_. This will also apply another
    # orthogonalization step to correct for any small errors in the matrices
    # that skipped the SVD step above.
    cdef double[:, :, :] cmatrix
    cmatrix = matrix
    cdef Py_ssize_t num_rotations = cmatrix.shape[0]
    cdef Py_ssize_t i, j, k
    cdef double[:] decision = _empty1(4)
    cdef int choice

    cdef double[:, :] quat = _empty2(num_rotations, 4)

    for ind in range(num_rotations):
        decision[0] = cmatrix[ind, 0, 0]
        decision[1] = cmatrix[ind, 1, 1]
        decision[2] = cmatrix[ind, 2, 2]
        decision[3] = cmatrix[ind, 0, 0] + cmatrix[ind, 1, 1] \
                    + cmatrix[ind, 2, 2]
        choice = _argmax4(decision)

        if choice != 3:
            i = choice
            j = (i + 1) % 3
            k = (j + 1) % 3

            quat[ind, i] = 1 - decision[3] + 2 * cmatrix[ind, i, i]
            quat[ind, j] = cmatrix[ind, j, i] + cmatrix[ind, i, j]
            quat[ind, k] = cmatrix[ind, k, i] + cmatrix[ind, i, k]
            quat[ind, 3] = cmatrix[ind, k, j] - cmatrix[ind, j, k]
        else:
            quat[ind, 0] = cmatrix[ind, 2, 1] - cmatrix[ind, 1, 2]
            quat[ind, 1] = cmatrix[ind, 0, 2] - cmatrix[ind, 2, 0]
            quat[ind, 2] = cmatrix[ind, 1, 0] - cmatrix[ind, 0, 1]
            quat[ind, 3] = 1 + decision[3]

        # normalize
        _normalize4(quat[ind])

    if is_single:
        return quat[0]
    return quat


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def from_rotvec(rotvec, degrees: cython.bint = False) -> double[:, :]:
    is_single = False
    rotvec = np.asarray(rotvec, dtype=float)
    if degrees:
        rotvec = np.deg2rad(rotvec)

    if rotvec.ndim not in [1, 2] or rotvec.shape[len(rotvec.shape)-1] != 3:
        raise ValueError("Expected `rot_vec` to have shape (3,) "
                            "or (N, 3), got {}".format(rotvec.shape))

    # If a single vector is given, convert it to a 2D 1 x 3 matrix but
    # set self._single to True so that we can return appropriate objects
    # in the `as_...` methods
    cdef double[:, :] crotvec
    if rotvec.shape == (3,):
        crotvec = rotvec[None, :]
        is_single = True
    else:
        crotvec = rotvec

    cdef Py_ssize_t num_rotations = crotvec.shape[0]
    cdef double angle, scale, angle2
    cdef double[:, :] quat = _empty2(num_rotations, 4)

    for ind in range(num_rotations):
        angle = _norm3(crotvec[ind, :])

        if angle <= 1e-3:  # small angle Taylor series expansion
            angle2 = angle * angle
            scale = 0.5 - angle2 / 48 + angle2 * angle2 / 3840
        else:  # large angle
            scale = sin(angle / 2) / angle

        quat[ind, 0] = scale * crotvec[ind, 0]
        quat[ind, 1] = scale * crotvec[ind, 1]
        quat[ind, 2] = scale * crotvec[ind, 2]
        quat[ind, 3] = cos(angle / 2)

    if is_single:
        return quat[0]
    return quat

@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def from_mrp(mrp) -> double[:, :]:
    is_single = False
    mrp = np.asarray(mrp, dtype=float)

    if mrp.ndim not in [1, 2] or mrp.shape[len(mrp.shape) - 1] != 3:
        raise ValueError("Expected `mrp` to have shape (3,) "
                            "or (N, 3), got {}".format(mrp.shape))

    # If a single vector is given, convert it to a 2D 1 x 3 matrix but
    # set self._single to True so that we can return appropriate objects
    # in the `as_...` methods
    cdef double[:, :] cmrp
    if mrp.shape == (3,):
        cmrp = mrp[None, :]
        is_single = True
    else:
        cmrp = mrp

    cdef Py_ssize_t num_rotations = cmrp.shape[0]
    cdef double[:, :] quat = _empty2(num_rotations, 4)
    cdef double mrp_squared_plus_1

    for ind in range(num_rotations):
        mrp_squared_plus_1 = 1 + _dot3(cmrp[ind, :], cmrp[ind, :])

        quat[ind, 0] = 2 * cmrp[ind, 0] / mrp_squared_plus_1
        quat[ind, 1] = 2 * cmrp[ind, 1] / mrp_squared_plus_1
        quat[ind, 2] = 2 * cmrp[ind, 2] / mrp_squared_plus_1
        quat[ind, 3] = (2 - mrp_squared_plus_1) / mrp_squared_plus_1

    if is_single:
        return quat[0]
    return quat
    


@cython.boundscheck(False)
@cython.wraparound(False)
def as_quat(quat: cython.double[:, :], normalize: cython.bint = True, copy: cython.bint = True, canonical: cython.bint = False, scalar_first: cython.bint = False):
    if quat.ndim != 2 or quat.shape[1] != 4 or quat.shape[0] == 0:
        raise ValueError(f"Expected `quat` to have shape (N, 4), got {quat.shape}.")
    # If a single quaternion is given, convert it to a 2D 1 x 4 matrix but
    # set self._single to True so that we can return appropriate objects
    # in the `to_...` methods
    # TODO: We always assume quat is a 2D array with shape (N, 4). Instead, we should broadcast over
    # all dimensions.

    cdef Py_ssize_t num_rotations = quat.shape[0]

    if canonical:
        _quat_canonical(quat)

    if scalar_first:
        quat = np.roll(quat, 1, axis=1)
    elif normalize or copy:
        quat = quat.copy()

    if normalize:
        for ind in range(num_rotations):
            if isnan(_normalize4(quat[ind, :])):
                raise ValueError("Found zero norm quaternions in `quat`.")

    return np.asarray(quat, dtype=float)

@cython.boundscheck(False)
@cython.wraparound(False)
def as_matrix(quat: cython.double[:, :]):

    cdef Py_ssize_t num_rotations = quat.shape[0]
    cdef double[:, :, :] matrix = _empty3(num_rotations, 3, 3)

    cdef double x, y, z, w, x2, y2, z2, w2
    cdef double xy, zw, xz, yw, yz, xw

    for ind in range(num_rotations):
        x = quat[ind, 0]
        y = quat[ind, 1]
        z = quat[ind, 2]
        w = quat[ind, 3]

        x2 = x * x
        y2 = y * y
        z2 = z * z
        w2 = w * w

        xy = x * y
        zw = z * w
        xz = x * z
        yw = y * w
        yz = y * z
        xw = x * w

        matrix[ind, 0, 0] = x2 - y2 - z2 + w2
        matrix[ind, 1, 0] = 2 * (xy + zw)
        matrix[ind, 2, 0] = 2 * (xz - yw)

        matrix[ind, 0, 1] = 2 * (xy - zw)
        matrix[ind, 1, 1] = - x2 + y2 - z2 + w2
        matrix[ind, 2, 1] = 2 * (yz + xw)

        matrix[ind, 0, 2] = 2 * (xz + yw)
        matrix[ind, 1, 2] = 2 * (yz - xw)
        matrix[ind, 2, 2] = - x2 - y2 + z2 + w2

    return np.asarray(matrix)


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def as_rotvec(quat: cython.double[:, :], degrees: cython.bint = False) -> double[:, :]:
    cdef Py_ssize_t num_rotations = len(quat)
    cdef double angle, scale, angle2
    cdef double[:, :] rotvec = _empty2(num_rotations, 3)
    cdef double[:] quat_single

    for ind in range(num_rotations):
        quat_single = quat[ind, :].copy()
        _quat_canonical_single(quat_single)  # w > 0 ensures that 0 <= angle <= pi

        angle = 2 * atan2(_norm3(quat_single), quat_single[3])

        if angle <= 1e-3:  # small angle Taylor series expansion
            angle2 = angle * angle
            scale = 2 + angle2 / 12 + 7 * angle2 * angle2 / 2880
        else:  # large angle
            scale = angle / sin(angle / 2)

        rotvec[ind, 0] = scale * quat_single[0]
        rotvec[ind, 1] = scale * quat_single[1]
        rotvec[ind, 2] = scale * quat_single[2]

    if degrees:
        rotvec = np.rad2deg(rotvec)

    return np.asarray(rotvec)


@cython.embedsignature(True)
def as_mrp(quat: cython.double[:, :]) -> cython.double[:, :]:
    cdef Py_ssize_t num_rotations = len(quat)
    cdef double[:, :] mrps = _empty2(num_rotations, 3)
    cdef int sign
    cdef double denominator

    for ind in range(num_rotations):

        # Ensure we are calculating the set of MRPs that correspond
        # to a rotation of <= 180
        sign = -1 if quat[ind, 3] < 0 else 1

        denominator = 1 + sign * quat[ind, 3]
        for i in range(3):
            mrps[ind, i] = sign * quat[ind, i] / denominator

    return np.asarray(mrps)


@cython.embedsignature(True)
def as_euler(quat: cython.double[:, :], seq: str, degrees: cython.bint = False):
    # Prepare axis sequence to call Euler angles conversion algorithm.
    if len(seq) != 3:
        raise ValueError("Expected 3 axes, got {}.".format(seq))

    intrinsic = (re.match(r'^[XYZ]{1,3}$', seq) is not None)
    extrinsic = (re.match(r'^[xyz]{1,3}$', seq) is not None)
    if not (intrinsic or extrinsic):
        raise ValueError("Expected axes from `seq` to be from "
                            "['x', 'y', 'z'] or ['X', 'Y', 'Z'], "
                            "got {}".format(seq))

    if any(seq[i] == seq[i+1] for i in range(2)):
        raise ValueError("Expected consecutive axes to be different, "
                            "got {}".format(seq))

    seq = seq.lower()
    if not extrinsic:
        seq = seq[::-1]

    # The algorithm assumes extrinsic frame transformations. The algorithm
    # in the paper is formulated for rotation quaternions, which are stored
    # directly by Rotation.
    # Adapt the algorithm for our case by reversing both axis sequence and
    # angles for intrinsic rotations when needed

    cdef const uchar[:] cseq = seq.encode()
    cdef int i = _elementary_basis_index(cseq[0])
    cdef int j = _elementary_basis_index(cseq[1])
    cdef int k = _elementary_basis_index(cseq[2])

    cdef bint symmetric = i == k
    if symmetric:
        k = 3 - i - j # get third axis

    # Step 0
    # Check if permutation is even (+1) or odd (-1)
    cdef int sign = (i - j) * (j - k) * (k - i) // 2

    cdef Py_ssize_t num_rotations = quat.shape[0]

    # some forward definitions
    cdef double a, b, c, d
    cdef double[:, :] angles = _empty2(num_rotations, 3)

    for ind in range(num_rotations):

        # Step 1
        # Permutate quaternion elements
        if symmetric:
            a = quat[ind, 3]
            b = quat[ind, i]
            c = quat[ind, j]
            d = quat[ind, k] * sign
        else:
            a = quat[ind, 3] - quat[ind, j]
            b = quat[ind, i] + quat[ind, k] * sign
            c = quat[ind, j] + quat[ind, 3]
            d = quat[ind, k] * sign - quat[ind, i]

        _get_angles(angles[ind], extrinsic, symmetric, sign, pi / 2, a, b, c, d)

    if degrees:
        angles = np.rad2deg(angles)
    return np.asarray(angles)


@cython.embedsignature(True)
def as_davenport(quat: double[:, :], axes, order: str, degrees: cython.bint = False) -> double[:, :]:
    if order in ['e', 'extrinsic']:
        extrinsic = True
    elif order in ['i', 'intrinsic']:
        extrinsic = False
    else:
        raise ValueError("order should be 'e'/'extrinsic' for extrinsic "
                            "sequences or 'i'/'intrinsic' for intrinsic "
                            "sequences, got {}".format(order))

    if len(axes) != 3:
        raise ValueError("Expected 3 axes, got {}.".format(len(axes)))

    axes = np.asarray(axes)

    if axes.shape[1] != 3:
        raise ValueError("Axes must be vectors of length 3.")

    n1, n2, n3 = axes

    # normalize axes
    n1 = n1 / np.linalg.norm(n1)
    n2 = n2 / np.linalg.norm(n2)
    n3 = n3 / np.linalg.norm(n3)

    if np.dot(n1, n2) >= 1e-7:
        raise ValueError("Consecutive axes must be orthogonal.")
    if np.dot(n2, n3) >= 1e-7:
        raise ValueError("Consecutive axes must be orthogonal.")

    angles = np.asarray(_compute_davenport_from_quat(
            quat, n1, n2, n3, extrinsic))

    if degrees:
        angles = np.rad2deg(angles)

    return angles


@cython.embedsignature(True)
@cython.boundscheck(False)
def inv(quat: double[:, :]) -> double[:, :]:
    cdef np.ndarray q_inv = np.array(quat, copy=True)
    q_inv[:, 0] *= -1
    q_inv[:, 1] *= -1
    q_inv[:, 2] *= -1
    return q_inv


@cython.embedsignature(True)
@_transition_to_rng('random_state', position_num=1)
def random(num=None, rng=None):
    rng = check_random_state(rng)

    if num is None:
        sample = rng.normal(size=4)
    else:
        sample = rng.normal(size=(num, 4))

    return sample


@cython.embedsignature(True)
def identity(num: cython.int | None = None) -> double[:, :]:
    if num is None:
        return np.array([0, 0, 0, 1])
    q = np.zeros((num, 4))
    q[:, 3] = 1
    return q


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def magnitude(quat: cython.double[:, :]) -> double[:, :]:
    cdef Py_ssize_t num_rotations = quat.shape[0]
    cdef double[:] angles = _empty1(num_rotations)

    for ind in range(num_rotations):
        angles[ind] = 2 * atan2(_norm3(quat[ind, :3]), abs(quat[ind, 3]))

    return np.asarray(angles)


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def approx_equal(quat: double[:, :], other: double[:, :], atol=None, degrees: cython.bint=False):
    if atol is None:
        if degrees:
            warnings.warn("atol must be set to use the degrees flag, "
                            "defaulting to 1e-8 radians.")
        atol = 1e-8  # radians
    elif degrees:
        atol = np.deg2rad(atol)

    q_result = _compose_quat(other, inv(quat))
    angles = magnitude(q_result)
    return angles < atol


@cython.embedsignature(True)
def mean(quat: double[:, :], weights=None):
    if weights is None:
        weights = np.ones(quat.shape[0])
    else:
        weights = np.asarray(weights)
        if weights.ndim != 1:
            raise ValueError("Expected `weights` to be 1 dimensional, got "
                                "shape {}.".format(weights.shape))
        if weights.shape[0] != quat.shape[0]:
            raise ValueError("Expected `weights` to have number of values "
                                "equal to number of rotations, got "
                                "{} values and {} rotations.".format(
                                weights.shape[0], quat.shape[0]))
        if np.any(weights < 0):
            raise ValueError("`weights` must be non-negative.")

    quat = np.asarray(quat)
    K = np.dot(weights * quat.T, quat)
    _, v = np.linalg.eigh(K)
    return v[:, -1]


@cython.embedsignature(True)
def reduce(quat: double[: :], left=None, right=None) -> tuple[double[:, :], int[:] | None, int[:] | None]:
    if left is None and right is None:
        reduced = quat
        return reduced, None, None
    elif right is None:
        right = identity(1)
    elif left is None:
        left = identity(1)

    # Levi-Civita tensor for triple product computations
    e = np.zeros((3, 3, 3))
    e[0, 1, 2] = e[1, 2, 0] = e[2, 0, 1] = 1
    e[0, 2, 1] = e[2, 1, 0] = e[1, 0, 2] = -1

    # We want to calculate the real components of q = l * p * r. It can
    # be shown that:
    #     qs = ls * ps * rs - ls * dot(pv, rv) - ps * dot(lv, rv)
    #          - rs * dot(lv, pv) - dot(cross(lv, pv), rv)
    # where ls and lv denote the scalar and vector components of l.

    def split_rotation(q):
        q = np.atleast_2d(q)
        return q[:, -1], q[:, :-1]

    p = quat
    ps, pv = split_rotation(p)
    ls, lv = split_rotation(left)
    rs, rv = split_rotation(right)

    qs = np.abs(np.einsum('i,j,k', ls, ps, rs) -
                np.einsum('i,jx,kx', ls, pv, rv) -
                np.einsum('ix,j,kx', lv, ps, rv) -
                np.einsum('ix,jx,k', lv, pv, rs) -
                np.einsum('xyz,ix,jy,kz', e, lv, pv, rv))
    qs = np.reshape(np.moveaxis(qs, 1, 0), (qs.shape[1], -1))

    # Find best indices from scalar components
    max_ind = np.argmax(np.reshape(qs, (len(qs), -1)), axis=1)
    left_best = max_ind // len(rv)
    right_best = max_ind % len(rv)

    # Reduce the rotation using the best indices
    left = left[left_best]
    right = right[right_best]
    reduced = _compose_quat(left, _compose_quat(p, right))

    return reduced, left_best, right_best


@cython.embedsignature(True)
def apply(quat: cython.double[:, :], vectors, inverse=False):
    vectors = np.asarray(vectors)
    if vectors.ndim > 2 or vectors.shape[-1] != 3:
        raise ValueError("Expected input of shape (3,) or (P, 3), "
                            "got {}.".format(vectors.shape))

    if vectors.shape == (3,):
        vectors = vectors[None, :]

    matrix = as_matrix(quat)

    n_vectors = vectors.shape[0]
    n_rotations = quat.shape[0]

    if n_vectors != 1 and n_rotations != 1 and n_vectors != n_rotations:
        raise ValueError("Expected equal numbers of rotations and vectors "
                            ", or a single rotation, or a single vector, got "
                            "{} rotations and {} vectors.".format(
                            n_rotations, n_vectors))

    if inverse:
        result = np.einsum('ikj,ik->ij', matrix, vectors)
    else:
        result = np.einsum('ijk,ik->ij', matrix, vectors)

    return result


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def setitem(quat: double[:, :], value: double[:, :], indexer):
    quat[indexer] = value
    return quat


@cython.embedsignature(True)
def align_vectors(a, b, weights=None, return_sensitivity: cython.bint = False):
    # Check input vectors
    a_original = np.array(a, dtype=float)
    b_original = np.array(b, dtype=float)
    a = np.atleast_2d(a_original)
    b = np.atleast_2d(b_original)
    if a.shape[-1] != 3:
        raise ValueError("Expected input `a` to have shape (3,) or "
                            "(N, 3), got {}".format(a_original.shape))
    if b.shape[-1] != 3:
        raise ValueError("Expected input `b` to have shape (3,) or "
                            "(N, 3), got {}".format(b_original.shape))
    if a.shape != b.shape:
        raise ValueError("Expected inputs `a` and `b` to have same shapes"
                            ", got {} and {} respectively.".format(
                                a_original.shape, b_original.shape))
    N = len(a)

    # Check weights
    if weights is None:
        weights = np.ones(N)
    else:
        weights = np.array(weights, dtype=float)
        if weights.ndim != 1:
            raise ValueError("Expected `weights` to be 1 dimensional, got "
                                "shape {}.".format(weights.shape))
        if N > 1 and (weights.shape[0] != N):
            raise ValueError("Expected `weights` to have number of values "
                                "equal to number of input vectors, got "
                                "{} values and {} vectors.".format(
                                weights.shape[0], N))
        if (weights < 0).any():
            raise ValueError("`weights` may not contain negative values")

    # For the special case of a single vector pair, we use the infinite
    # weight code path
    if N == 1:
        weight_is_inf = np.array([True])
    else:
        weight_is_inf = np.isposinf(weights)
    n_inf = np.sum(weight_is_inf)

    # Check for an infinite weight, which indicates that the corresponding
    # vector pair is the primary unmoving reference to which we align the
    # other secondary vectors
    if n_inf > 1:
        raise ValueError("Only one infinite weight is allowed")
    elif n_inf == 1:
        if return_sensitivity:
            raise ValueError("Cannot return sensitivity matrix with an "
                                "infinite weight or one vector pair")
        a_pri = a[weight_is_inf]
        b_pri = b[weight_is_inf]
        a_pri_norm = _norm3(a_pri[0])
        b_pri_norm = _norm3(b_pri[0])
        if a_pri_norm == 0 or b_pri_norm == 0:
            raise ValueError("Cannot align zero length primary vectors")
        a_pri /= a_pri_norm
        b_pri /= b_pri_norm

        # We first find the minimum angle rotation between the primary
        # vectors.
        cross = np.cross(b_pri[0], a_pri[0])
        cross_norm = _norm3(cross)
        theta = atan2(cross_norm, _dot3(a_pri[0], b_pri[0]))
        tolerance = 1e-3  # tolerance for small angle approximation (rad)
        q_flip = np.array([[0.0, 0.0, 0.0, 1.0]])
        if (np.pi - theta) < tolerance:
            # Near pi radians, the Taylor series approximation of x/sin(x)
            # diverges, so for numerical stability we flip pi and then
            # rotate back by the small angle pi - theta
            if cross_norm == 0:
                # For antiparallel vectors, cross = [0, 0, 0] so we need to
                # manually set an arbitrary orthogonal axis of rotation
                i = np.argmin(np.abs(a_pri[0]))
                r = np.zeros(3)
                r[i - 1], r[i - 2] = a_pri[0][i - 2], -a_pri[0][i - 1]
            else:
                r = cross  # Shortest angle orthogonal axis of rotation
            q_flip = from_rotvec((r / np.linalg.norm(r) * np.pi)[None, ...])
            theta = np.pi - theta
            cross = -cross
        if abs(theta) < tolerance:
            # Small angle Taylor series approximation for numerical stability
            theta2 = theta * theta
            r = cross * (1 + theta2 / 6 + theta2 * theta2 * 7 / 360)
        else:
            r = cross * theta / np.sin(theta)
        q_pri = _compose_quat(from_rotvec(r[None, :]), q_flip)

        if N == 1:
            # No secondary vectors, so we are done
            q_opt = q_pri
        else:
            a_sec = a[~weight_is_inf]
            b_sec = b[~weight_is_inf]
            weights_sec = weights[~weight_is_inf]

            # We apply the first rotation to the b vectors to align the
            # primary vectors, resulting in vectors c. The secondary
            # vectors must now be rotated about that primary vector to best
            # align c to a.
            c_sec = apply(q_pri, b_sec)

            # Calculate vector components for the angle calculation. The
            # angle phi to rotate a single 2D vector C to align to 2D
            # vector A in the xy plane can be found with the equation
            # phi = atan2(cross(C, A), dot(C, A))
            #     = atan2(|C|*|A|*sin(phi), |C|*|A|*cos(phi))
            # The below equations perform the same operation, but with the
            # 3D rotation restricted to the 2D plane normal to a_pri, where
            # the secondary vectors are projected into that plane. We then
            # find the composite angle from the weighted sum of the
            # axial components in that plane, minimizing the 2D alignment
            # problem.
            # Note that einsum('ij,ij->i', X, Y) is the row-wise dot
            # product of X and Y.
            sin_term = np.einsum('ij,ij->i', np.cross(c_sec, a_sec), a_pri)
            cos_term = (np.einsum('ij,ij->i', c_sec, a_sec)
                        - (np.einsum('ij,ij->i', c_sec, a_pri)
                            * np.einsum('ij,ij->i', a_sec, a_pri)))
            phi = atan2(np.sum(weights_sec * sin_term),
                        np.sum(weights_sec * cos_term))
            q_sec = from_rotvec((phi * a_pri[0])[None, :])

            # Compose these to get the optimal rotation
            q_opt = _compose_quat(q_sec, q_pri)

        # Calculated the root sum squared distance. We force the error to
        # be zero for the infinite weight vectors since they will align
        # exactly.
        weights_inf_zero = weights.copy()
        if N > 1 or np.isposinf(weights[0]):
            # Skip non-infinite weight single vectors pairs, we used the
            # infinite weight code path but don't want to zero that weight
            weights_inf_zero[weight_is_inf] = 0
        a_est = apply(q_opt, b)
        rssd = np.sqrt(np.sum(weights_inf_zero @ (a - a_est)**2))

        return q_opt[0, ...], rssd, None

    # If no infinite weights and multiple vectors, proceed with normal
    # algorithm
    # Note that einsum('ji,jk->ik', X, Y) is equivalent to np.dot(X.T, Y)
    B = np.einsum('ji,jk->ik', weights[:, None] * a, b)
    u, s, vh = np.linalg.svd(B)

    # Correct improper rotation if necessary (as in Kabsch algorithm)
    if np.linalg.det(u @ vh) < 0:
        s[-1] = -s[-1]
        u[:, -1] = -u[:, -1]

    C = np.dot(u, vh)

    if s[1] + s[2] < 1e-16 * s[0]:
        warnings.warn("Optimal rotation is not uniquely or poorly defined "
                        "for the given sets of vectors.")

    rssd = np.sqrt(max(
        np.sum(weights * np.sum(b ** 2 + a ** 2, axis=1)) - 2 * np.sum(s),
        0))

    if return_sensitivity:
        zeta = (s[0] + s[1]) * (s[1] + s[2]) * (s[2] + s[0])
        kappa = s[0] * s[1] + s[1] * s[2] + s[2] * s[0]
        with np.errstate(divide='ignore', invalid='ignore'):
            sensitivity = np.mean(weights) / zeta * (
                    kappa * np.eye(3) + np.dot(B, B.T))
        return from_matrix(C), rssd, sensitivity
    return from_matrix(C), rssd, None


@cython.embedsignature(True)
def pow(quat: double[:, :], n: float) -> double[:, :]:
    # Exact short-cuts
    if n == 0:
        return identity(quat.shape[0])
    elif n == -1:
        return inv(quat)
    elif n == 1:
        return quat
    # general scaling of rotation angle
    return from_rotvec(n * as_rotvec(quat))


@cython.embedsignature(True)
def from_davenport(axes, order: str, angles: double[:, :], degrees: cython.bint = False) -> double[:, :]:
    """Initialize from Davenport angles.

    Rotations in 3-D can be represented by a sequence of 3
    rotations around a sequence of axes.

    The three rotations can either be in a global frame of reference
    (extrinsic) or in a body centred frame of reference (intrinsic), which
    is attached to, and moves with, the object under rotation [1]_.

    For both Euler angles and Davenport angles, consecutive axes must
    be are orthogonal (``axis2`` is orthogonal to both ``axis1`` and
    ``axis3``). For Euler angles, there is an additional relationship
    between ``axis1`` or ``axis3``, with two possibilities:

        - ``axis1`` and ``axis3`` are also orthogonal (asymmetric sequence)
        - ``axis1 == axis3`` (symmetric sequence)

    For Davenport angles, this last relationship is relaxed [2]_, and only
    the consecutive orthogonal axes requirement is maintained.

    Parameters
    ----------
    axes : array_like, shape (3,) or ([1 or 2 or 3], 3)
        Axis of rotation, if one dimensional. If two dimensional, describes the
        sequence of axes for rotations, where each axes[i, :] is the ith
        axis. If more than one axis is given, then the second axis must be
        orthogonal to both the first and third axes.
    order : string
        If it is equal to 'e' or 'extrinsic', the sequence will be
        extrinsic. If it is equal to 'i' or 'intrinsic', sequence
        will be treated as intrinsic.
    angles : float or array_like, shape (N,) or (N, [1 or 2 or 3])
        Euler angles specified in radians (`degrees` is False) or degrees
        (`degrees` is True).
        For a single axis, `angles` can be:

        - a single value
        - array_like with shape (N,), where each `angle[i]`
            corresponds to a single rotation
        - array_like with shape (N, 1), where each `angle[i, 0]`
            corresponds to a single rotation

        For 2 and 3 axes, `angles` can be:

        - array_like with shape (W,) where `W` is the number of rows of
            `axes`, which corresponds to a single rotation with `W` axes
        - array_like with shape (N, W) where each `angle[i]`
            corresponds to a sequence of Davenport angles describing a
            single rotation

    degrees : bool, optional
        If True, then the given angles are assumed to be in degrees.
        Default is False.

    Returns
    -------
    rotation : `Rotation` instance
        Object containing the rotation represented by the sequence of
        rotations around given axes with given angles.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Euler_angles#Definition_by_intrinsic_rotations
    .. [2] Shuster, Malcolm & Markley, Landis. (2003). Generalization of
            the Euler Angles. Journal of the Astronautical Sciences. 51. 123-132. 10.1007/BF03546304.

    Examples
    --------
    >>> from scipy.spatial.transform import Rotation as R

    Davenport angles are a generalization of Euler angles, when we use the
    canonical basis axes:

    >>> ex = [1, 0, 0]
    >>> ey = [0, 1, 0]
    >>> ez = [0, 0, 1]

    Initialize a single rotation with a given axis sequence:

    >>> axes = [ez, ey, ex]
    >>> r = R.from_davenport(axes, 'extrinsic', [90, 0, 0], degrees=True)
    >>> r.as_quat().shape
    (4,)

    It is equivalent to Euler angles in this case:

    >>> r.as_euler('zyx', degrees=True)
    array([90.,  0., -0.])

    Initialize multiple rotations in one object:

    >>> r = R.from_davenport(axes, 'extrinsic', [[90, 45, 30], [35, 45, 90]], degrees=True)
    >>> r.as_quat().shape
    (2, 4)

    Using only one or two axes is also possible:

    >>> r = R.from_davenport([ez, ex], 'extrinsic', [[90, 45], [35, 45]], degrees=True)
    >>> r.as_quat().shape
    (2, 4)

    Non-canonical axes are possible, and they do not need to be normalized,
    as long as consecutive axes are orthogonal:

    >>> e1 = [2, 0, 0]
    >>> e2 = [0, 1, 0]
    >>> e3 = [1, 0, 1]
    >>> axes = [e1, e2, e3]
    >>> r = R.from_davenport(axes, 'extrinsic', [90, 45, 30], degrees=True)
    >>> r.as_quat()
    [ 0.701057,  0.430459, -0.092296,  0.560986]
    """
    if order in ['e', 'extrinsic']:
        extrinsic = True
    elif order in ['i', 'intrinsic']:
        extrinsic = False
    else:
        raise ValueError("order should be 'e'/'extrinsic' for extrinsic "
                            "sequences or 'i'/'intrinsic' for intrinsic "
                            "sequences, got {}".format(order))

    axes = np.asarray(axes)
    if axes.ndim == 1:
        axes = axes.reshape([1, 3])

    num_axes = len(axes)

    if axes.shape[1] != 3:
        raise ValueError("Axes must be vectors of length 3.")

    if num_axes < 1 or num_axes > 3:
        raise ValueError("Expected up to 3 axes, got {}".format(num_axes))

    # normalize axes
    norm = np.repeat(np.linalg.norm(axes, axis=1), 3)
    axes = axes / norm.reshape(num_axes, 3)

    if (num_axes > 1 and abs(np.dot(axes[0], axes[1])) >= 1e-7 or
        num_axes > 2 and abs(np.dot(axes[1], axes[2])) >= 1e-7):
        raise ValueError("Consecutive axes must be orthogonal.")

    angles, is_single = _format_angles(angles, degrees, num_axes)

    q = identity(len(angles))
    for i in range(num_axes):
        qi = from_rotvec(angles[:, i, np.newaxis] * axes[i])
        if extrinsic:
            q = compose_quat(qi, q)
        else:
            q = compose_quat(q, qi)

    return q[0] if is_single else q


cdef class Rotation:
    cdef double[:, :] _quat
    cdef bint _single

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __init__(self, quat, normalize=True, copy=True, scalar_first=False):
        self._single = False
        quat = np.asarray(quat, dtype=float)

        if (quat.ndim not in [1, 2]
            or quat.shape[len(quat.shape) - 1] != 4
            or quat.shape[0] == 0):
            raise ValueError("Expected `quat` to have shape (4,) or (N, 4), "
                             f"got {quat.shape}.")

        # If a single quaternion is given, convert it to a 2D 1 x 4 matrix but
        # set self._single to True so that we can return appropriate objects
        # in the `to_...` methods
        if quat.shape == (4,):
            quat = quat[None, :]
            self._single = True

        cdef Py_ssize_t num_rotations = quat.shape[0]

        if scalar_first:
            quat = np.roll(quat, -1, axis=1)
        elif normalize or copy:
            quat = quat.copy()

        if normalize:
            for ind in range(num_rotations):
                if isnan(_normalize4(quat[ind, :])):
                    raise ValueError("Found zero norm quaternions in `quat`.")

        self._quat = quat

    @cython.embedsignature(True)
    @classmethod
    def create_group(cls, group, axis='Z'):
        """Create a 3D rotation group.

        Parameters
        ----------
        group : string
            The name of the group. Must be one of 'I', 'O', 'T', 'Dn', 'Cn',
            where `n` is a positive integer. The groups are:

                * I: Icosahedral group
                * O: Octahedral group
                * T: Tetrahedral group
                * D: Dicyclic group
                * C: Cyclic group

        axis : integer
            The cyclic rotation axis. Must be one of ['X', 'Y', 'Z'] (or
            lowercase). Default is 'Z'. Ignored for groups 'I', 'O', and 'T'.

        Returns
        -------
        rotation : `Rotation` instance
            Object containing the elements of the rotation group.

        Notes
        -----
        This method generates rotation groups only. The full 3-dimensional
        point groups [PointGroups]_ also contain reflections.

        References
        ----------
        .. [PointGroups] `Point groups
           <https://en.wikipedia.org/wiki/Point_groups_in_three_dimensions>`_
           on Wikipedia.
        """
        return create_group(cls, group, axis=axis)


class Slerp:
    """Spherical Linear Interpolation of Rotations.

    The interpolation between consecutive rotations is performed as a rotation
    around a fixed axis with a constant angular velocity [1]_. This ensures
    that the interpolated rotations follow the shortest path between initial
    and final orientations.

    Parameters
    ----------
    times : array_like, shape (N,)
        Times of the known rotations. At least 2 times must be specified.
    rotations : `Rotation` instance
        Rotations to perform the interpolation between. Must contain N
        rotations.

    Methods
    -------
    __call__

    See Also
    --------
    Rotation

    Notes
    -----
    .. versionadded:: 1.2.0

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Slerp#Quaternion_Slerp

    Examples
    --------
    >>> from scipy.spatial.transform import Rotation as R
    >>> from scipy.spatial.transform import Slerp

    Setup the fixed keyframe rotations and times:

    >>> key_rots = R.random(5, random_state=2342345)
    >>> key_times = [0, 1, 2, 3, 4]

    Create the interpolator object:

    >>> slerp = Slerp(key_times, key_rots)

    Interpolate the rotations at the given times:

    >>> times = [0, 0.5, 0.25, 1, 1.5, 2, 2.75, 3, 3.25, 3.60, 4]
    >>> interp_rots = slerp(times)

    The keyframe rotations expressed as Euler angles:

    >>> key_rots.as_euler('xyz', degrees=True)
    array([[ 14.31443779, -27.50095894,  -3.7275787 ],
           [ -1.79924227, -24.69421529, 164.57701743],
           [146.15020772,  43.22849451, -31.34891088],
           [ 46.39959442,  11.62126073, -45.99719267],
           [-88.94647804, -49.64400082, -65.80546984]])

    The interpolated rotations expressed as Euler angles. These agree with the
    keyframe rotations at both endpoints of the range of keyframe times.

    >>> interp_rots.as_euler('xyz', degrees=True)
    array([[  14.31443779,  -27.50095894,   -3.7275787 ],
           [   4.74588574,  -32.44683966,   81.25139984],
           [  10.71094749,  -31.56690154,   38.06896408],
           [  -1.79924227,  -24.69421529,  164.57701743],
           [  11.72796022,   51.64207311, -171.7374683 ],
           [ 146.15020772,   43.22849451,  -31.34891088],
           [  68.10921869,   20.67625074,  -48.74886034],
           [  46.39959442,   11.62126073,  -45.99719267],
           [  12.35552615,    4.21525086,  -64.89288124],
           [ -30.08117143,  -19.90769513,  -78.98121326],
           [ -88.94647804,  -49.64400082,  -65.80546984]])

    """
    def __init__(self, times, rotations):
        if not isinstance(rotations, Rotation):
            raise TypeError("`rotations` must be a `Rotation` instance.")

        if rotations.single or len(rotations) == 1:
            raise ValueError("`rotations` must be a sequence of at least 2 rotations.")

        times = np.asarray(times)
        if times.ndim != 1:
            raise ValueError("Expected times to be specified in a 1 "
                             "dimensional array, got {} "
                             "dimensions.".format(times.ndim))

        if times.shape[0] != len(rotations):
            raise ValueError("Expected number of rotations to be equal to "
                             "number of timestamps given, got {} rotations "
                             "and {} timestamps.".format(
                                len(rotations), times.shape[0]))
        self.times = times
        self.timedelta = np.diff(times)

        if np.any(self.timedelta <= 0):
            raise ValueError("Times must be in strictly increasing order.")

        self.rotations = rotations[:-1]
        self.rotvecs = (self.rotations.inv() * rotations[1:]).as_rotvec()

    def __call__(self, times):
        """Interpolate rotations.

        Compute the interpolated rotations at the given `times`.

        Parameters
        ----------
        times : array_like
            Times to compute the interpolations at. Can be a scalar or
            1-dimensional.

        Returns
        -------
        interpolated_rotation : `Rotation` instance
            Object containing the rotations computed at given `times`.

        """
        # Clearly differentiate from self.times property
        compute_times = np.asarray(times)
        if compute_times.ndim > 1:
            raise ValueError("`times` must be at most 1-dimensional.")

        single_time = compute_times.ndim == 0
        compute_times = np.atleast_1d(compute_times)

        # side = 'left' (default) excludes t_min.
        ind = np.searchsorted(self.times, compute_times) - 1
        # Include t_min. Without this step, index for t_min equals -1
        ind[compute_times == self.times[0]] = 0
        if np.any(np.logical_or(ind < 0, ind > len(self.rotations) - 1)):
            raise ValueError("Interpolation times must be within the range "
                             "[{}, {}], both inclusive.".format(
                                self.times[0], self.times[-1]))

        alpha = (compute_times - self.times[ind]) / self.timedelta[ind]

        result = (self.rotations[ind] *
                  Rotation.from_rotvec(self.rotvecs[ind] * alpha[:, None]))

        if single_time:
            result = result[0]

        return result
