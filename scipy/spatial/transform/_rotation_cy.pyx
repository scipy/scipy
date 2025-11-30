# cython: cpow=True

import re
import warnings
import numpy as np
from scipy._lib._util import check_random_state, _transition_to_rng

cimport numpy as np
cimport cython
from cython.view cimport array
from libc.math cimport sqrt, sin, cos, atan2, acos, hypot, isnan, NAN, pi

np.import_array()

# utilities for empty array initialization
cdef inline double[:] _empty1(int n) noexcept:
    if n == 0:
        return array(shape=(1,), itemsize=sizeof(double), format=b"d")[:0]
    return array(shape=(n,), itemsize=sizeof(double), format=b"d")
cdef inline double[:, :] _empty2(int n1, int n2) noexcept :
    if n1 == 0:
        return array(shape=(1, n2), itemsize=sizeof(double), format=b"d")[:0]
    return array(shape=(n1, n2), itemsize=sizeof(double), format=b"d")
cdef inline double[:, :, :] _empty3(int n1, int n2, int n3) noexcept:
    if n1 == 0:
        return array(shape=(1, n2, n3), itemsize=sizeof(double), format=b"d")[:0]
    return array(shape=(n1, n2, n3), itemsize=sizeof(double), format=b"d")
cdef inline double[:, :] _zeros2(int n1, int n2) noexcept:
    if n1 == 0:
        return array(shape=(1, n2), itemsize=sizeof(double), format=b"d")[:0]
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
cdef inline int _get_angles(
    double[:] angles, bint extrinsic, bint symmetric, bint sign,
    double lamb, double a, double b, double c, double d) noexcept:
    # Returns 1 if a gimbal warning is detected and 0 otherwise (so that
    # warnings can be handled at a higher level)

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
        return 1
    return 0


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:, :] _compute_euler_from_quat(
    np.ndarray[double, ndim=2] quat,
    const uchar[:] seq,
    bint extrinsic,
    bint suppress_warnings
) noexcept:
    # The algorithm assumes extrinsic frame transformations. The algorithm
    # in the paper is formulated for rotation quaternions, which are stored
    # directly by Rotation.
    # Adapt the algorithm for our case by reversing both axis sequence and
    # angles for intrinsic rotations when needed

    if not extrinsic:
        seq = seq[::-1]

    cdef int i = _elementary_basis_index(seq[0])
    cdef int j = _elementary_basis_index(seq[1])
    cdef int k = _elementary_basis_index(seq[2])

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
    cdef int n_gimbal_warnings = 0

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

        n_gimbal_warnings += _get_angles(
            angles[ind], extrinsic, symmetric, sign, pi / 2, a, b, c, d
        )

    if n_gimbal_warnings > 0 and not suppress_warnings:
        warnings.warn("Gimbal lock detected. Setting third angle to zero "
                      "since it is not possible to uniquely determine "
                      "all angles.", stacklevel=2)

    return angles

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:, :] _compute_davenport_from_quat(
    np.ndarray[double, ndim=2] quat, np.ndarray[double, ndim=1] n1,
    np.ndarray[double, ndim=1] n2, np.ndarray[double, ndim=1] n3,
    bint extrinsic,
    bint suppress_warnings
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
    cdef int n_gimbal_warnings = 0

    for ind in range(num_rotations):
        _compose_quat_single(quat_lamb, quat[ind], quat_transformed)

        # Step 1
        # Permutate quaternion elements
        a = quat_transformed[3]
        b = _dot3(quat_transformed[:3], n1)
        c = _dot3(quat_transformed[:3], n2)
        d = _dot3(quat_transformed[:3], n_cross)

        n_gimbal_warnings += _get_angles(angles[ind], extrinsic, False, 1, lamb, a, b, c, d)

        if correct_set:
            angles[ind, 1] = -angles[ind, 1]

    if n_gimbal_warnings > 0 and not suppress_warnings:
        warnings.warn("Gimbal lock detected. Setting third angle to zero "
                      "since it is not possible to uniquely determine "
                      "all angles.", stacklevel=2)

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
cdef inline double[:, :] _compose_quat(
    const double[:, :] p, const double[:, :] q
) noexcept:
    cdef Py_ssize_t n = q.shape[0] if p.shape[0] == 1 else p.shape[0]
         
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

def compose_quat(p, q):
    """Composition of quaternions."""
    return _compose_quat(p, q)

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


@cython.boundscheck(False)
@cython.wraparound(False)
def from_quat(double[:, :] quat, bint normalize=True, bint copy=True, bint scalar_first=False):
    if quat.ndim != 2 or quat.shape[1] != 4:
        raise ValueError(f"Expected `quat` to have shape (N, 4), got {quat.shape}.")

    cdef Py_ssize_t num_rotations = quat.shape[0]

    if num_rotations > 0:  # Avoid 0-sized axis errors
        if scalar_first:
            quat = np.roll(quat, -1, axis=1)
        elif normalize or copy:
            quat = quat.copy()

        if normalize:
            for ind in range(num_rotations):
                if isnan(_normalize4(quat[ind, :])):
                    raise ValueError("Found zero norm quaternions in `quat`.")

    return np.asarray(quat, dtype=float)


@cython.boundscheck(False)
@cython.wraparound(False)
def from_euler(seq, angles, bint degrees=False):
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
    
    angles = np.asarray(angles, dtype=float)

    if angles.ndim > 2:  # The backend should never be called with these inputs
        raise ValueError("Cython backend is only compatible with up to 2D inputs")

    if degrees:
        angles = np.deg2rad(angles)

    is_single = angles.ndim < 2
    angles = np.atleast_2d(angles)  # Ensure 0, 1 and 2D arrays are 2D
    
    if angles.shape[1] != num_axes:
        raise ValueError("Expected last dimension of `angles` to match number of"
                         f" sequence axes specified, got {angles.shape[1]}.")

    quat = _elementary_quat_compose(seq.encode(), angles, intrinsic)

    if is_single:
        return np.asarray(quat, dtype=float)[0]
    return np.asarray(quat, dtype=float)


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def from_matrix(matrix):
    cdef int ind
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
    if is_single:
        return _from_matrix_orthogonal(matrix[0])
    return _from_matrix_orthogonal(matrix)


@cython.boundscheck(False)
@cython.wraparound(False)
def _from_matrix_orthogonal(matrix):
    is_single = False
    if matrix.shape == (3, 3):
        matrix = matrix[np.newaxis, :, :]
        is_single = True

    cdef double[:, :, :] cmatrix
    cmatrix = matrix
    cdef Py_ssize_t num_rotations = cmatrix.shape[0]
    cdef Py_ssize_t i, j, k
    cdef double[:] decision = _empty1(4)
    cdef int choice
    cdef int ind

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
        return np.asarray(quat, dtype=float)[0]
    return np.asarray(quat, dtype=float)


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def from_rotvec(rotvec, bint degrees=False):
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
        return np.asarray(quat, dtype=float)[0]
    return np.asarray(quat, dtype=float)


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def from_mrp(mrp):
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
    

@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def as_quat(double[:, :] quat, bint canonical=False, *, bint scalar_first=False) -> double[:, :]:
    q = np.array(quat, copy=True)
    if canonical:
        _quat_canonical(q)

    if scalar_first:
        q = np.roll(q, 1, axis=1)

    return q


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def as_matrix(double[:, :] quat) -> double[:, :, :]:

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
def as_rotvec(double[:, :] quat, bint degrees=False) -> double[:, :]:
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
@cython.boundscheck(False)
@cython.wraparound(False)
def as_mrp(double[:, :] quat) -> double[:, :]:
    cdef Py_ssize_t num_rotations = len(quat)
    cdef double[:, :] mrps = _empty2(num_rotations, 3)
    cdef int sign
    cdef double denominator

    # TODO: Check if broadcasting is possible similar to the xp backend
    for ind in range(num_rotations):

        # Ensure we are calculating the set of MRPs that correspond
        # to a rotation of <= 180
        sign = -1 if quat[ind, 3] < 0 else 1

        denominator = 1 + sign * quat[ind, 3]
        for i in range(3):
            mrps[ind, i] = sign * quat[ind, i] / denominator

    return np.asarray(mrps)


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def as_euler(
    double[:, :] quat,
    seq,
    bint degrees=False,
    bint suppress_warnings=False
) -> double[:, :]:
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
    # directly by Rotation. See:
    # Bernardes E, Viollet S (2022) Quaternion to Euler angles conversion: A
    # direct, general and computationally efficient method.
    # PLoS ONE 17(11): e0276302. https://doi.org/10.1371/journal.pone.0276302
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
    cdef int n_gimbal_warnings = 0

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

        n_gimbal_warnings += _get_angles(
            angles[ind], extrinsic, symmetric, sign, pi / 2, a, b, c, d
        )

    if n_gimbal_warnings > 0 and not suppress_warnings:
        warnings.warn("Gimbal lock detected. Setting third angle to zero "
                      "since it is not possible to uniquely determine "
                      "all angles.", stacklevel=2)

    if degrees:
        angles = np.rad2deg(angles)
    return np.asarray(angles)


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def as_davenport(
    double[:, :] quat,
    double[:, :] axes,
    order,
    bint degrees=False,
    bint suppress_warnings=False
) -> double[:, :]:
    cdef np.ndarray[double, ndim=2] q = np.asarray(quat)

    cdef bint extrinsic
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

    if axes.ndim != 2 or axes.shape[1] != 3:
        raise ValueError("Axes must be vectors of length 3.")

    cdef np.ndarray[double, ndim=1] n1, n2, n3
    n1, n2, n3 = np.asarray(axes)

    # normalize axes
    n1 = n1 / np.linalg.norm(n1)
    n2 = n2 / np.linalg.norm(n2)
    n3 = n3 / np.linalg.norm(n3)

    if np.dot(n1, n2) >= 1e-7:
        raise ValueError("Consecutive axes must be orthogonal.")
    if np.dot(n2, n3) >= 1e-7:
        raise ValueError("Consecutive axes must be orthogonal.")

    angles = np.asarray(_compute_davenport_from_quat(
        q, n1, n2, n3, extrinsic, suppress_warnings
    ))

    if degrees:
        angles = np.rad2deg(angles)

    return angles


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def inv(double[:, :] quat) -> double[:, :]:
    cdef np.ndarray[double, ndim=2] q_inv = np.array(quat, copy=True)
    q_inv[:, 0] *= -1
    q_inv[:, 1] *= -1
    q_inv[:, 2] *= -1
    return q_inv


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@_transition_to_rng('random_state', position_num=1)
def random(num: int | None = None, rng=None, shape: int | tuple[int, ...] | None = None):
    rng = check_random_state(rng)
    if num is None and shape is None:
        shape = ()
    elif num is not None:
        shape = (num,)
    elif isinstance(shape, int):
        shape = (shape,)
    elif not isinstance(shape, tuple):
        raise ValueError("`shape` must be an int or a tuple of ints or None.")
    return rng.normal(size=shape + (4,))


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def identity(num: int | None, shape: int | tuple[int, ...] | None = None):
    if num is None and shape is None:
        shape = ()
    elif num is not None:
        shape = (num,)
    elif isinstance(shape, int):
        shape = (shape,)
    elif not isinstance(shape, tuple):
        raise ValueError("`shape` must be an int or a tuple of ints or None.")
    q = np.zeros(shape + (4,), dtype=np.float64)
    q[..., 3] = 1
    return q


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def magnitude(double[:, :] quat) -> double[:, :]:
    cdef Py_ssize_t num_rotations = quat.shape[0]
    cdef double[:] angles = _empty1(num_rotations)

    for ind in range(num_rotations):
        angles[ind] = 2 * atan2(_norm3(quat[ind, :3]), abs(quat[ind, 3]))

    return np.asarray(angles)


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def approx_equal(double[:, :] quat, double[:, :] other, atol = None, bint degrees=False) -> bool[:]:
    len_self = len(quat)
    len_other = len(other)
    if not(len_self == 1 or len_other == 1 or len_self == len_other):
        raise ValueError(
            f"Expected broadcastable shapes in both rotations, got {len_self} "
            f"rotations in first and {len_other} rotations in second object."
        )

    if atol is None:
        if degrees:
            warnings.warn("atol must be set to use the degrees flag, "
                            "defaulting to 1e-8 radians.")
        atol = 1e-8  # radians
    elif degrees:
        atol = np.deg2rad(atol)

    angles = magnitude(compose_quat(other, inv(quat)))
    return angles < atol


@cython.embedsignature(True)
@cython.boundscheck(False)
def mean(double[:, :] quat, weights=None, axis=None):
    if quat.shape[0] == 0:
        raise ValueError("Mean of an empty rotation set is undefined.")
    # The Cython path assumes quat is Nx4, so axis has to be None, 0, -1, (0,), (-1,),
    # or (). The code path is unchanged for any of the options except (), where we
    # immediately return the quaternion
    if axis == ():
        return quat

    if axis is None:
        axis = (0,)
    if isinstance(axis, int):
        axis = (axis,)
    if not isinstance(axis, tuple):  # Must be tuple by now
        raise ValueError("`axis` must be None, int, or tuple of ints.")
    if min(axis) < -1 or max(axis) > 0:
        raise ValueError(
            f"axis {axis} is out of bounds for rotation with shape "
            f"{np.asarray(quat).shape[:-1]}."
        )
    # Axis must be 0 for the cython backend. Everything else should have raised an
    # error during validation.
    axis = (0,)

    if weights is None:
        weights = np.ones(quat.shape[0])
    else:
        weights = np.asarray(weights)
        if np.any(weights < 0):
            raise ValueError("`weights` must be non-negative.")
        if weights.ndim != 1:
            raise ValueError(f"Expected `weights` to be 1 dimensional, got "
                             f"{weights.shape}.")
        if weights.shape[0] != quat.shape[0]:
            raise ValueError("Expected `weights` to have number of values equal to "
                             f"number of rotations, got {weights.shape[0]} values and "
                             f"{quat.shape[0]} rotations.")

    quat = np.asarray(quat)
    K = np.dot(weights * quat.T, quat)
    _, v = np.linalg.eigh(K)
    return v[:, -1]


@cython.embedsignature(True)
@cython.boundscheck(False)
def reduce(double[:, :] quat, left=None, right=None):
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
@cython.boundscheck(False)
def apply(double[:, :] quat, double[:, :] vectors, bint inverse=False) -> double[:, :]:
    cdef Py_ssize_t n_vectors = len(vectors)
    cdef Py_ssize_t n_rotations = len(quat)

    if n_vectors != 1 and n_rotations != 1 and n_vectors != n_rotations:
        raise ValueError(
            f"Cannot broadcast {n_rotations} rotations to {n_vectors} vectors."
        )

    cdef np.ndarray matrix = as_matrix(quat)

    if inverse:
        matrix = np.swapaxes(matrix, -1, -2)

    if n_vectors == 1 and n_rotations == 1:
        return np.matmul(matrix, vectors[0])

    if n_rotations == 1:
        # Single rotation/many vectors, use matmul for speed: The axes argument
        # is such that the input arguments don't need to be transposed and the
        # output argument is continuous in memory.
        return np.matmul(matrix, vectors, axes=[(-2, -1), (-1, -2), (-1, -2)])[0]

    # for stacks of matrices einsum is faster
    return np.einsum('ijk,ik->ij', matrix, vectors)


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def setitem(quat: double[:, :], value: double[:, :], indexer):
    quat[indexer] = value
    return quat


@cython.embedsignature(True)
def align_vectors(a, b, weights=None, bint return_sensitivity=False):
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
    elif n_inf == 1:  # TODO: Refactor into _align_vectors_fixed as in the xp backend
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
        return _from_matrix_orthogonal(C), rssd, sensitivity
    return _from_matrix_orthogonal(C), rssd, None


@cython.embedsignature(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def pow(double[:, :] quat, n) -> double[:, :]:
    if isinstance(n, np.ndarray) and n.ndim != 0 and n.shape != (1,):
        raise ValueError("Array exponent must be a scalar")
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
@cython.boundscheck(False)
@cython.wraparound(False)
def from_davenport(axes, order: str, angles, bint degrees=False):
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

    angles = np.asarray(angles, dtype=float)
    if degrees:
        angles = np.deg2rad(angles)

    is_single = angles.ndim < 2
    if angles.ndim != 2:
        angles = np.atleast_2d(angles)

    if angles.shape[1] != axes.shape[0]:
        raise ValueError("Expected `angles` to match number of axes, got "
                         f"{angles.shape[1]} angles and {axes.shape[0]} axes.")

    q = identity(len(angles))
    for i in range(num_axes):
        qi = from_rotvec(angles[:, i, np.newaxis] * axes[i])
        if extrinsic:
            q = compose_quat(qi, q)
        else:
            q = compose_quat(q, qi)

    if is_single:
        return np.asarray(q, dtype=float)[0]
    return np.asarray(q, dtype=float)
