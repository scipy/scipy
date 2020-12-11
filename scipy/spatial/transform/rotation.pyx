import re
import warnings
import numpy as np
from scipy._lib._util import check_random_state
from ._rotation_groups import create_group

cimport numpy as np
cimport cython
from cython.view cimport array
from libc.math cimport sqrt, sin, cos, atan2, acos
from numpy.math cimport PI as pi, NAN, isnan # avoid MSVC error

# utilities for empty array initialization
cdef inline double[:] _empty1(int n):
    return array(shape=(n,), itemsize=sizeof(double), format=b"d")
cdef inline double[:, :] _empty2(int n1, int n2):
    return array(shape=(n1, n2), itemsize=sizeof(double), format=b"d")
cdef inline double[:, :, :] _empty3(int n1, int n2, int n3):
    return array(shape=(n1, n2, n3), itemsize=sizeof(double), format=b"d")
cdef inline double[:, :] _zeros2(int n1, int n2):
    cdef double[:, :] arr = array(shape=(n1, n2),
        itemsize=sizeof(double), format=b"d")
    arr[:, :] = 0
    return arr

# flat implementations of numpy functions
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline double[:] _cross3(const double[:] a, const double[:] b):
    cdef double[:] result = _empty1(3)
    result[0] = a[1]*b[2] - a[2]*b[1]
    result[1] = a[2]*b[0] - a[0]*b[2]
    result[2] = a[0]*b[1] - a[1]*b[0]
    return result

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline double _dot3(const double[:] a, const double[:] b) nogil:
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline double _norm3(const double[:] elems) nogil:
    return sqrt(_dot3(elems, elems))

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline double _normalize4(double[:] elems) nogil:
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
cdef inline int _argmax4(double[:] a) nogil:
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
cdef inline const double[:] _elementary_basis_vector(uchar axis):
    if axis == b'x': return _ex
    elif axis == b'y': return _ey
    elif axis == b'z': return _ez

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:, :] _compute_euler_from_matrix(
    np.ndarray[double, ndim=3] matrix, const uchar[:] seq, bint extrinsic=False
):
    # The algorithm assumes intrinsic frame transformations. The algorithm
    # in the paper is formulated for rotation matrices which are transposition
    # rotation matrices used within Rotation.
    # Adapt the algorithm for our case by
    # 1. Instead of transposing our representation, use the transpose of the
    #    O matrix as defined in the paper, and be careful to swap indices
    # 2. Reversing both axis sequence and angles for extrinsic rotations

    if extrinsic:
        seq = seq[::-1]

    cdef Py_ssize_t num_rotations = matrix.shape[0]

    # Step 0
    # Algorithm assumes axes as column vectors, here we use 1D vectors
    cdef const double[:] n1 = _elementary_basis_vector(seq[0])
    cdef const double[:] n2 = _elementary_basis_vector(seq[1])
    cdef const double[:] n3 = _elementary_basis_vector(seq[2])

    # Step 2
    cdef double sl = _dot3(_cross3(n1, n2), n3)
    cdef double cl = _dot3(n1, n3)

    # angle offset is lambda from the paper referenced in [2] from docstring of
    # `as_euler` function
    cdef double offset = atan2(sl, cl)
    cdef double[:, :] c_ = _empty2(3, 3)
    c_[0, :] = n2
    c_[1, :] = _cross3(n1, n2)
    c_[2, :] = n1
    cdef np.ndarray[double, ndim=2] c = np.asarray(c_)

    rot = np.array([
        [1, 0, 0],
        [0, cl, sl],
        [0, -sl, cl],
    ])

    # some forward definitions
    cdef double[:, :] angles = _empty2(num_rotations, 3)
    cdef double[:, :] matrix_trans # transformed matrix
    cdef double[:] _angles # accessor for each rotation
    cdef np.ndarray[double, ndim=2] res
    cdef double eps = 1e-7
    cdef bint safe1, safe2, safe, adjust

    for ind in range(num_rotations):
        _angles = angles[ind, :]

        # Step 3
        res = np.dot(c, matrix[ind, :, :])
        matrix_trans = np.dot(res, c.T.dot(rot))

        # Step 4
        # Ensure less than unit norm
        matrix_trans[2, 2] = min(matrix_trans[2, 2], 1)
        matrix_trans[2, 2] = max(matrix_trans[2, 2], -1)
        _angles[1] = acos(matrix_trans[2, 2])

        # Steps 5, 6
        safe1 = abs(_angles[1]) >= eps
        safe2 = abs(_angles[1] - <double>pi) >= eps
        safe = safe1 and safe2

        # Step 4 (Completion)
        _angles[1] += offset

        # 5b
        if safe:
            _angles[0] = atan2(matrix_trans[0, 2], -matrix_trans[1, 2])
            _angles[2] = atan2(matrix_trans[2, 0], matrix_trans[2, 1])

        if extrinsic:
            # For extrinsic, set first angle to zero so that after reversal we
            # ensure that third angle is zero
            # 6a
            if not safe:
                _angles[0] = 0
            # 6b
            if not safe1:
                _angles[2] = atan2(matrix_trans[1, 0] - matrix_trans[0, 1],
                                   matrix_trans[0, 0] + matrix_trans[1, 1])
            # 6c
            if not safe2:
                _angles[2] = -atan2(matrix_trans[1, 0] + matrix_trans[0, 1],
                                    matrix_trans[0, 0] - matrix_trans[1, 1])
        else:
            # For instrinsic, set third angle to zero
            # 6a
            if not safe:
                _angles[2] = 0
            # 6b
            if not safe1:
                _angles[0] = atan2(matrix_trans[1, 0] - matrix_trans[0, 1],
                                   matrix_trans[0, 0] + matrix_trans[1, 1])
            # 6c
            if not safe2:
                _angles[0] = atan2(matrix_trans[1, 0] + matrix_trans[0, 1],
                                   matrix_trans[0, 0] - matrix_trans[1, 1])

        # Step 7
        if seq[0] == seq[2]:
            # lambda = 0, so we can only ensure angle2 -> [0, pi]
            adjust = _angles[1] < 0 or _angles[1] > pi
        else:
            # lambda = + or - pi/2, so we can ensure angle2 -> [-pi/2, pi/2]
            adjust = _angles[1] < -pi / 2 or _angles[1] > pi / 2

        # Dont adjust gimbal locked angle sequences
        if adjust and safe:
            _angles[0] += pi
            _angles[1] = 2 * offset - _angles[1]
            _angles[2] -= pi

        for i in range(3):
            if _angles[i] < -pi:
                _angles[i] += 2 * pi
            elif _angles[i] > pi:
                _angles[i] -= 2 * pi

        # Step 8
        if not safe:
            warnings.warn("Gimbal lock detected. Setting third angle to zero "
                          "since it is not possible to uniquely determine "
                          "all angles.")

    # Reverse role of extrinsic and intrinsic rotations, but let third angle be
    # zero for gimbal locked cases
    if extrinsic:
        angles = angles[:, ::-1]
    return angles

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline void _compose_quat_single( # calculate p * q into r
    const double[:] p, const double[:] q, double[:] r
):
    cdef double[:] cross = _cross3(p[:3], q[:3])

    r[0] = p[3]*q[0] + q[3]*p[0] + cross[0]
    r[1] = p[3]*q[1] + q[3]*p[1] + cross[1]
    r[2] = p[3]*q[2] + q[3]*p[2] + cross[2]
    r[3] = p[3]*q[3] - p[0]*q[0] - p[1]*q[1] - p[2]*q[2]

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline double[:, :] _compose_quat(
    const double[:, :] p, const double[:, :] q
):
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
):
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
):
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

cdef class Rotation(object):
    """Rotation in 3 dimensions.

    This class provides an interface to initialize from and represent rotations
    with:

    - Quaternions
    - Rotation Matrices
    - Rotation Vectors
    - Modified Rodrigues Parameters
    - Euler Angles

    The following operations on rotations are supported:

    - Application on vectors
    - Rotation Composition
    - Rotation Inversion
    - Rotation Indexing

    Indexing within a rotation is supported since multiple rotation transforms
    can be stored within a single `Rotation` instance.

    To create `Rotation` objects use ``from_...`` methods (see examples below).
    ``Rotation(...)`` is not supposed to be instantiated directly.

    Attributes
    ----------
    single

    Methods
    -------
    __len__
    from_quat
    from_matrix
    from_rotvec
    from_mrp
    from_euler
    as_quat
    as_matrix
    as_rotvec
    as_mrp
    as_euler
    apply
    __mul__
    inv
    magnitude
    mean
    reduce
    create_group
    __getitem__
    identity
    random
    align_vectors

    See Also
    --------
    Slerp

    Notes
    -----
    .. versionadded: 1.2.0

    Examples
    --------
    >>> from scipy.spatial.transform import Rotation as R

    A `Rotation` instance can be initialized in any of the above formats and
    converted to any of the others. The underlying object is independent of the
    representation used for initialization.

    Consider a counter-clockwise rotation of 90 degrees about the z-axis. This
    corresponds to the following quaternion (in scalar-last format):

    >>> r = R.from_quat([0, 0, np.sin(np.pi/4), np.cos(np.pi/4)])

    The rotation can be expressed in any of the other formats:

    >>> r.as_matrix()
    array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
    [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
    [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
    >>> r.as_rotvec()
    array([0.        , 0.        , 1.57079633])
    >>> r.as_euler('zyx', degrees=True)
    array([90.,  0.,  0.])

    The same rotation can be initialized using a rotation matrix:

    >>> r = R.from_matrix([[0, -1, 0],
    ...                    [1, 0, 0],
    ...                    [0, 0, 1]])

    Representation in other formats:

    >>> r.as_quat()
    array([0.        , 0.        , 0.70710678, 0.70710678])
    >>> r.as_rotvec()
    array([0.        , 0.        , 1.57079633])
    >>> r.as_euler('zyx', degrees=True)
    array([90.,  0.,  0.])

    The rotation vector corresponding to this rotation is given by:

    >>> r = R.from_rotvec(np.pi/2 * np.array([0, 0, 1]))

    Representation in other formats:

    >>> r.as_quat()
    array([0.        , 0.        , 0.70710678, 0.70710678])
    >>> r.as_matrix()
    array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
           [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
           [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
    >>> r.as_euler('zyx', degrees=True)
    array([90.,  0.,  0.])

    The ``from_euler`` method is quite flexible in the range of input formats
    it supports. Here we initialize a single rotation about a single axis:

    >>> r = R.from_euler('z', 90, degrees=True)

    Again, the object is representation independent and can be converted to any
    other format:

    >>> r.as_quat()
    array([0.        , 0.        , 0.70710678, 0.70710678])
    >>> r.as_matrix()
    array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
           [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
           [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
    >>> r.as_rotvec()
    array([0.        , 0.        , 1.57079633])

    It is also possible to initialize multiple rotations in a single instance
    using any of the `from_...` functions. Here we initialize a stack of 3
    rotations using the ``from_euler`` method:

    >>> r = R.from_euler('zyx', [
    ... [90, 0, 0],
    ... [0, 45, 0],
    ... [45, 60, 30]], degrees=True)

    The other representations also now return a stack of 3 rotations. For
    example:

    >>> r.as_quat()
    array([[0.        , 0.        , 0.70710678, 0.70710678],
           [0.        , 0.38268343, 0.        , 0.92387953],
           [0.39190384, 0.36042341, 0.43967974, 0.72331741]])

    Applying the above rotations onto a vector:

    >>> v = [1, 2, 3]
    >>> r.apply(v)
    array([[-2.        ,  1.        ,  3.        ],
           [ 2.82842712,  2.        ,  1.41421356],
           [ 2.24452282,  0.78093109,  2.89002836]])

    A `Rotation` instance can be indexed and sliced as if it were a single
    1D array or list:

    >>> r.as_quat()
    array([[0.        , 0.        , 0.70710678, 0.70710678],
           [0.        , 0.38268343, 0.        , 0.92387953],
           [0.39190384, 0.36042341, 0.43967974, 0.72331741]])
    >>> p = r[0]
    >>> p.as_matrix()
    array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
           [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
           [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
    >>> q = r[1:3]
    >>> q.as_quat()
    array([[0.        , 0.38268343, 0.        , 0.92387953],
           [0.39190384, 0.36042341, 0.43967974, 0.72331741]])

    In fact it can be converted to numpy.array:

    >>> r_array = np.asarray(r)
    >>> r_array.shape
    (3,)
    >>> r_array[0].as_matrix()
    array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
           [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
           [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])

    Multiple rotations can be composed using the ``*`` operator:

    >>> r1 = R.from_euler('z', 90, degrees=True)
    >>> r2 = R.from_rotvec([np.pi/4, 0, 0])
    >>> v = [1, 2, 3]
    >>> r2.apply(r1.apply(v))
    array([-2.        , -1.41421356,  2.82842712])
    >>> r3 = r2 * r1 # Note the order
    >>> r3.apply(v)
    array([-2.        , -1.41421356,  2.82842712])

    Finally, it is also possible to invert rotations:

    >>> r1 = R.from_euler('z', [90, 45], degrees=True)
    >>> r2 = r1.inv()
    >>> r2.as_euler('zyx', degrees=True)
    array([[-90.,   0.,   0.],
           [-45.,   0.,   0.]])

    These examples serve as an overview into the `Rotation` class and highlight
    major functionalities. For more thorough examples of the range of input and
    output formats supported, consult the individual method's examples.

    """
    cdef double[:, :] _quat
    cdef bint _single

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __init__(self, quat, normalize=True, copy=True):
        self._single = False
        quat = np.asarray(quat, dtype=float)

        if quat.ndim not in [1, 2] or quat.shape[len(quat.shape) - 1] != 4:
            raise ValueError("Expected `quat` to have shape (4,) or (N x 4), "
                             "got {}.".format(quat.shape))

        # If a single quaternion is given, convert it to a 2D 1 x 4 matrix but
        # set self._single to True so that we can return appropriate objects
        # in the `to_...` methods
        if quat.shape == (4,):
            quat = quat[None, :]
            self._single = True

        cdef Py_ssize_t num_rotations = quat.shape[0]
        if normalize:
            self._quat = quat.copy()
            for ind in range(num_rotations):
                if isnan(_normalize4(self._quat[ind, :])):
                    raise ValueError("Found zero norm quaternions in `quat`.")
        else:
            self._quat = quat.copy() if copy else quat

    @property
    def single(self):
        """Whether this instance represents a single rotation."""
        return self._single

    def __len__(self):
        """Number of rotations contained in this object.

        Multiple rotations can be stored in a single instance.

        Returns
        -------
        length : int
            Number of rotations stored in object.

        Raises
        ------
        TypeError if the instance was created as a single rotation.
        """
        if self._single:
            raise TypeError("Single rotation has no len().")

        return self._quat.shape[0]

    @classmethod
    def from_quat(cls, quat):
        """Initialize from quaternions.

        3D rotations can be represented using unit-norm quaternions [1]_.

        Parameters
        ----------
        quat : array_like, shape (N, 4) or (4,)
            Each row is a (possibly non-unit norm) quaternion in scalar-last
            (x, y, z, w) format. Each quaternion will be normalized to unit
            norm.

        Returns
        -------
        rotation : `Rotation` instance
            Object containing the rotations represented by input quaternions.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Initialize a single rotation:

        >>> r = R.from_quat([1, 0, 0, 0])
        >>> r.as_quat()
        array([1., 0., 0., 0.])
        >>> r.as_quat().shape
        (4,)

        Initialize multiple rotations in a single object:

        >>> r = R.from_quat([
        ... [1, 0, 0, 0],
        ... [0, 0, 0, 1]
        ... ])
        >>> r.as_quat()
        array([[1., 0., 0., 0.],
               [0., 0., 0., 1.]])
        >>> r.as_quat().shape
        (2, 4)

        It is also possible to have a stack of a single rotation:

        >>> r = R.from_quat([[0, 0, 0, 1]])
        >>> r.as_quat()
        array([[0., 0., 0., 1.]])
        >>> r.as_quat().shape
        (1, 4)

        Quaternions are normalized before initialization.

        >>> r = R.from_quat([0, 0, 1, 1])
        >>> r.as_quat()
        array([0.        , 0.        , 0.70710678, 0.70710678])
        """
        return cls(quat, normalize=True)

    @classmethod
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def from_matrix(cls, matrix):
        """Initialize from rotation matrix.

        Rotations in 3 dimensions can be represented with 3 x 3 proper
        orthogonal matrices [1]_. If the input is not proper orthogonal,
        an approximation is created using the method described in [2]_.

        Parameters
        ----------
        matrix : array_like, shape (N, 3, 3) or (3, 3)
            A single matrix or a stack of matrices, where ``matrix[i]`` is
            the i-th matrix.

        Returns
        -------
        rotation : `Rotation` instance
            Object containing the rotations represented by the rotation
            matrices.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
        .. [2] F. Landis Markley, "Unit Quaternion from Rotation Matrix",
               Journal of guidance, control, and dynamics vol. 31.2, pp.
               440-442, 2008.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Initialize a single rotation:

        >>> r = R.from_matrix([
        ... [0, -1, 0],
        ... [1, 0, 0],
        ... [0, 0, 1]])
        >>> r.as_matrix().shape
        (3, 3)

        Initialize multiple rotations in a single object:

        >>> r = R.from_matrix([
        ... [
        ...     [0, -1, 0],
        ...     [1, 0, 0],
        ...     [0, 0, 1],
        ... ],
        ... [
        ...     [1, 0, 0],
        ...     [0, 0, -1],
        ...     [0, 1, 0],
        ... ]])
        >>> r.as_matrix().shape
        (2, 3, 3)

        If input matrices are not special orthogonal (orthogonal with
        determinant equal to +1), then a special orthogonal estimate is stored:

        >>> a = np.array([
        ... [0, -0.5, 0],
        ... [0.5, 0, 0],
        ... [0, 0, 0.5]])
        >>> np.linalg.det(a)
        0.12500000000000003
        >>> r = R.from_matrix(a)
        >>> matrix = r.as_matrix()
        >>> matrix
        array([[-0.38461538, -0.92307692,  0.        ],
               [ 0.92307692, -0.38461538,  0.        ],
               [ 0.        ,  0.        ,  1.        ]])
        >>> np.linalg.det(matrix)
        1.0000000000000002

        It is also possible to have a stack containing a single rotation:

        >>> r = R.from_matrix([[
        ... [0, -1, 0],
        ... [1, 0, 0],
        ... [0, 0, 1]]])
        >>> r.as_matrix()
        array([[[ 0., -1.,  0.],
                [ 1.,  0.,  0.],
                [ 0.,  0.,  1.]]])
        >>> r.as_matrix().shape
        (1, 3, 3)

        Notes
        -----
        This function was called from_dcm before.

        .. versionadded:: 1.4.0
        """
        is_single = False
        matrix = np.asarray(matrix, dtype=float)

        if (matrix.ndim not in [2, 3] or
            matrix.shape[len(matrix.shape)-2:] != (3, 3)):
            raise ValueError("Expected `matrix` to have shape (3, 3) or "
                             "(N, 3, 3), got {}".format(matrix.shape))

        # If a single matrix is given, convert it to 3D 1 x 3 x 3 matrix but
        # set self._single to True so that we can return appropriate objects in
        # the `to_...` methods
        cdef double[:, :, :] cmatrix
        if matrix.shape == (3, 3):
            cmatrix = matrix[None, :, :]
            is_single = True
        else:
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
            return cls(quat[0], normalize=False, copy=False)
        else:
            return cls(quat, normalize=False, copy=False)

    @classmethod
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def from_rotvec(cls, rotvec):
        """Initialize from rotation vectors.

        A rotation vector is a 3 dimensional vector which is co-directional to
        the axis of rotation and whose norm gives the angle of rotation (in
        radians) [1]_.

        Parameters
        ----------
        rotvec : array_like, shape (N, 3) or (3,)
            A single vector or a stack of vectors, where `rot_vec[i]` gives
            the ith rotation vector.

        Returns
        -------
        rotation : `Rotation` instance
            Object containing the rotations represented by input rotation
            vectors.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation#Rotation_vector

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Initialize a single rotation:

        >>> r = R.from_rotvec(np.pi/2 * np.array([0, 0, 1]))
        >>> r.as_rotvec()
        array([0.        , 0.        , 1.57079633])
        >>> r.as_rotvec().shape
        (3,)

        Initialize multiple rotations in one object:

        >>> r = R.from_rotvec([
        ... [0, 0, np.pi/2],
        ... [np.pi/2, 0, 0]])
        >>> r.as_rotvec()
        array([[0.        , 0.        , 1.57079633],
               [1.57079633, 0.        , 0.        ]])
        >>> r.as_rotvec().shape
        (2, 3)

        It is also possible to have a stack of a single rotaton:

        >>> r = R.from_rotvec([[0, 0, np.pi/2]])
        >>> r.as_rotvec().shape
        (1, 3)

        """
        is_single = False
        rotvec = np.asarray(rotvec, dtype=float)

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

            if angle <= 1e-3:  # small angle
                angle2 = angle * angle
                scale = 0.5 - angle2 / 48 + angle2 * angle2 / 3840
            else:  # large angle
                scale = sin(angle / 2) / angle

            quat[ind, 0] = scale * crotvec[ind, 0]
            quat[ind, 1] = scale * crotvec[ind, 1]
            quat[ind, 2] = scale * crotvec[ind, 2]
            quat[ind, 3] = cos(angle / 2)

        if is_single:
            return cls(quat[0], normalize=False, copy=False)
        else:
            return cls(quat, normalize=False, copy=False)

    @classmethod
    def from_euler(cls, seq, angles, degrees=False):
        """Initialize from Euler angles.

        Rotations in 3-D can be represented by a sequence of 3
        rotations around a sequence of axes. In theory, any three axes spanning
        the 3-D Euclidean space are enough. In practice, the axes of rotation are
        chosen to be the basis vectors.

        The three rotations can either be in a global frame of reference
        (extrinsic) or in a body centred frame of reference (intrinsic), which
        is attached to, and moves with, the object under rotation [1]_.

        Parameters
        ----------
        seq : string
            Specifies sequence of axes for rotations. Up to 3 characters
            belonging to the set {'X', 'Y', 'Z'} for intrinsic rotations, or
            {'x', 'y', 'z'} for extrinsic rotations. Extrinsic and intrinsic
            rotations cannot be mixed in one function call.
        angles : float or array_like, shape (N,) or (N, [1 or 2 or 3])
            Euler angles specified in radians (`degrees` is False) or degrees
            (`degrees` is True).
            For a single character `seq`, `angles` can be:

            - a single value
            - array_like with shape (N,), where each `angle[i]`
              corresponds to a single rotation
            - array_like with shape (N, 1), where each `angle[i, 0]`
              corresponds to a single rotation

            For 2- and 3-character wide `seq`, `angles` can be:

            - array_like with shape (W,) where `W` is the width of
              `seq`, which corresponds to a single rotation with `W` axes
            - array_like with shape (N, W) where each `angle[i]`
              corresponds to a sequence of Euler angles describing a single
              rotation

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

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Initialize a single rotation along a single axis:

        >>> r = R.from_euler('x', 90, degrees=True)
        >>> r.as_quat().shape
        (4,)

        Initialize a single rotation with a given axis sequence:

        >>> r = R.from_euler('zyx', [90, 45, 30], degrees=True)
        >>> r.as_quat().shape
        (4,)

        Initialize a stack with a single rotation around a single axis:

        >>> r = R.from_euler('x', [90], degrees=True)
        >>> r.as_quat().shape
        (1, 4)

        Initialize a stack with a single rotation with an axis sequence:

        >>> r = R.from_euler('zyx', [[90, 45, 30]], degrees=True)
        >>> r.as_quat().shape
        (1, 4)

        Initialize multiple elementary rotations in one object:

        >>> r = R.from_euler('x', [90, 45, 30], degrees=True)
        >>> r.as_quat().shape
        (3, 4)

        Initialize multiple rotations in one object:

        >>> r = R.from_euler('zyx', [[90, 45, 30], [35, 45, 90]], degrees=True)
        >>> r.as_quat().shape
        (2, 4)

        """
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
                                 "of axes specified, got {} for shape".format(
                                 angles.shape))

            if angles.ndim == 1:
                # (1, num_axes)
                angles = angles[None, :]
                is_single = True

        # By now angles should have shape (num_rot, num_axes)
        # sanity check
        if angles.ndim != 2 or angles.shape[-1] != num_axes:
            raise ValueError("Expected angles to have shape (num_rotations, "
                             "num_axes), got {}.".format(angles.shape))

        quat = _elementary_quat_compose(seq.encode(), angles, intrinsic)

        if is_single:
            return cls(quat[0], normalize=False, copy=False)
        else:
            return cls(quat, normalize=False, copy=False)

    @classmethod
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def from_mrp(cls, mrp):
        """Initialize from Modified Rodrigues Parameters (MRPs).

        MRPs are a 3 dimensional vector co-directional to the axis of rotation and whose
        magnitude is equal to ``tan(theta / 4)``, where ``theta`` is the angle of rotation
        (in radians) [1]_.

        MRPs have a singuarity at 360 degrees which can be avoided by ensuring the angle of
        rotation does not exceed 180 degrees, i.e. switching the direction of the rotation when
        it is past 180 degrees.

        Parameters
        ----------
        mrp : array_like, shape (N, 3) or (3,)
            A single vector or a stack of vectors, where `mrp[i]` gives
            the ith set of MRPs.

        Returns
        -------
        rotation : `Rotation` instance
            Object containing the rotations represented by input MRPs.

        References
        ----------
        .. [1] Shuster, M. D. "A Survery of Attitude Representations",
               The Journal of Astronautical Sciences, Vol. 41, No.4, 1993,
               pp. 475-476

        Notes
        -----

        .. versionadded:: 1.6.0

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Initialize a single rotation:

        >>> r = R.from_mrp([0, 0, 1])
        >>> r.as_euler('xyz', degrees=True)
        array([0.        , 0.        , 180.      ])
        >>> r.as_euler('xyz').shape
        (3,)

        Initialize multiple rotations in one object:

        >>> r = R.from_mrp([
        ... [0, 0, 1],
        ... [1, 0, 0]])
        >>> r.as_euler('xyz', degrees=True)
        array([[0.        , 0.        , 180.      ],
               [180.0     , 0.        , 0.        ]])
        >>> r.as_euler('xyz').shape
        (2, 3)

        It is also possible to have a stack of a single rotation:

        >>> r = R.from_mrp([[0, 0, np.pi/2]])
        >>> r.as_euler('xyz').shape
        (1, 3)

        """
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
            return cls(quat[0], normalize=False, copy=False)
        else:
            return cls(quat, normalize=False, copy=False)

    def as_quat(self):
        """Represent as quaternions.

        Rotations in 3 dimensions can be represented using unit norm
        quaternions [1]_. The mapping from quaternions to rotations is
        two-to-one, i.e. quaternions ``q`` and ``-q``, where ``-q`` simply
        reverses the sign of each component, represent the same spatial
        rotation. The returned value is in scalar-last (x, y, z, w) format.

        Returns
        -------
        quat : `numpy.ndarray`, shape (4,) or (N, 4)
            Shape depends on shape of inputs used for initialization.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Represent a single rotation:

        >>> r = R.from_matrix([[0, -1, 0],
        ...                    [1, 0, 0],
        ...                    [0, 0, 1]])
        >>> r.as_quat()
        array([0.        , 0.        , 0.70710678, 0.70710678])
        >>> r.as_quat().shape
        (4,)

        Represent a stack with a single rotation:

        >>> r = R.from_quat([[0, 0, 0, 1]])
        >>> r.as_quat().shape
        (1, 4)

        Represent multiple rotations in a single object:

        >>> r = R.from_rotvec([[np.pi, 0, 0], [0, 0, np.pi/2]])
        >>> r.as_quat().shape
        (2, 4)

        """
        if self._single:
            return np.array(self._quat[0], copy=True)
        else:
            return np.array(self._quat, copy=True)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def as_matrix(self):
        """Represent as rotation matrix.

        3D rotations can be represented using rotation matrices, which
        are 3 x 3 real orthogonal matrices with determinant equal to +1 [1]_.

        Returns
        -------
        matrix : ndarray, shape (3, 3) or (N, 3, 3)
            Shape depends on shape of inputs used for initialization.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Represent a single rotation:

        >>> r = R.from_rotvec([0, 0, np.pi/2])
        >>> r.as_matrix()
        array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
               [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
               [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
        >>> r.as_matrix().shape
        (3, 3)

        Represent a stack with a single rotation:

        >>> r = R.from_quat([[1, 1, 0, 0]])
        >>> r.as_matrix()
        array([[[ 0.,  1.,  0.],
                [ 1.,  0.,  0.],
                [ 0.,  0., -1.]]])
        >>> r.as_matrix().shape
        (1, 3, 3)

        Represent multiple rotations:

        >>> r = R.from_rotvec([[np.pi/2, 0, 0], [0, 0, np.pi/2]])
        >>> r.as_matrix()
        array([[[ 1.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                [ 0.00000000e+00,  2.22044605e-16, -1.00000000e+00],
                [ 0.00000000e+00,  1.00000000e+00,  2.22044605e-16]],
               [[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
                [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
                [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]]])
        >>> r.as_matrix().shape
        (2, 3, 3)

        Notes
        -----
        This function was called as_dcm before.

        .. versionadded:: 1.4.0
        """
        cdef double[:, :] quat = self._quat
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

        ret = np.asarray(matrix)
        if self._single:
            return ret[0]
        else:
            return ret

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def as_rotvec(self):
        """Represent as rotation vectors.

        A rotation vector is a 3 dimensional vector which is co-directional to
        the axis of rotation and whose norm gives the angle of rotation (in
        radians) [1]_.

        Returns
        -------
        rotvec : ndarray, shape (3,) or (N, 3)
            Shape depends on shape of inputs used for initialization.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation#Rotation_vector

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Represent a single rotation:

        >>> r = R.from_euler('z', 90, degrees=True)
        >>> r.as_rotvec()
        array([0.        , 0.        , 1.57079633])
        >>> r.as_rotvec().shape
        (3,)

        Represent a stack with a single rotation:

        >>> r = R.from_quat([[0, 0, 1, 1]])
        >>> r.as_rotvec()
        array([[0.        , 0.        , 1.57079633]])
        >>> r.as_rotvec().shape
        (1, 3)

        Represent multiple rotations in a single object:

        >>> r = R.from_quat([[0, 0, 1, 1], [1, 1, 0, 1]])
        >>> r.as_rotvec()
        array([[0.        , 0.        , 1.57079633],
               [1.35102172, 1.35102172, 0.        ]])
        >>> r.as_rotvec().shape
        (2, 3)

        """
        cdef Py_ssize_t num_rotations = len(self._quat)
        cdef double angle, scale, angle2
        cdef double[:, :] rotvec = _empty2(num_rotations, 3)
        cdef double[:] quat

        for ind in range(num_rotations):
            if self._quat[ind, 3] < 0:  # w > 0 to ensure 0 <= angle <= pi
                quat = self._quat[ind, :].copy()
                for i in range(4):
                    quat[i] *= -1
            else:
                quat = self._quat[ind, :]

            angle = 2 * atan2(_norm3(quat), quat[3])

            if angle <= 1e-3:  # small angle
                angle2 = angle * angle
                scale = 2 + angle2 / 12 + 7 * angle2 * angle2 / 2880
            else:  # large angle
                scale = angle / sin(angle / 2)

            rotvec[ind, 0] = scale * quat[0]
            rotvec[ind, 1] = scale * quat[1]
            rotvec[ind, 2] = scale * quat[2]

        if self._single:
            return np.asarray(rotvec[0])
        else:
            return np.asarray(rotvec)

    def as_euler(self, seq, degrees=False):
        """Represent as Euler angles.

        Any orientation can be expressed as a composition of 3 elementary
        rotations. Once the axis sequence has been chosen, Euler angles define
        the angle of rotation around each respective axis [1]_.

        The algorithm from [2]_ has been used to calculate Euler angles for the
        rotation about a given sequence of axes.

        Euler angles suffer from the problem of gimbal lock [3]_, where the
        representation loses a degree of freedom and it is not possible to
        determine the first and third angles uniquely. In this case,
        a warning is raised, and the third angle is set to zero. Note however
        that the returned angles still represent the correct rotation.

        Parameters
        ----------
        seq : string, length 3
            3 characters belonging to the set {'X', 'Y', 'Z'} for intrinsic
            rotations, or {'x', 'y', 'z'} for extrinsic rotations [1]_.
            Adjacent axes cannot be the same.
            Extrinsic and intrinsic rotations cannot be mixed in one function
            call.
        degrees : boolean, optional
            Returned angles are in degrees if this flag is True, else they are
            in radians. Default is False.

        Returns
        -------
        angles : ndarray, shape (3,) or (N, 3)
            Shape depends on shape of inputs used to initialize object.
            The returned angles are in the range:

            - First angle belongs to [-180, 180] degrees (both inclusive)
            - Third angle belongs to [-180, 180] degrees (both inclusive)
            - Second angle belongs to:

                - [-90, 90] degrees if all axes are different (like xyz)
                - [0, 180] degrees if first and third axes are the same
                  (like zxz)

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Euler_angles#Definition_by_intrinsic_rotations
        .. [2] Malcolm D. Shuster, F. Landis Markley, "General formula for
               extraction the Euler angles", Journal of guidance, control, and
               dynamics, vol. 29.1, pp. 215-221. 2006
        .. [3] https://en.wikipedia.org/wiki/Gimbal_lock#In_applied_mathematics

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Represent a single rotation:

        >>> r = R.from_rotvec([0, 0, np.pi/2])
        >>> r.as_euler('zxy', degrees=True)
        array([90.,  0.,  0.])
        >>> r.as_euler('zxy', degrees=True).shape
        (3,)

        Represent a stack of single rotation:

        >>> r = R.from_rotvec([[0, 0, np.pi/2]])
        >>> r.as_euler('zxy', degrees=True)
        array([[90.,  0.,  0.]])
        >>> r.as_euler('zxy', degrees=True).shape
        (1, 3)

        Represent multiple rotations in a single object:

        >>> r = R.from_rotvec([
        ... [0, 0, np.pi/2],
        ... [0, -np.pi/3, 0],
        ... [np.pi/4, 0, 0]])
        >>> r.as_euler('zxy', degrees=True)
        array([[ 90.,   0.,   0.],
               [  0.,   0., -60.],
               [  0.,  45.,   0.]])
        >>> r.as_euler('zxy', degrees=True).shape
        (3, 3)

        """
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

        matrix = self.as_matrix()
        if matrix.ndim == 2:
            matrix = matrix[None, :, :]
        angles = np.asarray(_compute_euler_from_matrix(
            matrix, seq.encode(), extrinsic))
        if degrees:
            angles = np.rad2deg(angles)

        return angles[0] if self._single else angles

    def as_mrp(self):
        """Represent as Modified Rodrigues Parameters (MRPs).

        MRPs are a 3 dimensional vector co-directional to the axis of rotation and whose
        magnitude is equal to ``tan(theta / 4)``, where ``theta`` is the angle of rotation
        (in radians) [1]_.

        MRPs have a singuarity at 360 degrees which can be avoided by ensuring the angle of
        rotation does not exceed 180 degrees, i.e. switching the direction of the rotation when
        it is past 180 degrees. This function will always return MRPs corresponding to a rotation
        of less than or equal to 180 degrees.

        Returns
        -------
        mrps : ndarray, shape (3,) or (N, 3)
            Shape depends on shape of inputs used for initialization.

        References
        ----------
        .. [1] Shuster, M. D. "A Survery of Attitude Representations",
               The Journal of Astronautical Sciences, Vol. 41, No.4, 1993,
               pp. 475-476

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Represent a single rotation:

        >>> r = R.from_rotvec([0, 0, np.pi])
        >>> r.as_mrp()
        array([0.        , 0.        , 1.         ])
        >>> r.as_mrp().shape
        (3,)

        Represent a stack with a single rotation:

        >>> r = R.from_euler('xyz', [[180, 0, 0]], degrees=True)
        >>> r.as_mrp()
        array([[1.       , 0.        , 0.         ]])
        >>> r.as_mrp().shape
        (1, 3)

        Represent multiple rotations:

        >>> r = R.from_rotvec([[np.pi/2, 0, 0], [0, 0, np.pi/2]])
        >>> r.as_mrp()
        array([[0.41421356, 0.        , 0.        ],
               [0.        , 0.        , 0.41421356]])
        >>> r.as_mrp().shape
        (2, 3)

        Notes
        -----

        .. versionadded:: 1.6.0
        """
        cdef Py_ssize_t num_rotations = len(self._quat)
        cdef double[:, :] mrps = _empty2(num_rotations, 3)
        cdef int sign
        cdef double denominator

        for ind in range(num_rotations):

            # Ensure we are calculating the set of MRPs that correspond
            # to a rotation of <= 180
            sign = -1 if self._quat[ind, 3] < 0 else 1

            denominator = 1 + sign * self._quat[ind, 3]
            for i in range(3):
                mrps[ind, i] = sign * self._quat[ind, i] / denominator

        if self._single:
            return np.asarray(mrps[0])
        else:
            return np.asarray(mrps)

    def apply(self, vectors, inverse=False):
        """Apply this rotation to a set of vectors.

        If the original frame rotates to the final frame by this rotation, then
        its application to a vector can be seen in two ways:

            - As a projection of vector components expressed in the final frame
              to the original frame.
            - As the physical rotation of a vector being glued to the original
              frame as it rotates. In this case the vector components are
              expressed in the original frame before and after the rotation.

        In terms of rotation matricies, this application is the same as
        ``self.as_matrix().dot(vectors)``.

        Parameters
        ----------
        vectors : array_like, shape (3,) or (N, 3)
            Each `vectors[i]` represents a vector in 3D space. A single vector
            can either be specified with shape `(3, )` or `(1, 3)`. The number
            of rotations and number of vectors given must follow standard numpy
            broadcasting rules: either one of them equals unity or they both
            equal each other.
        inverse : boolean, optional
            If True then the inverse of the rotation(s) is applied to the input
            vectors. Default is False.

        Returns
        -------
        rotated_vectors : ndarray, shape (3,) or (N, 3)
            Result of applying rotation on input vectors.
            Shape depends on the following cases:

                - If object contains a single rotation (as opposed to a stack
                  with a single rotation) and a single vector is specified with
                  shape ``(3,)``, then `rotated_vectors` has shape ``(3,)``.
                - In all other cases, `rotated_vectors` has shape ``(N, 3)``,
                  where ``N`` is either the number of rotations or vectors.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Single rotation applied on a single vector:

        >>> vector = np.array([1, 0, 0])
        >>> r = R.from_rotvec([0, 0, np.pi/2])
        >>> r.as_matrix()
        array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
               [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
               [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
        >>> r.apply(vector)
        array([2.22044605e-16, 1.00000000e+00, 0.00000000e+00])
        >>> r.apply(vector).shape
        (3,)

        Single rotation applied on multiple vectors:

        >>> vectors = np.array([
        ... [1, 0, 0],
        ... [1, 2, 3]])
        >>> r = R.from_rotvec([0, 0, np.pi/4])
        >>> r.as_matrix()
        array([[ 0.70710678, -0.70710678,  0.        ],
               [ 0.70710678,  0.70710678,  0.        ],
               [ 0.        ,  0.        ,  1.        ]])
        >>> r.apply(vectors)
        array([[ 0.70710678,  0.70710678,  0.        ],
               [-0.70710678,  2.12132034,  3.        ]])
        >>> r.apply(vectors).shape
        (2, 3)

        Multiple rotations on a single vector:

        >>> r = R.from_rotvec([[0, 0, np.pi/4], [np.pi/2, 0, 0]])
        >>> vector = np.array([1,2,3])
        >>> r.as_matrix()
        array([[[ 7.07106781e-01, -7.07106781e-01,  0.00000000e+00],
                [ 7.07106781e-01,  7.07106781e-01,  0.00000000e+00],
                [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]],
               [[ 1.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                [ 0.00000000e+00,  2.22044605e-16, -1.00000000e+00],
                [ 0.00000000e+00,  1.00000000e+00,  2.22044605e-16]]])
        >>> r.apply(vector)
        array([[-0.70710678,  2.12132034,  3.        ],
               [ 1.        , -3.        ,  2.        ]])
        >>> r.apply(vector).shape
        (2, 3)

        Multiple rotations on multiple vectors. Each rotation is applied on the
        corresponding vector:

        >>> r = R.from_euler('zxy', [
        ... [0, 0, 90],
        ... [45, 30, 60]], degrees=True)
        >>> vectors = [
        ... [1, 2, 3],
        ... [1, 0, -1]]
        >>> r.apply(vectors)
        array([[ 3.        ,  2.        , -1.        ],
               [-0.09026039,  1.11237244, -0.86860844]])
        >>> r.apply(vectors).shape
        (2, 3)

        It is also possible to apply the inverse rotation:

        >>> r = R.from_euler('zxy', [
        ... [0, 0, 90],
        ... [45, 30, 60]], degrees=True)
        >>> vectors = [
        ... [1, 2, 3],
        ... [1, 0, -1]]
        >>> r.apply(vectors, inverse=True)
        array([[-3.        ,  2.        ,  1.        ],
               [ 1.09533535, -0.8365163 ,  0.3169873 ]])

        """
        vectors = np.asarray(vectors)
        if vectors.ndim > 2 or vectors.shape[-1] != 3:
            raise ValueError("Expected input of shape (3,) or (P, 3), "
                             "got {}.".format(vectors.shape))

        single_vector = False
        if vectors.shape == (3,):
            single_vector = True
            vectors = vectors[None, :]

        matrix = self.as_matrix()
        if self._single:
            matrix = matrix[None, :, :]

        n_vectors = vectors.shape[0]
        n_rotations = len(self._quat)

        if n_vectors != 1 and n_rotations != 1 and n_vectors != n_rotations:
            raise ValueError("Expected equal numbers of rotations and vectors "
                             ", or a single rotation, or a single vector, got "
                             "{} rotations and {} vectors.".format(
                                n_rotations, n_vectors))

        if inverse:
            result = np.einsum('ikj,ik->ij', matrix, vectors)
        else:
            result = np.einsum('ijk,ik->ij', matrix, vectors)

        if self._single and single_vector:
            return result[0]
        else:
            return result

    def __mul__(Rotation self, Rotation other):
        """Compose this rotation with the other.

        If `p` and `q` are two rotations, then the composition of 'q followed
        by p' is equivalent to `p * q`. In terms of rotation matrices,
        the composition can be expressed as
        ``p.as_matrix().dot(q.as_matrix())``.

        Parameters
        ----------
        other : `Rotation` instance
            Object containing the rotations to be composed with this one. Note
            that rotation compositions are not commutative, so ``p * q`` is
            different from ``q * p``.

        Returns
        -------
        composition : `Rotation` instance
            This function supports composition of multiple rotations at a time.
            The following cases are possible:

            - Either ``p`` or ``q`` contains a single rotation. In this case
              `composition` contains the result of composing each rotation in
              the other object with the single rotation.
            - Both ``p`` and ``q`` contain ``N`` rotations. In this case each
              rotation ``p[i]`` is composed with the corresponding rotation
              ``q[i]`` and `output` contains ``N`` rotations.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Composition of two single rotations:

        >>> p = R.from_quat([0, 0, 1, 1])
        >>> q = R.from_quat([1, 0, 0, 1])
        >>> p.as_matrix()
        array([[ 0., -1.,  0.],
               [ 1.,  0.,  0.],
               [ 0.,  0.,  1.]])
        >>> q.as_matrix()
        array([[ 1.,  0.,  0.],
               [ 0.,  0., -1.],
               [ 0.,  1.,  0.]])
        >>> r = p * q
        >>> r.as_matrix()
        array([[0., 0., 1.],
               [1., 0., 0.],
               [0., 1., 0.]])

        Composition of two objects containing equal number of rotations:

        >>> p = R.from_quat([[0, 0, 1, 1], [1, 0, 0, 1]])
        >>> q = R.from_rotvec([[np.pi/4, 0, 0], [-np.pi/4, 0, np.pi/4]])
        >>> p.as_quat()
        array([[0.        , 0.        , 0.70710678, 0.70710678],
               [0.70710678, 0.        , 0.        , 0.70710678]])
        >>> q.as_quat()
        array([[ 0.38268343,  0.        ,  0.        ,  0.92387953],
               [-0.37282173,  0.        ,  0.37282173,  0.84971049]])
        >>> r = p * q
        >>> r.as_quat()
        array([[ 0.27059805,  0.27059805,  0.65328148,  0.65328148],
               [ 0.33721128, -0.26362477,  0.26362477,  0.86446082]])

        """
        len_self = len(self._quat)
        len_other = len(other._quat)
        if not(len_self == 1 or len_other == 1 or len_self == len_other):
            raise ValueError("Expected equal number of rotations in both "
                             "or a single rotation in either object, "
                             "got {} rotations in first and {} rotations in "
                             "second object.".format(
                                len(self), len(other)))
        result = _compose_quat(self._quat, other._quat)
        if self._single and other._single:
            result = result[0]
        return self.__class__(result, normalize=True, copy=False)

    def inv(self):
        """Invert this rotation.

        Composition of a rotation with its inverse results in an identity
        transformation.

        Returns
        -------
        inverse : `Rotation` instance
            Object containing inverse of the rotations in the current instance.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Inverting a single rotation:

        >>> p = R.from_euler('z', 45, degrees=True)
        >>> q = p.inv()
        >>> q.as_euler('zyx', degrees=True)
        array([-45.,   0.,   0.])

        Inverting multiple rotations:

        >>> p = R.from_rotvec([[0, 0, np.pi/3], [-np.pi/4, 0, 0]])
        >>> q = p.inv()
        >>> q.as_rotvec()
        array([[-0.        , -0.        , -1.04719755],
               [ 0.78539816, -0.        , -0.        ]])

        """
        cdef np.ndarray quat = np.array(self._quat, copy=True)
        quat[:, -1] *= -1
        if self._single:
            quat = quat[0]
        return self.__class__(quat, copy=False)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def magnitude(self):
        """Get the magnitude(s) of the rotation(s).

        Returns
        -------
        magnitude : ndarray or float
            Angle(s) in radians, float if object contains a single rotation
            and ndarray if object contains multiple rotations.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> r = R.from_quat(np.eye(4))
        >>> r.magnitude()
        array([3.14159265, 3.14159265, 3.14159265, 0.        ])

        Magnitude of a single rotation:

        >>> r[0].magnitude()
        3.141592653589793
        """
        cdef double[:, :] quat = self._quat
        cdef Py_ssize_t num_rotations = quat.shape[0]
        cdef double[:] angles = _empty1(num_rotations)

        for ind in range(num_rotations):
            angles[ind] = 2 * atan2(_norm3(quat[ind, :3]), abs(quat[ind, 3]))

        if self._single:
            return angles[0]
        else:
            return np.asarray(angles)

    def mean(self, weights=None):
        """Get the mean of the rotations.

        Parameters
        ----------
        weights : array_like shape (N,), optional
            Weights describing the relative importance of the rotations. If
            None (default), then all values in `weights` are assumed to be
            equal.

        Returns
        -------
        mean : `Rotation` instance
            Object containing the mean of the rotations in the current
            instance.

        Notes
        -----
        The mean used is the chordal L2 mean (also called the projected or
        induced arithmetic mean). If ``p`` is a set of rotations with mean
        ``m``, then ``m`` is the rotation which minimizes
        ``(weights[:, None, None] * (p.as_matrix() - m.as_matrix())**2).sum()``.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> r = R.from_euler('zyx', [[0, 0, 0],
        ...                          [1, 0, 0],
        ...                          [0, 1, 0],
        ...                          [0, 0, 1]], degrees=True)
        >>> r.mean().as_euler('zyx', degrees=True)
        array([0.24945696, 0.25054542, 0.24945696])
        """
        if weights is None:
            weights = np.ones(len(self))
        else:
            weights = np.asarray(weights)
            if weights.ndim != 1:
                raise ValueError("Expected `weights` to be 1 dimensional, got "
                                 "shape {}.".format(weights.shape))
            if weights.shape[0] != len(self):
                raise ValueError("Expected `weights` to have number of values "
                                 "equal to number of rotations, got "
                                 "{} values and {} rotations.".format(
                                    weights.shape[0], len(self)))
            if np.any(weights < 0):
                raise ValueError("`weights` must be non-negative.")

        quat = np.asarray(self._quat)
        K = np.dot(weights * quat.T, quat)
        l, v = np.linalg.eigh(K)
        return self.__class__(v[:, -1], normalize=False)

    def reduce(self, left=None, right=None, return_indices=False):
        """Reduce this rotation with the provided rotation groups.

        Reduction of a rotation ``p`` is a transformation of the form
        ``q = l * p * r``, where ``l`` and ``r`` are chosen from `left` and
        `right` respectively, such that rotation ``q`` has the smallest
        magnitude.

        If `left` and `right` are rotation groups representing symmetries of
        two objects rotated by ``p``, then ``q`` is the rotation of the
        smallest magnitude to align these objects considering their symmetries.

        Parameters
        ----------
        left : `Rotation` instance, optional
            Object containing the left rotation(s). Default value (None)
            corresponds to the identity rotation.
        right : `Rotation` instance, optional
            Object containing the right rotation(s). Default value (None)
            corresponds to the identity rotation.
        return_indices : bool, optional
            Whether to return the indices of the rotations from `left` and
            `right` used for reduction.

        Returns
        -------
        reduced : `Rotation` instance
            Object containing reduced rotations.
        left_best, right_best: integer ndarray
            Indices of elements from `left` and `right` used for reduction.
        """
        if left is None and right is None:
            reduced = self.__class__(self._quat, normalize=False, copy=True)
            if return_indices:
                return reduced, None, None
            else:
                return reduced
        elif right is None:
            right = Rotation.identity()
        elif left is None:
            left = Rotation.identity()

        # Levi-Civita tensor for triple product computations
        e = np.zeros((3, 3, 3))
        e[0, 1, 2] = e[1, 2, 0] = e[2, 0, 1] = 1
        e[0, 2, 1] = e[2, 1, 0] = e[1, 0, 2] = -1

        # We want to calculate the real components of q = l * p * r. It can
        # be shown that:
        #     qs = ls * ps * rs - ls * dot(pv, rv) - ps * dot(lv, rv)
        #          - rs * dot(lv, pv) - dot(cross(lv, pv), rv)
        # where ls and lv denote the scalar and vector components of l.

        def split_rotation(R):
            q = np.atleast_2d(R.as_quat())
            return q[:, -1], q[:, :-1]

        p = self
        ps, pv = split_rotation(p)
        ls, lv = split_rotation(left)
        rs, rv = split_rotation(right)

        qs = np.abs(np.einsum('i,j,k', ls, ps, rs) -
                    np.einsum('i,jx,kx', ls, pv, rv) -
                    np.einsum('ix,j,kx', lv, ps, rv) -
                    np.einsum('ix,jx,k', lv, pv, rs) -
                    np.einsum('xyz,ix,jy,kz', e, lv, pv, rv))
        qs = np.reshape(np.rollaxis(qs, 1), (qs.shape[1], -1))

        # Find best indices from scalar components
        max_ind = np.argmax(np.reshape(qs, (len(qs), -1)), axis=1)
        left_best = max_ind // len(rv)
        right_best = max_ind % len(rv)

        if not left.single:
            left = left[left_best]
        if not right.single:
            right = right[right_best]

        # Reduce the rotation using the best indices
        reduced = left * p * right
        if self._single:
            # Reduce the rotation using the best indices
            reduced = self.__class__(reduced.as_quat()[0], normalize=False)
            left_best = left_best[0]
            right_best = right_best[0]

        if return_indices:
            if left is None:
                left_best = None
            if right is None:
                right_best = None
            return reduced, left_best, right_best
        else:
            return reduced

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

    def __getitem__(self, indexer):
        """Extract rotation(s) at given index(es) from object.

        Create a new `Rotation` instance containing a subset of rotations
        stored in this object.

        Parameters
        ----------
        indexer : index, slice, or index array
            Specifies which rotation(s) to extract. A single indexer must be
            specified, i.e. as if indexing a 1 dimensional array or list.

        Returns
        -------
        rotation : `Rotation` instance
            Contains
                - a single rotation, if `indexer` is a single index
                - a stack of rotation(s), if `indexer` is a slice, or and index
                  array.

        Raises
        ------
        TypeError if the instance was created as a single rotation.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> r = R.from_quat([
        ... [1, 1, 0, 0],
        ... [0, 1, 0, 1],
        ... [1, 1, -1, 0]])
        >>> r.as_quat()
        array([[ 0.70710678,  0.70710678,  0.        ,  0.        ],
               [ 0.        ,  0.70710678,  0.        ,  0.70710678],
               [ 0.57735027,  0.57735027, -0.57735027,  0.        ]])

        Indexing using a single index:

        >>> p = r[0]
        >>> p.as_quat()
        array([0.70710678, 0.70710678, 0.        , 0.        ])

        Array slicing:

        >>> q = r[1:3]
        >>> q.as_quat()
        array([[ 0.        ,  0.70710678,  0.        ,  0.70710678],
               [ 0.57735027,  0.57735027, -0.57735027,  0.        ]])

        """
        if self._single:
            raise TypeError("Single rotation is not subscriptable.")

        return self.__class__(np.asarray(self._quat)[indexer], normalize=False)

    @classmethod
    def identity(cls, num=None):
        """Get identity rotation(s).

        Composition with the identity rotation has no effect.

        Parameters
        ----------
        num : int or None, optional
            Number of identity rotations to generate. If None (default), then a
            single rotation is generated.

        Returns
        -------
        identity : Rotation object
            The identity rotation.
        """
        if num is None:
            q = [0, 0, 0, 1]
        else:
            q = np.zeros((num, 4))
            q[:, 3] = 1
        return cls(q, normalize=False)

    @classmethod
    def random(cls, num=None, random_state=None):
        """Generate uniformly distributed rotations.

        Parameters
        ----------
        num : int or None, optional
            Number of random rotations to generate. If None (default), then a
            single rotation is generated.
        random_state : int, RandomState instance or None, optional
            Accepts an integer as a seed for the random generator or a
            RandomState object. If None (default), uses global `numpy.random`
            random state.

        Returns
        -------
        random_rotation : `Rotation` instance
            Contains a single rotation if `num` is None. Otherwise contains a
            stack of `num` rotations.

        Notes
        -----
        This function is optimized for efficiently sampling random rotation
        matrices in three dimensions. For generating random rotation matrices
        in higher dimensions, see `scipy.stats.special_ortho_group`.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Sample a single rotation:

        >>> R.random(random_state=1234).as_euler('zxy', degrees=True)
        array([-110.5976185 ,   55.32758512,   76.3289269 ])

        Sample a stack of rotations:

        >>> R.random(5, random_state=1234).as_euler('zxy', degrees=True)
        array([[-110.5976185 ,   55.32758512,   76.3289269 ],
               [ -91.59132005,  -14.3629884 ,  -93.91933182],
               [  25.23835501,   45.02035145, -121.67867086],
               [ -51.51414184,  -15.29022692, -172.46870023],
               [ -81.63376847,  -27.39521579,    2.60408416]])

        See Also
        --------
        scipy.stats.special_ortho_group

       """
        random_state = check_random_state(random_state)

        if num is None:
            sample = random_state.normal(size=4)
        else:
            sample = random_state.normal(size=(num, 4))

        return cls(sample)

    @classmethod
    def align_vectors(cls, a, b, weights=None, return_sensitivity=False):
        """Estimate a rotation to optimally align two sets of vectors.

        Find a rotation between frames A and B which best aligns a set of
        vectors `a` and `b` observed in these frames. The following loss
        function is minimized to solve for the rotation matrix
        :math:`C`:

        .. math::

            L(C) = \\frac{1}{2} \\sum_{i = 1}^{n} w_i \\lVert \\mathbf{a}_i -
            C \\mathbf{b}_i \\rVert^2 ,

        where :math:`w_i`'s are the `weights` corresponding to each vector.

        The rotation is estimated with Kabsch algorithm [1]_.

        Parameters
        ----------
        a : array_like, shape (N, 3)
            Vector components observed in initial frame A. Each row of `a`
            denotes a vector.
        b : array_like, shape (N, 3)
            Vector components observed in another frame B. Each row of `b`
            denotes a vector.
        weights : array_like shape (N,), optional
            Weights describing the relative importance of the vector
            observations. If None (default), then all values in `weights` are
            assumed to be 1.
        return_sensitivity : bool, optional
            Whether to return the sensitivity matrix. See Notes for details.
            Default is False.

        Returns
        -------
        estimated_rotation : `Rotation` instance
            Best estimate of the rotation that transforms `b` to `a`.
        rmsd : float
            Root mean square distance (weighted) between the given set of
            vectors after alignment. It is equal to ``sqrt(2 * minimum_loss)``,
            where ``minimum_loss`` is the loss function evaluated for the
            found optimal rotation.
        sensitivity_matrix : ndarray, shape (3, 3)
            Sensitivity matrix of the estimated rotation estimate as explained
            in Notes. Returned only when `return_sensitivity` is True.

        Notes
        -----
        This method can also compute the sensitivity of the estimated rotation
        to small perturbations of the vector measurements. Specifically we
        consider the rotation estimate error as a small rotation vector of
        frame A. The sensitivity matrix is proportional to the covariance of
        this rotation vector assuming that the vectors in `a` was measured with
        errors significantly less than their lengths. To get the true
        covariance matrix, the returned sensitivity matrix must be multiplied
        by harmonic mean [3]_ of variance in each observation. Note that
        `weights` are supposed to be inversely proportional to the observation
        variances to get consistent results. For example, if all vectors are
        measured with the same accuracy of 0.01 (`weights` must be all equal),
        then you should multiple the sensitivity matrix by 0.01**2 to get the
        covariance.

        Refer to [2]_ for more rigorous discussion of the covariance
        estimation.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Kabsch_algorithm
        .. [2] F. Landis Markley,
                "Attitude determination using vector observations: a fast
                optimal matrix algorithm", Journal of Astronautical Sciences,
                Vol. 41, No.2, 1993, pp. 261-280.
        .. [3] https://en.wikipedia.org/wiki/Harmonic_mean
        """
        a = np.asarray(a)
        if a.ndim != 2 or a.shape[-1] != 3:
            raise ValueError("Expected input `a` to have shape (N, 3), "
                             "got {}".format(a.shape))
        b = np.asarray(b)
        if b.ndim != 2 or b.shape[-1] != 3:
            raise ValueError("Expected input `b` to have shape (N, 3), "
                             "got {}.".format(b.shape))

        if a.shape != b.shape:
            raise ValueError("Expected inputs `a` and `b` to have same shapes"
                             ", got {} and {} respectively.".format(
                                a.shape, b.shape))

        if weights is None:
            weights = np.ones(len(b))
        else:
            weights = np.asarray(weights)
            if weights.ndim != 1:
                raise ValueError("Expected `weights` to be 1 dimensional, got "
                                 "shape {}.".format(weights.shape))
            if weights.shape[0] != b.shape[0]:
                raise ValueError("Expected `weights` to have number of values "
                                 "equal to number of input vectors, got "
                                 "{} values and {} vectors.".format(
                                    weights.shape[0], b.shape[0]))

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

        rmsd = np.sqrt(max(
            np.sum(weights * np.sum(b ** 2 + a ** 2, axis=1)) - 2 * np.sum(s),
            0))

        if return_sensitivity:
            zeta = (s[0] + s[1]) * (s[1] + s[2]) * (s[2] + s[0])
            kappa = s[0] * s[1] + s[1] * s[2] + s[2] * s[0]
            with np.errstate(divide='ignore', invalid='ignore'):
                sensitivity = np.mean(weights) / zeta * (
                        kappa * np.eye(3) + np.dot(B, B.T))
            return cls.from_matrix(C), rmsd, sensitivity
        else:
            return cls.from_matrix(C), rmsd


class Slerp(object):
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
        if rotations.single:
            raise ValueError("`rotations` must be a sequence of rotations.")

        if len(rotations) == 1:
                raise ValueError("`rotations` must contain at least 2 "
                                 "rotations.")

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
