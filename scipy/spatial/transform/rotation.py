from __future__ import division, print_function, absolute_import

import re
import warnings
import numpy as np
import scipy.linalg


_AXIS_TO_IND = {'x': 0, 'y': 1, 'z': 2}


def _elementary_basis_vector(axis):
    b = np.zeros(3)
    b[_AXIS_TO_IND[axis]] = 1
    return b


def _compute_euler_from_dcm(dcm, seq, extrinsic=False):
    # The algorithm assumes intrinsic frame transformations. For representation
    # the paper uses transformation matrices, which are transpose of the
    # direction cosine matrices used by our Rotation class.
    # Adapt the algorithm for our case by
    # 1. Instead of transposing our representation, use the transpose of the
    #    O matrix as defined in the paper, and be careful to swap indices
    # 2. Reversing both axis sequence and angles for extrinsic rotations

    if extrinsic:
        seq = seq[::-1]

    if dcm.ndim == 2:
        dcm = dcm[None, :, :]
    num_rotations = dcm.shape[0]

    # Step 0
    # Algorithm assumes axes as column vectors, here we use 1D vectors
    n1 = _elementary_basis_vector(seq[0])
    n2 = _elementary_basis_vector(seq[1])
    n3 = _elementary_basis_vector(seq[2])

    # Step 2
    sl = np.dot(np.cross(n1, n2), n3)
    cl = np.dot(n1, n3)

    # angle offset is lambda from the paper referenced in [2] from docstring of
    # `as_euler` function
    offset = np.arctan2(sl, cl)
    c = np.vstack((n2, np.cross(n1, n2), n1))

    # Step 3
    rot = np.array([
        [1, 0, 0],
        [0, cl, sl],
        [0, -sl, cl],
    ])
    res = np.einsum('...ij,...jk->...ik', c, dcm)
    dcm_transformed = np.einsum('...ij,...jk->...ik', res, c.T.dot(rot))

    # Step 4
    angles = np.empty((num_rotations, 3))
    # Ensure less than unit norm
    positive_unity = dcm_transformed[:, 2, 2] > 1
    negative_unity = dcm_transformed[:, 2, 2] < -1
    dcm_transformed[positive_unity, 2, 2] = 1
    dcm_transformed[negative_unity, 2, 2] = -1
    angles[:, 1] = np.arccos(dcm_transformed[:, 2, 2])

    # Steps 5, 6
    eps = 1e-7
    safe1 = (np.abs(angles[:, 1]) >= eps)
    safe2 = (np.abs(angles[:, 1] - np.pi) >= eps)

    # Step 4 (Completion)
    angles[:, 1] += offset

    # 5b
    safe_mask = np.logical_and(safe1, safe2)
    angles[safe_mask, 0] = np.arctan2(dcm_transformed[safe_mask, 0, 2],
                                      -dcm_transformed[safe_mask, 1, 2])
    angles[safe_mask, 2] = np.arctan2(dcm_transformed[safe_mask, 2, 0],
                                      dcm_transformed[safe_mask, 2, 1])

    if extrinsic:
        # For extrinsic, set first angle to zero so that after reversal we
        # ensure that third angle is zero
        # 6a
        angles[~safe_mask, 0] = 0
        # 6b
        angles[~safe1, 2] = np.arctan2(
            dcm_transformed[~safe1, 1, 0] - dcm_transformed[~safe1, 0, 1],
            dcm_transformed[~safe1, 0, 0] + dcm_transformed[~safe1, 1, 1]
        )
        # 6c
        angles[~safe2, 2] = -np.arctan2(
            dcm_transformed[~safe2, 1, 0] + dcm_transformed[~safe2, 0, 1],
            dcm_transformed[~safe2, 0, 0] - dcm_transformed[~safe2, 1, 1]
        )
    else:
        # For instrinsic, set third angle to zero
        # 6a
        angles[~safe_mask, 2] = 0
        # 6b
        angles[~safe1, 0] = np.arctan2(
            dcm_transformed[~safe1, 1, 0] - dcm_transformed[~safe1, 0, 1],
            dcm_transformed[~safe1, 0, 0] + dcm_transformed[~safe1, 1, 1]
        )
        # 6c
        angles[~safe2, 0] = np.arctan2(
            dcm_transformed[~safe2, 1, 0] + dcm_transformed[~safe2, 0, 1],
            dcm_transformed[~safe2, 0, 0] - dcm_transformed[~safe2, 1, 1]
        )

    # Step 7
    if seq[0] == seq[2]:
        # lambda = 0, so we can only ensure angle2 -> [0, pi]
        adjust_mask = np.logical_or(angles[:, 1] < 0, angles[:, 1] > np.pi)
    else:
        # lambda = + or - pi/2, so we can ensure angle2 -> [-pi/2, pi/2]
        adjust_mask = np.logical_or(angles[:, 1] < -np.pi / 2,
                                    angles[:, 1] > np.pi / 2)

    # Dont adjust gimbal locked angle sequences
    adjust_mask = np.logical_and(adjust_mask, safe_mask)

    angles[adjust_mask, 0] += np.pi
    angles[adjust_mask, 1] = 2 * offset - angles[adjust_mask, 1]
    angles[adjust_mask, 2] -= np.pi

    angles[angles < -np.pi] += 2 * np.pi
    angles[angles > np.pi] -= 2 * np.pi

    # Step 8
    if not np.all(safe_mask):
        warnings.warn("Gimbal lock detected. Setting third angle to zero since"
                      " it is not possible to uniquely determine all angles.")

    # Reverse role of extrinsic and intrinsic rotations, but let third angle be
    # zero for gimbal locked cases
    if extrinsic:
        angles = angles[:, ::-1]
    return angles


def _make_elementary_quat(axis, angles):
    quat = np.zeros((angles.shape[0], 4))

    quat[:, 3] = np.cos(angles / 2)
    quat[:, _AXIS_TO_IND[axis]] = np.sin(angles / 2)
    return quat


def _compose_quat(p, q):
    product = np.empty((max(p.shape[0], q.shape[0]), 4))
    product[:, 3] = p[:, 3] * q[:, 3] - np.sum(p[:, :3] * q[:, :3], axis=1)
    product[:, :3] = (p[:, None, 3] * q[:, :3] + q[:, None, 3] * p[:, :3] +
                      np.cross(p[:, :3], q[:, :3]))
    return product


def _elementary_quat_compose(seq, angles, intrinsic=False):
    result = _make_elementary_quat(seq[0], angles[:, 0])

    for idx, axis in enumerate(seq[1:], start=1):
        if intrinsic:
            result = _compose_quat(
                result,
                _make_elementary_quat(axis, angles[:, idx]))
        else:
            result = _compose_quat(
                _make_elementary_quat(axis, angles[:, idx]),
                result)
    return result


class Rotation(object):
    """Rotation in 3 dimensions.

    This class provides an interface to initialize from and represent rotations
    with:

        - Quaternions
        - Direction Cosine Matrices
        - Rotation Vectors
        - Euler angles

    The following operations on rotations are supported:

        - Application on vectors
        - Rotation Composition
        - Rotation Inversion
        - Rotation Indexing

    Indexing within a rotation is supported since multiple rotation transforms
    can be stored within a single `Rotation` instance.

    To create `Rotation` objects use `from_...` classmethods, `__init__` is not
    supposed to be used directly.

    Methods
    -------
    __len__
    from_quat
    from_dcm
    from_rotvec
    from_euler
    as_quat
    as_dcm
    as_rotvec
    as_euler
    apply
    __mul__
    inv
    __getitem__

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

    >>> r.as_dcm()
    array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
    [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
    [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
    >>> r.as_rotvec()
    array([0.        , 0.        , 1.57079633])
    >>> r.as_euler('zyx', degrees=True)
    array([90.,  0.,  0.])

    The same rotation can be initialized using a direction cosine matrix:

    >>> r = R.from_dcm(np.array([
    ... [0, -1, 0],
    ... [1, 0, 0],
    ... [0, 0, 1]]))

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
    >>> r.as_dcm()
    array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
           [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
           [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
    >>> r.as_euler('zyx', degrees=True)
    array([90.,  0.,  0.])

    The `from_euler` function is quite flexible in the range of input formats
    it supports. Here we initialize a single rotation about a single axis:

    >>> r = R.from_euler('z', 90, degrees=True)

    Again, the object is representation independent and can be converted to any
    other format:

    >>> r.as_quat()
    array([0.        , 0.        , 0.70710678, 0.70710678])
    >>> r.as_dcm()
    array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
           [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
           [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
    >>> r.as_rotvec()
    array([0.        , 0.        , 1.57079633])

    It is also possible to initialize multiple rotations in a single instance
    using any of the `from_...` functions. Here we initialize a stack of 3
    rotations using the `from_euler` function:

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
    >>> p.as_dcm()
    array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
           [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
           [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
    >>> q = r[1:3]
    >>> q.as_quat()
    array([[0.        , 0.38268343, 0.        , 0.92387953],
           [0.39190384, 0.36042341, 0.43967974, 0.72331741]])

    Multiple rotations can be composed using the `*` operator:

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
    >>> r2 = r.inv()
    >>> r2.as_euler('zyx', degrees=True)
    array([[-90.,   0.,   0.],
           [-45.,   0.,   0.]])

    These examples serve as an overview into the `Rotation` class and highlight
    major functionalities. For more thorough examples of the range of input and
    output formats supported, consult the individual method's examples.
    """
    def __init__(self, quat, normalized=False, copy=True):
        self._single = False
        quat = np.asarray(quat, dtype=float)

        if quat.ndim not in [1, 2] or quat.shape[-1] != 4:
            raise ValueError("Expected `quat` to have shape (4,) or (N x 4), "
                             "got {}.".format(quat.shape))

        # If a single quaternion is given, convert it to a 2D 1 x 4 matrix but
        # set self._single to True so that we can return appropriate objects
        # in the `to_...` methods
        if quat.shape == (4,):
            quat = quat[None, :]
            self._single = True

        if normalized:
            self._quat = quat.copy() if copy else quat
        else:
            self._quat = quat.copy()
            norms = scipy.linalg.norm(quat, axis=1)

            zero_norms = norms == 0
            if zero_norms.any():
                raise ValueError("Found zero norm quaternions in `quat`.")

            # Ensure norm is broadcasted along each column.
            self._quat[~zero_norms] /= norms[~zero_norms][:, None]

    def __len__(self):
        """Number of rotations contained in this object.

        Multiple rotations can be stored in a single instance.

        Returns
        -------
        length : int
            Number of rotations stored in object.
        """
        return self._quat.shape[0]

    @classmethod
    def from_quat(cls, quat, normalized=False):
        """Initialize from quaternions.

        3D rotations can be represented using unit-norm quaternions [1]_.

        Parameters
        ----------
        quat : array_like, shape (N, 4) or (4,)
            Each row is a (possibly non-unit norm) quaternion in scalar-last
            (x, y, z, w) format.
        normalized : boolean, optional
            If `False`, input quaternions are normalized to unit norm before
            being stored. If `True`, quaternions are assumed to already have
            unit norm and are stored as given. Default is `False`.

        Returns
        -------
        rotation : `Rotation` instance
            Object containing the rotations represented by input quaternions.

        References
        ----------
        .. [1] `Quaternions and Spatial Rotation
               <https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation>`_

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

        By default, quaternions are normalized before initialization.

        >>> r = R.from_quat([0, 0, 1, 1])
        >>> r.as_quat()
        array([0.        , 0.        , 0.70710678, 0.70710678])

        If unit norms are ensured, skip the normalization step.

        >>> r = R.from_quat([0, 0, 1, 0], normalized=True)
        >>> r.as_quat()
        array([0., 0., 1., 0.])
        """

        return cls(quat, normalized)

    @classmethod
    def from_dcm(cls, dcm):
        """Initialize from direction cosine matrices.

        Rotations in 3 dimensions can be represented using 3 x 3 proper
        orthogonal matrices [1]_. If the input is not proper orthogonal,
        an approximation is created using the method described in [2]_.

        Parameters
        ----------
        dcm : array_like, shape (N, 3, 3) or (3, 3)
            A single matrix or a stack of matrices, where `dcm[i]` is the i-th
            matrix.

        Returns
        -------
        rotation : `Rotation` instance
            Object containing the rotations represented by the input direction
            cosine matrices.

        References
        ----------
        .. [1] `Direction Cosine Matrix
                <https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions>`_
        .. [2] F. Landis Markley, `Unit Quaternion from Rotation Matrix
               <https://arc.aiaa.org/doi/abs/10.2514/1.31730>`_

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Initialize a single rotation:

        >>> r = R.from_dcm([
        ... [0, -1, 0],
        ... [1, 0, 0],
        ... [0, 0, 1]])
        >>> r.as_dcm().shape
        (3, 3)

        Initialize multiple rotations in a single object:

        >>> r = R.from_dcm([
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
        >>> r.as_dcm().shape
        (2, 3, 3)

        If input matrices are not special orthogonal (orthogonal with
        determinant equal to +1), then a special orthogonal estimate is stored:

        >>> a = np.array([
        ... [0, -0.5, 0],
        ... [0.5, 0, 0],
        ... [0, 0, 0.5]])
        >>> np.linalg.det(a)
        0.12500000000000003
        >>> dcm = r.as_dcm()
        >>> dcm
        array([[-0.38461538, -0.92307692,  0.        ],
               [ 0.92307692, -0.38461538,  0.        ],
               [ 0.        ,  0.        ,  1.        ]])
        >>> np.linalg.det(dcm)
        1.0000000000000002

        It is also possible to have a stack containing a single rotation:

        >>> r = R.from_dcm([[
        ... [0, -1, 0],
        ... [1, 0, 0],
        ... [0, 0, 1]]])
        >>> r.as_dcm()
        array([[[ 0., -1.,  0.],
                [ 1.,  0.,  0.],
                [ 0.,  0.,  1.]]])
        >>> r.as_dcm().shape
        (1, 3, 3)
        """
        is_single = False
        dcm = np.asarray(dcm, dtype=float)

        if dcm.ndim not in [2, 3] or dcm.shape[-2:] != (3, 3):
            raise ValueError("Expected `dcm` to have shape (3, 3) or "
                             "(N, 3, 3), got {}".format(dcm.shape))

        # If a single dcm is given, convert it to 3D 1 x 3 x 3 matrix but set
        # self._single to True so that we can return appropriate objects in
        # the `to_...` methods
        if dcm.shape == (3, 3):
            dcm = dcm.reshape((1, 3, 3))
            is_single = True

        num_rotations = dcm.shape[0]

        decision_matrix = np.empty((num_rotations, 4))
        decision_matrix[:, :3] = dcm.diagonal(axis1=1, axis2=2)
        decision_matrix[:, -1] = decision_matrix[:, :3].sum(axis=1)
        choices = decision_matrix.argmax(axis=1)

        quat = np.empty((num_rotations, 4))

        ind = np.nonzero(choices != 3)[0]
        i = choices[ind]
        j = (i + 1) % 3
        k = (j + 1) % 3

        quat[ind, i] = 1 - decision_matrix[ind, -1] + 2 * dcm[ind, i, i]
        quat[ind, j] = dcm[ind, j, i] + dcm[ind, i, j]
        quat[ind, k] = dcm[ind, k, i] + dcm[ind, i, k]
        quat[ind, 3] = dcm[ind, k, j] - dcm[ind, j, k]

        ind = np.nonzero(choices == 3)[0]
        quat[ind, 0] = dcm[ind, 2, 1] - dcm[ind, 1, 2]
        quat[ind, 1] = dcm[ind, 0, 2] - dcm[ind, 2, 0]
        quat[ind, 2] = dcm[ind, 1, 0] - dcm[ind, 0, 1]
        quat[ind, 3] = 1 + decision_matrix[ind, -1]

        quat /= np.linalg.norm(quat, axis=1)[:, None]

        if is_single:
            return cls(quat[0], normalized=True, copy=False)
        else:
            return cls(quat, normalized=True, copy=False)

    @classmethod
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
        .. [1] `Rotation Vectors
                <https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation#Rotation_vector>`_

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

        if rotvec.ndim not in [1, 2] or rotvec.shape[-1] != 3:
            raise ValueError("Expected `rot_vec` to have shape (3,) "
                             "or (N, 3), got {}".format(rotvec.shape))

        # If a single vector is given, convert it to a 2D 1 x 3 matrix but
        # set self._single to True so that we can return appropriate objects
        # in the `as_...` methods
        if rotvec.shape == (3,):
            rotvec = rotvec[None, :]
            is_single = True

        num_rotations = rotvec.shape[0]

        norms = np.linalg.norm(rotvec, axis=1)
        small_angle = (norms <= 1e-3)
        large_angle = ~small_angle

        scale = np.empty(num_rotations)
        scale[small_angle] = (0.5 - norms[small_angle] ** 2 / 48 +
                              norms[small_angle] ** 4 / 3840)
        scale[large_angle] = (np.sin(norms[large_angle] / 2) /
                              norms[large_angle])

        quat = np.empty((num_rotations, 4))
        quat[:, :3] = scale[:, None] * rotvec
        quat[:, 3] = np.cos(norms / 2)

        if is_single:
            return cls(quat[0], normalized=True, copy=False)
        else:
            return cls(quat, normalized=True, copy=False)

    @classmethod
    def from_euler(cls, seq, angles, degrees=False):
        """Initialize from Euler angles.

        Rotations in 3 dimensions can be represented by a sequece of 3
        rotations around a sequence of axes. In theory, any three axes spanning
        the 3D Euclidean space are enough. In practice the axes of rotation are
        chosen to be the basis vectors.

        The three rotations can either be in a global frame of reference
        (extrinsic) or in a body centred frame of refernce (intrinsic), which
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

        degrees : boolean, optional
            If True, then the given angles are assumed to be in degrees.
            Default is False.

        Returns
        -------
        rotation : `Rotation` instance
            Object containing the rotation represented by the sequence of
            rotations around given axes with given angles.

        References
        ----------
        .. [1] `Euler angle definitions
                <https://en.wikipedia.org/wiki/Euler_angles#Definition_by_intrinsic_rotations>`_

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
                                 "of axes specified, got {} for shape").format(
                                 angles.shape)

            if angles.ndim == 1:
                # (1, num_axes)
                angles = angles[None, :]
                is_single = True

        # By now angles should have shape (num_rot, num_axes)
        # sanity check
        if angles.ndim != 2 or angles.shape[-1] != num_axes:
            raise ValueError("Expected angles to have shape (num_rotations, "
                             "num_axes), got {}.".format(angles.shape))

        quat = _elementary_quat_compose(seq, angles, intrinsic)
        return cls(quat[0] if is_single else quat, normalized=True, copy=False)

    def as_quat(self):
        """Represent as quaternions.

        Rotations in 3 dimensions can be represented using unit norm
        quaternions [1]_. The mapping from quaternions to rotations is
        two-to-one, i.e. quaternions `q` and `-q`, where `-q` simply reverses
        the sign of each component, represent the same spatial rotation.

        Returns
        -------
        quat : `numpy.ndarray`, shape (4,) or (N, 4)
            Shape depends on shape of inputs used for initialization.

        References
        ----------
        .. [1] `Quaternions and Spatial Rotation
               <https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation>`_

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Represent a single rotation:

        >>> r = R.from_dcm([
        ... [0, -1, 0],
        ... [1, 0, 0],
        ... [0, 0, 1]])
        >>> r.as_quat()
        array([0.        , 0.        , 0.70710678, 0.70710678])
        >>> r.as_quat().shape
        (4,)

        Represent a stack with a single rotation:

        >>> r = R.from_quat([[0, 0, 0, 1]])
        >>> r.as_quat().shape
        (1, 4)

        Represent multiple rotaions in a single object:

        >>> r = R.from_rotvec([[np.pi, 0, 0], [0, 0, np.pi/2]])
        >>> r.as_quat().shape
        (2, 4)
        """
        if self._single:
            return self._quat[0].copy()
        else:
            return self._quat.copy()

    def as_dcm(self):
        """Represent as direction cosine matrices.

        3D rotations can be represented using direction cosine matrices, which
        are 3 x 3 real orthogonal matrices with determinant equal to +1 [1]_.

        Returns
        -------
        dcm : `numpy.ndarray`, shape (3, 3) or (N, 3, 3)
            Shape depends on shape of inputs used for initialization.

        References
        ----------
        .. [1] `Direction Cosine Matrix
                <https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions>`_

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Represent a single rotation:

        >>> r = R.from_rotvec([0, 0, np.pi/2])
        >>> r.as_dcm()
        array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
               [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
               [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
        >>> r.as_dcm().shape
        (3, 3)

        Represent a stack with a single rotation:

        >>> r = R.from_quat([[1, 1, 0, 0]])
        >>> r.as_dcm()
        array([[[ 0.,  1.,  0.],
                [ 1.,  0.,  0.],
                [ 0.,  0., -1.]]])
        >>> r.as_dcm().shape
        (1, 3, 3)

        Represent multiple rotations:

        >>> r = R.from_rotvec([[np.pi/2, 0, 0], [0, 0, np.pi/2]])
        >>> r.as_dcm()
        array([[[ 1.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                [ 0.00000000e+00,  2.22044605e-16, -1.00000000e+00],
                [ 0.00000000e+00,  1.00000000e+00,  2.22044605e-16]],
               [[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
                [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
                [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]]])
        >>> r.as_dcm().shape
        (2, 3, 3)
        """

        x = self._quat[:, 0]
        y = self._quat[:, 1]
        z = self._quat[:, 2]
        w = self._quat[:, 3]

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

        num_rotations = len(self)
        dcm = np.empty((num_rotations, 3, 3))

        dcm[:, 0, 0] = x2 - y2 - z2 + w2
        dcm[:, 1, 0] = 2 * (xy + zw)
        dcm[:, 2, 0] = 2 * (xz - yw)

        dcm[:, 0, 1] = 2 * (xy - zw)
        dcm[:, 1, 1] = - x2 + y2 - z2 + w2
        dcm[:, 2, 1] = 2 * (yz + xw)

        dcm[:, 0, 2] = 2 * (xz + yw)
        dcm[:, 1, 2] = 2 * (yz - xw)
        dcm[:, 2, 2] = - x2 - y2 + z2 + w2

        if self._single:
            return dcm[0]
        else:
            return dcm

    def as_rotvec(self):
        """Represent as rotation vectors.

        A rotation vector is a 3 dimensional vector which is co-directional to
        the axis of rotation and whose norm gives the angle of rotation (in
        radians) [1]_.

        Returns
        -------
        rotvec : `numpy.ndarray`, shape (3,) or (N, 3)
            Shape depends on shape of inputs used for initialization.

        References
        ----------
        .. [1] `Rotation Vectors
                <https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation#Rotation_vector>`_

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
        quat = self._quat.copy()
        # w > 0 to ensure 0 <= angle <= pi
        quat[quat[:, 3] < 0] *= -1

        angle = 2 * np.arctan2(np.linalg.norm(quat[:, :3], axis=1), quat[:, 3])

        small_angle = (angle <= 1e-3)
        large_angle = ~small_angle

        num_rotations = len(self)
        scale = np.empty(num_rotations)
        scale[small_angle] = (2 + angle[small_angle] ** 2 / 12 +
                              7 * angle[small_angle] ** 4 / 2880)
        scale[large_angle] = (angle[large_angle] /
                              np.sin(angle[large_angle] / 2))

        rotvec = scale[:, None] * quat[:, :3]

        if self._single:
            return rotvec[0]
        else:
            return rotvec

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
        angles : `numpy.ndarray`, shape (3,) or (N, 3)
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
        .. [1] `Euler angle definitions
                <https://en.wikipedia.org/wiki/Euler_angles#Definition_by_intrinsic_rotations>`_
        .. [2] Malcolm D. Shuster, F. Landis Markley
                `General Formula for Euler Angles
                <https://arc.aiaa.org/doi/abs/10.2514/1.16622>`_
        .. [3] `Gimbal lock
                <https://en.wikipedia.org/wiki/Gimbal_lock#In_applied_mathematics>`_

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

        angles = _compute_euler_from_dcm(self.as_dcm(), seq, extrinsic)
        if degrees:
            angles = np.rad2deg(angles)

        return angles[0] if self._single else angles

    def apply(self, vectors, inverse=False):
        """Apply this rotation to a set of vectors.

        If the original frame rotates to the final frame by this rotation, then
        its application to a vector can be seen in two ways:

            - As a projection of vector components expressed in the final frame
              to the original frame.
            - As the physical rotation of a vector being glued to the original
              frame as it rotates. In this case the vector components are
              expressed in the original frame before and after the rotation.

        In terms of DCMs, this application is the same as
        `self.as_dcm().dot(vectors)`.

        Parameters
        ----------
        vectors : array_like, shape (3,) or (N, 3)
            Each `vectors[i]` represents a vector in 3D space. A single vector
            can either be specified with shape `(3, )` or `(1, 3)`. The number
            of rotations and number of vectors given must follow standard numpy
            broadcasting rules: either one of them equals unity or they both
            equal each other.
        inverse : boolean, optional
            If `inverse` is `True` then the inverse of the rotation(s) is
            applied to the input vectors. Default is `False`.

        Returns
        -------
        rotated_vectors : `numpy.ndarray`, shape (3,) or (N, 3)
            Result of applying rotation on input vectors.
            Shape depends on the following cases:

                - If object contains a single rotation (as opposed to a stack
                  with a single rotation) and a single vector is specified with
                  shape `(3,)`, then `output` has shape `(3,)`.
                - In all other cases, `output` has shape `(N, 3)`, where `N` is
                  either the number of rotations or vectors.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Single rotation applied on a single vector:

        >>> vector = np.array([1, 0, 0])
        >>> r = R.from_rotvec([0, 0, np.pi/2])
        >>> r.as_dcm()
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
        >>> r.as_dcm()
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
        >>> r.as_dcm()
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

        dcm = self.as_dcm()
        if self._single:
            dcm = dcm[None, :, :]

        n_vectors = vectors.shape[0]
        n_rotations = len(self)

        if n_vectors != 1 and n_rotations != 1 and n_vectors != n_rotations:
            raise ValueError("Expected equal numbers of rotations and vectors "
                             ", or a single rotation, or a single vector, got "
                             "{} rotations and {} vectors.".format(
                                n_rotations, n_vectors))

        if inverse:
            result = np.einsum('ikj,ik->ij', dcm, vectors)
        else:
            result = np.einsum('ijk,ik->ij', dcm, vectors)

        if self._single and single_vector:
            return result[0]
        else:
            return result

    def __mul__(self, other):
        """Compose this rotation with the other.

        If `p` and `q` are two rotations, then the composition of 'q followed
        by p' is equivalent to `p * q`. In terms of DCMs, the composition can
        be expressed as `p.as_dcm().dot(q.as_dcm())`.

        Parameters
        ----------
        other : `Rotation` instance
            Object containing the rotaions to be composed with this one. Note
            that rotation compositions are not commutative, so `p * q` is
            different from `q * p`.

        Returns
        -------
        composition : `Rotation` instance
            This function supports composition of multiple rotations at a time.
            The following cases are possible:

            - Either `p` or `q` contains a single rotation. In this case
              `output` contains the result of composing each rotation in the
              other object with the single rotation.
            - Both `p` and `q` contain `N` rotations. In this case each
              rotation `p[i]` is composed with the corresponding rotation
              `q[i]` and `output` contains `N` rotations.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Composition of two single rotations:

        >>> p = R.from_quat([0, 0, 1, 1])
        >>> q = R.from_quat([1, 0, 0, 1])
        >>> p.as_dcm()
        array([[ 0., -1.,  0.],
               [ 1.,  0.,  0.],
               [ 0.,  0.,  1.]])
        >>> q.as_dcm()
        array([[ 1.,  0.,  0.],
               [ 0.,  0., -1.],
               [ 0.,  1.,  0.]])
        >>> r = p * q
        >>> r.as_dcm()
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
        if not(len(self) == 1 or len(other) == 1 or len(self) == len(other)):
            raise ValueError("Expected equal number of rotations in both "
                             "or a single rotation in either object, "
                             "got {} rotations in first and {} rotations in "
                             "second object.".format(
                                len(self), len(other)))
        result = _compose_quat(self._quat, other._quat)
        if self._single and other._single:
            result = result[0]
        return self.__class__(result, normalized=True, copy=False)

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
        quat = self._quat.copy()
        quat[:, -1] *= -1
        if self._single:
            quat = quat[0]
        return self.__class__(quat, normalized=True, copy=False)

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
        return self.__class__(self._quat[indexer], normalized=True)
