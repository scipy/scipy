from __future__ import division, print_function, absolute_import

import numpy as np
import scipy.linalg
import re
import warnings


AXIS_TO_IND = {'x': 0, 'y': 1, 'z': 2}


def _elementary_basis_vector(axis):
    b = np.zeros(3)
    b[AXIS_TO_IND[axis]] = 1
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
    num_rotations = angles.shape[0]
    quat = np.zeros((num_rotations, 4))

    quat[:, 3] = np.cos(angles / 2)
    quat[:, AXIS_TO_IND[axis]] = np.sin(angles / 2)
    return quat


def _compose_quat(p, q):
    product = np.empty((max(p.shape[0], q.shape[0]), 4))
    # Scalar part of result
    product[:, 3] = p[:, 3] * q[:, 3] - np.sum(p[:, :3] * q[:, :3], axis=1)
    # Vector part of result
    product[:, :3] = (p[:, None, 3] * q[:, :3] + q[:, None, 3] * p[:, :3] +
                      np.cross(p[:, :3], q[:, :3]))
    return product


def _elementary_quat_compose(seq, angles, intrinsic=False):
    # Initialize result to first axis
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

    This class will include initializers from different representations,
    converters and some useful algorithms such as quaternion slerp and
    rotation estimation.

    For initializing Rotations usage of `from_...` methods such as
    `from_quaternion` is recommended instead of using `__init__`.

    Methods
    -------
    from_quaternion
    as_quaternion
    from_dcm
    as_dcm
    from_rotvec
    as_rotvec
    from_euler
    as_euler
    inv
    __mul__
    apply
    __getitem__
    """
    def __init__(self, quat, normalized=False, copy=True):
        self._single = False
        # Try to convert to numpy array
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
        if not normalized:
            self._quat = quat.copy()
            # L2 norm of each row
            norms = scipy.linalg.norm(quat, axis=1)

            # Raise ValueError for zero (eps?) norm quaternions
            zero_norms = norms == 0
            if zero_norms.any():
                raise ValueError("Found zero norm quaternions in `quat`.")

            # Normalize each quaternion, ensuring norm is broadcasted along
            # each column.
            self._quat[~zero_norms] /= norms[~zero_norms][:, None]

    @classmethod
    def from_quaternion(cls, quat, normalized=False):
        """Initialize Rotation from quaternions.

        This classmethod returns a `Rotation` object from the input quaternions
        If `normalized` is `True`, then the quaternions are assumed to have
        unit norm, else the quaternions are normalized before the object is
        created.

        Parameters
        ----------
        quat : array_like, shape (N, 4) or (4,)
            Each row is a (possibly non-unit norm) quaternion in scalar-last
            (x, y, z, w) format.
        normalized : boolean, optional
            If this flag is `True`, then it is assumed that the input
            quaternions all have unit norm and are not normalized again.
            Default is False.
        """

        return cls(quat, normalized)

    def as_quaternion(self):
        """Return the quaternion representation of the Rotation.

        This function returns a numpy array of shape (4,) or (N x 4) depending
        on the input that was used to initialize the object.
        """
        if self._single:
            return self._quat[0]
        else:
            return self._quat

    def as_dcm(self):
        """Return the direction cosine matrix representation of the Rotation.

        This function returns a numpy.ndarray of shape (3, 3) or (N, 3, 3)
        depending on the input that was used to initialize the object.
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

        num_rotations = self._quat.shape[0]
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

    @classmethod
    def from_dcm(cls, dcm):
        """Initialize rotation from direction cosine matrix.

        This classmethod return a `Rotation` object from the input direction
        cosine matrices. If the input matrix is not orthogonal, this method
        creates an approximation using the algorithm described in [1]_.

        Parameters
        ----------
        dcm : array_like, shape (N, 3, 3) or (3, 3)
            A single matrix or a stack of matrices, where `dcm[i]` is the i-th
            matrix.

        References
        ----------
        .. [1] F. Landis Markley, `Unit Quaternion from Rotation Matrix
               <https://arc.aiaa.org/doi/abs/10.2514/1.31730>`_
        """
        is_single = False
        dcm = np.asarray(dcm, dtype=float)

        if dcm.ndim not in [2, 3] or (dcm.shape[-2], dcm.shape[-1]) != (3, 3):
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
            return cls(quat[0], normalized=True)
        else:
            return cls(quat, normalized=True)

    @classmethod
    def from_rotvec(cls, rotvec):
        """Initialize class from rotation vector.

        A rotation vector is a 3 dimensional vector which is co-directional to
        the axis of rotaion and whose norm gives the angle of rotation (in
        radians).

        Parameters
        ----------
        rotvec : array_like, shape (N, 3) or (3,)
            A single vector or a stack of vectors, where `rot_vec[i]` gives
            the ith rotation vector.
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
        # Use the Taylor expansion of sin(x/2) / x for small angles
        scale[small_angle] = (0.5 - norms[small_angle] ** 2 / 48 +
                              norms[small_angle] ** 4 / 3840)
        scale[large_angle] = (np.sin(norms[large_angle] / 2) /
                              norms[large_angle])

        quat = np.empty((num_rotations, 4))
        quat[:, :3] = scale[:, None] * rotvec
        quat[:, 3] = np.cos(norms / 2)

        if is_single:
            return cls(quat[0], normalized=True)
        else:
            return cls(quat, normalized=True)

    def as_rotvec(self):
        """Return the rotation vector representation of the Rotation.

        This function returns a numpy.ndarray of shape (3,) or (N, 3)
        depending on the input that was used to initialize the object.

        A rotation vector is a 3 dimensional vector which is co-directional to
        the axis of rotation and whose norm gives the angle of rotation (in
        radians).
        """
        quat = self._quat.copy()
        # w > 0 to ensure 0 <= angle <= pi
        quat[quat[:, 3] < 0] *= -1

        angle = 2 * np.arctan2(np.linalg.norm(quat[:, :3], axis=1), quat[:, 3])

        small_angle = (angle <= 1e-3)
        large_angle = ~small_angle

        num_rotations = quat.shape[0]
        scale = np.empty(num_rotations)
        # Use the Taylor expansion of x / sin(x/2) for small angles
        scale[small_angle] = (2 + angle[small_angle] ** 2 / 12 +
                              7 * angle[small_angle] ** 4 / 2880)
        scale[large_angle] = (angle[large_angle] /
                              np.sin(angle[large_angle] / 2))

        rotvec = scale[:, None] * quat[:, :3]

        if self._single:
            return rotvec[0]
        else:
            return rotvec

    @classmethod
    def from_euler(cls, seq, angles, degrees=False):
        """Initialize rotation from Euler angles.

        Parameters
        ----------
        seq : string
            Up to 3 characters belonging to the set {'X', 'Y', 'Z'} for
            intrinsic rotations, or {'x', 'y', 'z'} for extrinsic
            rotations [1]_. Extrinsic and intrinsic rotations cannot be mixed
            in one function call.
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
            If True, then the given angles are taken to be in degrees. Default
            is False.

        References
        ----------
        .. [1] `Euler angle definitions
                <https://en.wikipedia.org/wiki/Euler_angles#Definition_by_intrinsic_rotations>`_
        """
        num_axes = len(seq)
        if num_axes < 1 or num_axes > 3:
            raise ValueError("Expected axis specification to be a non-empty "
                             "string of upto 3 characters, got {}".format(seq))

        intrinsic = (re.compile('^[XYZ]{1,3}$').match(seq) is not None)
        extrinsic = (re.compile('^[xyz]{1,3}$').match(seq) is not None)
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
        return cls(quat[0] if is_single else quat, normalized=True)

    def as_euler(self, seq, degrees=False):
        """Return the Euler angles representation of the Rotation.

        This function returns a numpy.ndarray of shape (N, 3) or (3,) depending
        on how the object was initialized. The algorithm presented in [2]_ has
        been adapted for our use to extract Euler angles. The paper presents an
        algorithm for extracting Euler angles from transformation matrices, as
        opposed to rotation matrices. Thus, the matrix representation used in
        the paper is a transpose of the direction cosine matrix representation
        returned by the `as_dcm` function.

        The returned angles are in the range:

            - First angle belongs to [-180, 180] degrees (both inclusive)
            - Third angle belongs to [-180, 180] degrees (both inclusive)
            - Second angle belongs to:

                - [-90, 90] degrees if all axes are different (like xyz)
                - [0, 180] degrees if first and third axes are the same
                  (like zxz)

        Euler angles suffer from the problem of gimbal lock [3]_. In this case,
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

        References
        ----------
        .. [1] `Euler angle definitions
                <https://en.wikipedia.org/wiki/Euler_angles#Definition_by_intrinsic_rotations>`_
        .. [2] Malcolm D. Shuster, F. Landis Markley
                `General Formula for Euler Angles
                <https://arc.aiaa.org/doi/abs/10.2514/1.16622>`_
        .. [3] `Gimbal lock
                <https://en.wikipedia.org/wiki/Gimbal_lock#In_applied_mathematics>`_
        """
        if len(seq) != 3:
            raise ValueError("Expected 3 axes, got {}.".format(seq))

        intrinsic = (re.compile('^[XYZ]{1,3}$').match(seq) is not None)
        extrinsic = (re.compile('^[xyz]{1,3}$').match(seq) is not None)
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

    def inv(self):
        """Returns the inverse of the current rotation.

        This function returns a new `Rotation` instance containing the inverse
        of all the rotations in the current instance.
        """
        quat = self._quat.copy()
        quat[:, -1] *= -1
        if self._single:
            quat = quat[0]
        return self.__class__(quat, normalized=True)

    def __mul__(self, other):
        """Compose this rotation with the other.

        If `p` and `q` are two rotations, then the composition of 'q followed
        by p' is equivalent to `p * q`. In terms of DCMs, the composition can
        be expressed as `p.as_dcm().dot(q.as_dcm())`.

        This function supports composition of multiple rotations at a time:

            - Either `p` or `q` contains a single rotation. In this case the
              returned object contains the result of composing each rotation in
              the other object with the single rotation.
            - Both `p` and `q` contain `N` rotations. In this case each
              rotation `p[i]` is composed with each rotation `q[i]` and the
              returned object contains `N` rotations.
        """
        if not(self._quat.shape[0] == 1 or other._quat.shape[0] == 1 or
               self._quat.shape[0] == other._quat.shape[0]):
            raise ValueError("Expected equal number of rotations in both "
                             "or a single rotation in either object, "
                             "got {} rotations in first and {} rotations in "
                             "second object.".format(
                                self._quat.shape[0], other._quat.shape[0]))
        result = _compose_quat(self._quat, other._quat)
        if self._single and other._single:
            result = result[0]
        return self.__class__(result, normalized=True)

    def apply(self, vectors, inverse=False):
        """Apply this rotation on a set of vectors.

        Rotates `vectors` by the rotation(s) represented in the object.
        If the original frame rotates to the final frame by this rotation, then
        its application to a vector can be seen in two ways:

            - As a projection of vector components expressed in the final frame
              to the original frame.
            - As the physical rotation of a vector being glued to the original
              frame as it rotates. In this case the vector components are
              expressed in the original frame before and after rotation.

        In terms of DCMs, this application is the same as
        `self.as_dcm().dot(vectors)`.

        The number of rotations and number of vectors given must follow
        standard numpy broadcasting rules: either one of them equals unity or
        they both equal each other.

        Returns a `numpy.ndarray` of shape `(3,)` if object contains a single
        rotation (as opposed to a stack with a single rotation) and a single
        vector is specified with shape `(3,)`. In all other cases, the returned
        array has shape `(N, 3)`.

        Parameters
        ----------
        vectors : array_like, shape (3,) or (P, 3)
            Each `vectors[i]` represents a vector in 3D space. A single vector
            can either be specified with shape `(3, )` or `(1, 3)`.
        inverse : boolean, optional
            If `inverse` is `True` then the inverse of the rotation(s) is
            applied to the input vectors. Default is `False`.
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
        n_rotations = self._quat.shape[0]

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

    def __getitem__(self, indexer):
        """Extract rotation at given index(es) from object.

        This function returns a new `Rotation` instance containing:

            - a single rotation, if `indexer` is a single index
            - a stack of rotation(s), if `indexer` is a slice, or an index
              array. In cases where a single index is ultimately specified, the
              stack will contain a single rotation.

        A single indexer must be specified. The semantics for this function are
        identical to that of numpy arrays and lists.

        Parameters
        ----------
        indexer : index, slice, or index array
            Specifies which rotation(s) to extract.
        """
        # __init__ now copies by default
        return self.__class__(self._quat[indexer], normalized=True)
