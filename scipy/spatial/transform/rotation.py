from __future__ import division, print_function, absolute_import

import numpy as np
import scipy.linalg
import re


AXIS_TO_IND = {'x': 0, 'y': 1, 'z': 2}


def _elementary_basis_vector(char):
    b = np.zeros(3)
    b[AXIS_TO_IND[char]] = 1
    return b


def _make_euler_from_dcm(dcm, seq, extrinsic=False):
    # The algorithm assumes intrinsic frame transformations. Their
    # ELEMENTARY dcms are transpose of what we return.
    # Adapt the algorithm for our case by
    # 1. Transposing our dcm representation
    # 2. Reversing both axis sequence and angles for extrinsic rotations

    seq = seq[::-1] if extrinsic else seq

    if dcm.ndim == 2:
        dcm = dcm[None, :, :]
    num_rotations = dcm.shape[0]

    # Step 0
    # Algorithms dcms are transpose of ours
    dcm = np.transpose(dcm, (0, 2, 1))
    # Algorithm assumes axes as column vectors, here we use 1D vectors
    n1 = _elementary_basis_vector(seq[0])
    n2 = _elementary_basis_vector(seq[1])
    n3 = _elementary_basis_vector(seq[2])

    # Step 2
    lamb = np.arctan2(np.dot(np.cross(n1, n2), n3), np.dot(n1, n3))
    c = np.empty((3, 3))
    c[0] = n2
    c[1] = np.cross(n1, n2)
    c[2] = n1

    # Step 3
    cl = np.cos(lamb)
    sl = np.sin(lamb)
    rt = np.array([
        [1, 0, 0],
        [0, cl, -sl],
        [0, sl, cl]
    ])
    rtc = rt.dot(c)
    ct = c.T
    o = np.empty_like(dcm)
    for ind in range(num_rotations):
        o[ind] = rtc.dot(dcm[ind]).dot(ct)

    # Step 4
    angle2 = lamb + np.arccos(o[:, 2, 2])

    # Steps 5, 6
    eps = np.finfo(float).resolution  # ~1e-15
    safe1 = (np.abs(angle2 - lamb) >= eps)
    safe2 = (np.abs(angle2 - lamb - np.pi) >= eps)

    angle1 = np.empty(num_rotations)
    angle3 = np.empty(num_rotations)

    # 5b
    safe_mask = np.logical_and(safe1, safe2)
    angle1[safe_mask] = np.arctan2(o[safe_mask, 2, 0], -o[safe_mask, 2, 1])
    angle3[safe_mask] = np.arctan2(o[safe_mask, 0, 2], o[safe_mask, 1, 2])

    # 6a
    angle3[~safe_mask] = 0
    # 6b
    angle1[~safe1] = np.arctan2(o[~safe1, 0, 1] - o[~safe1, 1, 0],
                                o[~safe1, 0, 0] + o[~safe1, 1, 1])
    # 6c
    angle1[~safe2] = np.arctan2(o[~safe2, 0, 1] + o[~safe2, 1, 0],
                                o[~safe2, 0, 0] - o[~safe2, 1, 1])

    # Step 7
    # python modulo operator works correctly for negative numbers
    adjust_mask = np.logical_or(angle2 < 0, angle2 > np.pi)
    angle1[adjust_mask] = (angle1[adjust_mask] + np.pi) % (2 * np.pi)
    angle2[adjust_mask] = (2 * lamb - angle2[adjust_mask]) % (2 * np.pi)
    angle3[adjust_mask] = (angle3[adjust_mask] - np.pi) % (2 * np.pi)

    # Step 8
    # TODO: if any observability flags are poor, possibly raise a UserWarning?

    # Reverse role of extrinsic and intrinsic rotations
    return np.column_stack((angle3, angle2, angle1) if extrinsic else
                           (angle1, angle2, angle3))


def _make_elementary_quat(axis, angles):
    num_rotations = angles.shape[0]
    quat = np.zeros((num_rotations, 4))

    quat[:, 3] = np.cos(angles / 2)
    quat[:, AXIS_TO_IND[axis]] = np.sin(angles / 2)
    return quat


def _compose_quat(p, q):
    # p and q should have same shape (N, 4)
    product = np.empty_like(p)
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
    """
    def __init__(self, quat, normalized=False):
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
            self._quat = quat
        else:
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
        on how the object was initialized.

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
        """
        if len(seq) != 3:
            raise ValueError("Expected 3 axes, got {}.".format(seq))

        intrinsic = (re.compile('^[XYZ]{1,3}$').match(seq) is not None)
        extrinsic = (re.compile('^[xyz]{1,3}$').match(seq) is not None)
        if not (intrinsic or extrinsic):
            raise ValueError("Expected axes from `seq` to be from ['x', 'y', "
                             "'z'] or ['X', 'Y', 'Z'], got {}".format(seq))

        if any(seq[i] == seq[i+1] for i in range(2)):
            raise ValueError("Expected consecutive axes to be different, "
                             "got {}".format(seq))

        seq = seq.lower()

        angles = _make_euler_from_dcm(self.as_dcm(), seq, extrinsic)
        if degrees:
            angles = np.rad2deg(angles)

        return angles[0] if self._single else angles
