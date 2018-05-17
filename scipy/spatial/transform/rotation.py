from __future__ import division, print_function, absolute_import

import numpy as np
import scipy.linalg
import warnings


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
    as_dcm
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
    def from_dcm(cls, mat, scalar=True):
        # For now, mat must be a 3D matrix
        num_rotations = mat.shape[0]
        if scalar:
            quat = np.empty((num_rotations, 4))
            for rot_num in range(num_rotations):
                A = mat[rot_num]
                trA = A.trace()
                choice = np.argmax(np.array([A[0, 0], A[1, 1], A[2, 2], trA]))

                if choice in [0, 1, 2]:
                    quat[rot_num, choice] = 1 - trA + 2 * A[choice, choice]
                    next_ind = (choice + 1) % 3
                    next_next = (next_ind + 1) % 3

                    quat[rot_num, next_ind] = A[choice, next_ind] + A[next_ind, choice]
                    quat[rot_num, next_next] = A[choice, next_next] + A[next_next, choice]

                    quat[rot_num, 3] = A[next_ind, next_next] - A[next_next, next_ind]
                else:
                    quat[rot_num] = np.array([
                        A[1, 2] - A[2, 1],
                        A[2, 0] - A[0, 2],
                        A[0, 1] - A[1, 0],
                        1 + trA])
        else:
            decision_matrix = np.zeros((num_rotations, 4))
            decision_matrix[:, :3] = mat.diagonal(axis1=1, axis2=2)
            decision_matrix[:, -1] = decision_matrix.sum(axis=1)
            choices = decision_matrix.argmax(axis=1)

            all_quat = np.empty((num_rotations, 4, 4))

            all_quat[:, 0, 0] = 1 - decision_matrix[:, -1] + 2 * mat[:, 0, 0]
            all_quat[:, 0, 1] = mat[:, 0, 1] + mat[:, 1, 0]
            all_quat[:, 0, 2] = mat[:, 0, 2] + mat[:, 2, 0]
            all_quat[:, 0, 3] = mat[:, 1, 2] - mat[:, 2, 1]

            all_quat[:, 1, 0] = mat[:, 1, 0] + mat[:, 0, 1]
            all_quat[:, 1, 1] = 1 - decision_matrix[:, -1] + 2 * mat[:, 1, 1]
            all_quat[:, 1, 2] = mat[:, 1, 2] + mat[:, 2, 1]
            all_quat[:, 1, 3] = mat[:, 2, 0] - mat[:, 0, 2]

            all_quat[:, 2, 0] = mat[:, 2, 0] + mat[:, 0, 2]
            all_quat[:, 2, 1] = mat[:, 2, 1] + mat[:, 1, 2]
            all_quat[:, 2, 2] = 1 - decision_matrix[:, -1] + 2 * mat[:, 2, 2]
            all_quat[:, 2, 3] = mat[:, 0, 1] - mat[:, 1, 0]

            all_quat[:, 3, 0] = mat[:, 1, 2] - mat[:, 2, 1]
            all_quat[:, 3, 1] = mat[:, 2, 0] - mat[:, 0, 2]
            all_quat[:, 3, 2] = mat[:, 0, 1] - mat[:, 1, 0]
            all_quat[:, 3, 3] = 1 + decision_matrix[:, -1]

            quat = np.array([np.take(all_quat[i], choices[i], axis=0) for i in range(num_rotations)])

        # For testing purposes only. Return cls(quat) in final implementation.
        # Creating new objects in python is time consuming. We do not want to
        # count that time in the benchmarks.
        return quat
