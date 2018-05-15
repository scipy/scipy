from __future__ import division, print_function, absolute_import

import numpy as np
import scipy.linalg
import warnings

class Rotation(object):
    """Rotation in 3 dimensions.

    This class will include initializers from different representations, converters
    and some useful algorithms such as quaternion slerp and rotation estimation.

    For initializing Rotations usage of `from_...` methods such as `from_quaternion`
    is recommended instead of using `__init__`.

    Methods
    -------
    from_quaternion
    """
    def __init__(self, quat, normalized=False):
        self._single = False
        # Try to convert to numpy array
        quat = np.asarray(quat, dtype=float)

        if quat.ndim not in [1, 2]:
            raise ValueError("Input quat should be of shape (4,) or (N x 4).")

        # If a single quaternion is given, convert it to a 2D 1 x 4 matrix but
        # set self._single to True so that we can return appropriate objects
        # in the `to_...` methods
        if quat.ndim == 1:
            quat = quat[None, :]
            self._single = True

        # Each row should have 4 numbers
        if quat.shape[1] != 4:
            raise ValueError("A quaternion should have 4 numbers.")

        self._quat = quat

        if not normalized:
            # L2 norm of each row
            norms = scipy.linalg.norm(quat, axis=1)

            # Warn user for zero (eps?) norm and set to identity quaternion,
            # which (0,0,0,1) in (x,y,z,w) format
            zero_norms = norms == 0
            if zero_norms.any():
                warnings.warn("Found zero norm quaternions in input, replacing with identity quaternions.")
                self._quat[zero_norms] = np.array([0, 0, 0, 1])
            # Normalize each quaternion, ensuring norm is broadcasted along
            # each column.
            self._quat[~zero_norms] /= norms[~zero_norms][:, None]

    @classmethod
    def from_quaternion(cls, quat, normalized=False):
        """Initialize Rotation from quaternions.

        This classmethod returns a `Rotation` object from the input quaternions.
        If `normalized` is `True`, then the quaternions are assumed to have
        unit norm, else the quaternions are normalized before the object is
        created.

        Parameters
        ----------
        quat : array_like, shape (N, 4) or (4,)
            Each row is a (possibly non-unit norm) quaternion in scalar-last
            (x, y, z, w) format.

        normalized : boolean, optional
            If this flag is `True`, then it is assumed that the input quaternions
            all have unit norm and are not normalized again. Default is False.
        """

        return cls(quat, normalized)
