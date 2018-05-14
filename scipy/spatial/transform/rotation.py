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
        quat = np.asarray(quat)

        # If a single quaternion is given, convert it to a 2D 1 x 4 matrix but
        # set self._single to True so that we can return appropriate objects
        # in the `to_...` methods
        if quat.shape == (4,):
            quat = quat[None, :]
            self._single = True

        # Each row should have 4 numbers
        if quat.shape[1] != 4:
            raise ValueError("A quaternion should have 4 numbers in (x,y,z,w) format")

        normalized_quat = quat
        if not normalized:
            # L2 norm of each row
            norms = scipy.linalg.norm(quat, axis=1)

            # Warn user for zero (eps?) norm and set to identity quaternion, which
            # is (0,0,0,1) in (x,y,z,w) format
            bool_indices = norms == 0
            if np.any(bool_indices == True):
                warnings.warn("Found zero norm quaternions in input, replacing with identity quaternions")
                quat[bool_indices] = np.array([0, 0, 0, 1])
                norms[bool_indices] = 1.0

            # Normalize each quaternion
            # In case quat is [4 x 4], ensure norm is broadcasted along each column
            normalised_quat = quat / norms[:, None]
        self._quat = normalised_quat

    @classmethod
    def from_quaternion(cls, quat, normalized=False):
        """Initialize Rotation from possibly unnormalized quaternions.

        This classmethod normalizes the input quaternions and returns a
        `Rotation` object.

        Parameters
        ----------
        quat : array_like, shape (N, 4) or (4,)
               Each row is a (possibly non-unit norm) quaternion in scalar-last
               (x, y, z, w) format.

        normalized : boolean, optional, default = False
                     If this flag is `True`, then it is assumed that the input
                     quaternions all have unit norm and are not normalized again.

        """

        return cls(quat, normalized)
