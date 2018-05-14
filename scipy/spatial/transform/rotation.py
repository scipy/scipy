from __future__ import division, print_function, absolute_import

__all__ = ['Rotation']

import numpy as np
import scipy.linalg

class Rotation(object):
    """Rotation in 3 dimensions.

    This class will include initializers from different representations, converters
    and some useful algorithms such as quaternion slerp and rotation estimation.

    """
    def __init__(self, quaternions):
        """Initialize class from normalized quaternions.

        Parameters
        ----------
        quaternions : (N, 4) numpy.array
                      Each row is a unit-norm quaternion stored in scalar-last
                      (x, y, z, w) format.

        """
        self._quat = quaternions

    @classmethod
    def from_quaternion(cls, quat):
        """Initialize Rotation from possibly unnormalized quaternions

        This classmethod normalizes the input quaternions and returns a
        `Rotation` object.

        Parameters
        ----------
        quat : (N, 4) or (4,) numpy.array
               Each row is a (possibly non-unit norm) quaternion in scalar-last
               (x, y, z, w) format.

        """

        # If a single quaternion is given, convert it to a 2D 1 x 4 matrix
        if quat.shape == (4,):
            quat = quat[None, :]

        # Each row should have 4 numbers
        assert(quat.shape[1] == 4)

        # L2 norm of each row
        norms = scipy.linalg.norm(quat, axis=1)

        # Divide each row by its norm.
        # In case quat is [4 x 4], ensure norm is broadcasted along each column
        normalised_quat = quat / norms[:, None]

        return cls(normalised_quat)
