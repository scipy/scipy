from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_array_almost_equal
from scipy.spatial.transform import Rotation


def test_generic_quat_matrix():
    x = np.array([[3,4,0,0], [5, 12, 0, 0]])
    r = Rotation.from_quaternion(x)
    expected_quat = x / np.array([[5], [13]])
    assert_array_almost_equal(r._quat, expected_quat)

def test_from_single_quaternion():
    x = np.array([3,4,0,0])
    r = Rotation.from_quaternion(x)
    expected_quat = x[None, :] / 5
    assert_array_almost_equal(r._quat, expected_quat)

def test_from_square_quat_matrix():
    # Ensure proper norm array broadcasting
    x = np.array([
        [3, 0, 0, 4],
        [5, 0, 12, 0],
        [0, 0, 0, 1],
        [0, 0, 0, -1]
        ])
    r = Rotation.from_quaternion(x)
    expected_quat = x / np.array([[5], [13], [1], [1]])
    assert_array_almost_equal(r._quat, expected_quat)
