from __future__ import division, print_function, absolute_import

import pytest

import numpy as np
from numpy.testing import assert_equal, assert_array_almost_equal
from scipy.spatial.transform import Rotation


def test_generic_quat_matrix():
    x = np.array([[3, 4, 0, 0], [5, 12, 0, 0]])
    r = Rotation.from_quaternion(x)
    expected_quat = x / np.array([[5], [13]])
    assert_array_almost_equal(r.as_quaternion(), expected_quat)


def test_from_single_1d_quaternion():
    x = np.array([3, 4, 0, 0])
    r = Rotation.from_quaternion(x)
    expected_quat = x / 5
    assert_array_almost_equal(r.as_quaternion(), expected_quat)


def test_from_single_2d_quaternion():
    x = np.array([[3, 4, 0, 0]])
    r = Rotation.from_quaternion(x)
    expected_quat = x / 5
    assert_array_almost_equal(r.as_quaternion(), expected_quat)


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
    assert_array_almost_equal(r.as_quaternion(), expected_quat)


def test_malformed_1d_from_quaternion():
    with pytest.raises(ValueError):
        Rotation.from_quaternion(np.array([1, 2, 3]))


def test_malformed_2d_from_quaternion():
    with pytest.raises(ValueError):
        Rotation.from_quaternion(np.array([
            [1, 2, 3, 4, 5],
            [4, 5, 6, 7, 8]
            ]))


def test_zero_norms_from_quaternion():
    x = np.array([
            [3, 4, 0, 0],
            [0, 0, 0, 0],
            [5, 0, 12, 0]
            ])
    with pytest.raises(ValueError):
        r = Rotation.from_quaternion(x)


def test_as_dcm_single_1d_quaternion():
    quat = [0, 0, 0, 1]
    mat = Rotation.from_quaternion(quat).as_dcm()
    # mat.shape == (3,3) due to 1d input
    assert_array_almost_equal(mat, np.eye(3))


def test_as_dcm_single_2d_quaternion():
    quat = [[0, 0, 1, 1]]
    mat = Rotation.from_quaternion(quat).as_dcm()
    assert_equal(mat.shape, (1, 3, 3))
    expected_mat = np.array([
        [0, -1, 0],
        [1, 0, 0],
        [0, 0, 1]
        ])
    assert_array_almost_equal(mat[0], expected_mat)


def test_as_dcm_from_square_input():
    quats = [
            [0, 0, 1, 1],
            [0, 1, 0, 1],
            [0, 0, 0, 1],
            [0, 0, 0, -1]
            ]
    mat = Rotation.from_quaternion(quats).as_dcm()
    assert_equal(mat.shape, (4, 3, 3))

    expected0 = np.array([
        [0, -1, 0],
        [1, 0, 0],
        [0, 0, 1]
        ])
    assert_array_almost_equal(mat[0], expected0)

    expected1 = np.array([
        [0, 0, 1],
        [0, 1, 0],
        [-1, 0, 0]
        ])
    assert_array_almost_equal(mat[1], expected1)

    assert_array_almost_equal(mat[2], np.eye(3))
    assert_array_almost_equal(mat[3], np.eye(3))


def test_as_dcm_from_generic_input():
    quats = [
            [0, 0, 1, 1],
            [0, 1, 0, 1],
            [1, 2, 3, 4]
            ]
    mat = Rotation.from_quaternion(quats).as_dcm()
    assert_equal(mat.shape, (3, 3, 3))

    expected0 = np.array([
        [0, -1, 0],
        [1, 0, 0],
        [0, 0, 1]
        ])
    assert_array_almost_equal(mat[0], expected0)

    expected1 = np.array([
        [0, 0, 1],
        [0, 1, 0],
        [-1, 0, 0]
        ])
    assert_array_almost_equal(mat[1], expected1)

    expected2 = np.array([
        [0.4, -2, 2.2],
        [2.8, 1, 0.4],
        [-1, 2, 2]
        ]) / 3
    assert_array_almost_equal(mat[2], expected2)
