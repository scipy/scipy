from __future__ import division, print_function, absolute_import

import pytest

import numpy as np
from numpy.testing import assert_equal, assert_array_almost_equal
from numpy.testing import assert_allclose
from scipy.spatial.transform import Rotation
from scipy.stats import special_ortho_group


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
        Rotation.from_quaternion(x)


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


def test_from_single_2d_dcm():
    dcm = [
            [0, 0, 1],
            [1, 0, 0],
            [0, 1, 0]
            ]
    expected_quat = [0.5, 0.5, 0.5, 0.5]
    assert_array_almost_equal(
            Rotation.from_dcm(dcm).as_quaternion(),
            expected_quat)


def test_from_single_3d_dcm():
    dcm = np.array([
        [0, 0, 1],
        [1, 0, 0],
        [0, 1, 0]
        ]).reshape((1, 3, 3))
    expected_quat = np.array([0.5, 0.5, 0.5, 0.5]).reshape((1, 4))
    assert_array_almost_equal(
            Rotation.from_dcm(dcm).as_quaternion(),
            expected_quat)


def test_from_dcm_calculation():
    expected_quat = np.array([1, 1, 6, 1]) / np.sqrt(39)
    dcm = np.array([
            [-0.8974359, -0.2564103, 0.3589744],
            [0.3589744, -0.8974359, 0.2564103],
            [0.2564103, 0.3589744, 0.8974359]
            ])
    assert_array_almost_equal(
            Rotation.from_dcm(dcm).as_quaternion(),
            expected_quat)
    assert_array_almost_equal(
            Rotation.from_dcm(dcm.reshape((1, 3, 3))).as_quaternion(),
            expected_quat.reshape((1, 4)))


def test_dcm_calculation_pipeline():
    dcm = special_ortho_group.rvs(3, size=10, random_state=0)
    assert_array_almost_equal(Rotation.from_dcm(dcm).as_dcm(), dcm)


def test_from_dcm_ortho_output():
    np.random.seed(0)
    dcm = np.random.random((100, 3, 3))
    ortho_dcm = Rotation.from_dcm(dcm).as_dcm()

    mult_result = np.einsum('...ij,...jk->...ik', ortho_dcm,
                            ortho_dcm.transpose((0, 2, 1)))

    eye3d = np.zeros((100, 3, 3))
    for i in range(3):
        eye3d[:, i, i] = 1.0

    assert_array_almost_equal(mult_result, eye3d)


def test_from_1d_single_rotvec():
    rotvec = [1, 0, 0]
    expected_quat = np.array([0.4794255, 0, 0, 0.8775826])
    result = Rotation.from_rotvec(rotvec)
    assert_array_almost_equal(result.as_quaternion(), expected_quat)


def test_from_2d_single_rotvec():
    rotvec = [[1, 0, 0]]
    expected_quat = np.array([[0.4794255, 0, 0, 0.8775826]])
    result = Rotation.from_rotvec(rotvec)
    assert_array_almost_equal(result.as_quaternion(), expected_quat)


def test_from_generic_rotvec():
    rotvec = [
            [1, 2, 2],
            [1, -1, 0.5],
            [0, 0, 0]
            ]
    expected_quat = np.array([
        [0.3324983, 0.6649967, 0.6649967, 0.0707372],
        [0.4544258, -0.4544258, 0.2272129, 0.7316889],
        [0, 0, 0, 1]
        ])
    assert_array_almost_equal(
            Rotation.from_rotvec(rotvec).as_quaternion(),
            expected_quat)


def test_from_rotvec_small_angle():
    rotvec = np.array([
        [5e-4 / np.sqrt(3), -5e-4 / np.sqrt(3), 5e-4 / np.sqrt(3)],
        [0.2, 0.3, 0.4],
        [0, 0, 0]
        ])

    quat = Rotation.from_rotvec(rotvec).as_quaternion()
    # cos(theta/2) ~~ 1 for small theta
    assert_allclose(quat[0, 3], 1)
    # sin(theta/2) / theta ~~ 0.5 for small theta
    assert_allclose(quat[0, :3], rotvec[0] * 0.5)

    assert_allclose(quat[1, 3], 0.9639685)
    assert_allclose(
            quat[1, :3],
            np.array([
                0.09879603932153465,
                0.14819405898230198,
                0.19759207864306931
                ]))

    assert_equal(quat[2], np.array([0, 0, 0, 1]))


def test_malformed_1d_from_rotvec():
    with pytest.raises(ValueError, match='Expected `rot_vec` to have shape'):
        Rotation.from_rotvec([1, 2])


def test_malformed_2d_from_rotvec():
    with pytest.raises(ValueError, match='Expected `rot_vec` to have shape'):
        Rotation.from_rotvec([
            [1, 2, 3, 4],
            [5, 6, 7, 8]
            ])


def test_as_generic_rotvec():
    quat = np.array([
            [1, 2, -1, 0.5],
            [1, -1, 1, 0.0003],
            [0, 0, 0, 1]
            ])
    quat /= np.linalg.norm(quat, axis=1)[:, None]

    rotvec = Rotation.from_quaternion(quat).as_rotvec()
    angle = np.linalg.norm(rotvec, axis=1)

    assert_allclose(quat[:, 3], np.cos(angle/2))
    assert_allclose(np.cross(rotvec, quat[:, :3]), np.zeros((3, 3)))


def test_as_rotvec_single_1d_input():
    quat = np.array([1, 2, -3, 2])
    expected_rotvec = np.array([0.5772381, 1.1544763, -1.7317144])

    actual_rotvec = Rotation.from_quaternion(quat).as_rotvec()

    assert_equal(actual_rotvec.shape, (3,))
    assert_allclose(actual_rotvec, expected_rotvec)


def test_as_rotvec_single_2d_input():
    quat = np.array([[1, 2, -3, 2]])
    expected_rotvec = np.array([[0.5772381, 1.1544763, -1.7317144]])

    actual_rotvec = Rotation.from_quaternion(quat).as_rotvec()

    assert_equal(actual_rotvec.shape, (1, 3))
    assert_allclose(actual_rotvec, expected_rotvec)


def test_rotvec_calc_pipeline():
    # Include small angles
    rotvec = np.array([
        [0, 0, 0],
        [1, -1, 2],
        [-3e-4, 3.5e-4, 7.5e-5]
        ])
    assert_allclose(Rotation.from_rotvec(rotvec).as_rotvec(), rotvec)


def test_from_euler_single_rotation():
    quat = Rotation.from_euler('z', 90, degrees=True).as_quaternion()
    expected_quat = np.array([0, 0, 1, 1]) / np.sqrt(2)
    assert_allclose(quat, expected_quat)


def test_single_intrinsic_extrinsic_rotation():
    ext = Rotation.from_euler('z', 90, degrees=True).as_dcm()
    int = Rotation.from_euler('Z', 90, degrees=True).as_dcm()
    assert_allclose(ext, int)


def test_from_euler_rotation_order():
    # Intrinsic rotation is same as extrinsic with order reversed
    np.random.seed(0)
    a = np.random.randint(low=0, high=180, size=(6, 3))
    b = a[:, np.argsort([2, 1, 0])]
    x = Rotation.from_euler('xyz', a, degrees=True).as_quaternion()
    y = Rotation.from_euler('ZYX', b, degrees=True).as_quaternion()
    assert_allclose(x, y)


def test_from_euler_elementary_extrinsic_rotation():
    # Simple test to check if extrinsic rotations are implemented correctly
    dcm = Rotation.from_euler('zx', [90, 90], degrees=True).as_dcm()
    expected_dcm = np.array([
        [0, -1, 0],
        [0, 0, -1],
        [1, 0, 0]
    ])
    assert_array_almost_equal(dcm, expected_dcm)
