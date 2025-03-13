import pytest

import numpy as np
from numpy.testing import assert_equal, assert_array_almost_equal
from numpy.testing import assert_allclose
from scipy.spatial.transform import Rotation, Slerp
from scipy.stats import special_ortho_group
from itertools import permutations
from scipy._lib._array_api import (
    xp_assert_equal,
    is_numpy,
    is_array_api_strict,
    is_torch,
)
import scipy._lib.array_api_extra as xpx

import pickle
import copy


def basis_vec(axis):
    if axis == "x":
        return [1, 0, 0]
    elif axis == "y":
        return [0, 1, 0]
    elif axis == "z":
        return [0, 0, 1]


def test_init_non_array():
    Rotation((0, 0, 0, 1))
    Rotation([0, 0, 0, 1])


def test_init(xp):
    x = xp.asarray([0, 0, 0, 1])
    Rotation(x)


def test_generic_quat_matrix(xp):
    x = xp.asarray([[3, 4, 0, 0], [5, 12, 0, 0]], dtype=xp.float64)
    r = Rotation.from_quat(x)
    expected_quat = x / xp.asarray([[5], [13.0]])
    assert_array_almost_equal(r.as_quat(), expected_quat)


def test_from_single_1d_quaternion(xp):
    x = xp.asarray([3.0, 4, 0, 0])
    r = Rotation.from_quat(x)
    expected_quat = x / 5
    assert_array_almost_equal(r.as_quat(), expected_quat)


def test_from_single_2d_quaternion(xp):
    x = xp.asarray([[3, 4, 0, 0]])
    r = Rotation.from_quat(x)
    expected_quat = xp.astype(x, r._quat.dtype) / 5.0
    assert_array_almost_equal(r.as_quat(), expected_quat)


def test_from_quat_scalar_first(xp):
    rng = np.random.RandomState(0)

    r = Rotation.from_quat([1, 0, 0, 0], scalar_first=True)
    assert_allclose(r.as_matrix(), xp.eye(3), rtol=1e-15, atol=1e-16)

    r = Rotation.from_quat(
        xp.tile(xp.asarray([1, 0, 0, 0]), (10, 1)), scalar_first=True
    )
    assert_allclose(
        r.as_matrix(), xp.tile(xp.eye(3), (10, 1, 1)), rtol=1e-15, atol=1e-16
    )

    q = xp.asarray(rng.randn(100, 4))
    q /= xp.linalg.vector_norm(q, axis=1)[:, None]
    for i in range(q.shape[0]):  # Array API conforming loop
        qi = q[i, ...]
        r = Rotation.from_quat(qi, scalar_first=True)
        assert_allclose(xp.roll(r.as_quat(), 1), qi, rtol=1e-15)

    r = Rotation.from_quat(q, scalar_first=True)
    assert_allclose(xp.roll(r.as_quat(), 1, axis=1), q, rtol=1e-15)


def test_from_quat_jax_compile():
    pytest.importorskip("jax")
    import jax

    q = jax.numpy.array([0.0, 0.0, 0.0, 1.0])
    r = jax.block_until_ready(jax.jit(Rotation.from_quat)(q))
    assert isinstance(r, Rotation)


def test_as_quat_scalar_first(xp):
    rng = np.random.RandomState(0)

    r = Rotation.from_euler("xyz", xp.zeros(3))
    assert_allclose(r.as_quat(scalar_first=True), [1, 0, 0, 0], rtol=1e-15, atol=1e-16)

    r = Rotation.from_euler("xyz", np.zeros((10, 3)))
    assert_allclose(
        r.as_quat(scalar_first=True),
        xp.tile(xp.asarray([1, 0, 0, 0]), (10, 1)),
        rtol=1e-15,
        atol=1e-16,
    )

    q = xp.asarray(rng.randn(100, 4))
    q /= xp.linalg.vector_norm(q, axis=1)[:, None]
    for i in range(q.shape[0]):  # Array API conforming loop
        qi = q[i, ...]
        r = Rotation.from_quat(qi)
        assert_allclose(r.as_quat(scalar_first=True), xp.roll(qi, 1), rtol=1e-15)

        assert_allclose(
            r.as_quat(canonical=True, scalar_first=True),
            xp.roll(r.as_quat(canonical=True), 1),
            rtol=1e-15,
        )

    r = Rotation.from_quat(q)
    assert_allclose(r.as_quat(scalar_first=True), xp.roll(q, 1, axis=1), rtol=1e-15)

    assert_allclose(
        r.as_quat(canonical=True, scalar_first=True),
        xp.roll(r.as_quat(canonical=True), 1, axis=1),
        rtol=1e-15,
    )


def test_from_square_quat_matrix(xp):
    # Ensure proper norm array broadcasting
    x = xp.asarray(
        [
            [3.0, 0, 0, 4],
            [5, 0, 12, 0],
            [0, 0, 0, 1],
            [-1, -1, -1, 1],
            [0, 0, 0, -1],  # Check double cover
            [-1, -1, -1, -1],  # Check double cover
        ]
    )
    r = Rotation.from_quat(x)
    expected_quat = x / xp.asarray([[5.0], [13], [1], [2], [1], [2]])
    assert_array_almost_equal(r.as_quat(), expected_quat)


def test_quat_double_to_canonical_single_cover(xp):
    x = xp.asarray(
        [[-1.0, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1], [-1, -1, -1, -1]]
    )
    r = Rotation.from_quat(x)
    expected_quat = xp.abs(x) / xp.linalg.vector_norm(x, axis=1)[:, None]
    assert_allclose(r.as_quat(canonical=True), expected_quat)


def test_quat_double_cover(xp):
    # See the Rotation.from_quat() docstring for scope of the quaternion
    # double cover property.
    # Check from_quat and as_quat(canonical=False)
    q = xp.asarray([0.0, 0, 0, -1])
    r = Rotation.from_quat(q)
    xp_assert_equal(q, r.as_quat(canonical=False))

    # Check composition and inverse
    q = xp.asarray([1.0, 0, 0, 1]) / np.sqrt(2)  # 90 deg rotation about x
    r = Rotation.from_quat(q)
    r3 = r * r * r
    assert_allclose(r.as_quat(canonical=False) * np.sqrt(2), [1, 0, 0, 1])
    assert_allclose(r.inv().as_quat(canonical=False) * np.sqrt(2), [-1, 0, 0, 1])
    assert_allclose(r3.as_quat(canonical=False) * np.sqrt(2), [1, 0, 0, -1])
    assert_allclose(r3.inv().as_quat(canonical=False) * np.sqrt(2), [-1, 0, 0, -1])

    # More sanity checks
    assert_allclose((r * r.inv()).as_quat(canonical=False), [0, 0, 0, 1], atol=2e-16)
    assert_allclose((r3 * r3.inv()).as_quat(canonical=False), [0, 0, 0, 1], atol=2e-16)
    assert_allclose((r * r3).as_quat(canonical=False), [0, 0, 0, -1], atol=2e-16)
    assert_allclose(
        (r.inv() * r3.inv()).as_quat(canonical=False), [0, 0, 0, -1], atol=2e-16
    )


def test_as_quat_jax_compile():
    pytest.importorskip("jax")
    import jax

    r = Rotation.from_quat(jax.numpy.array([1, 0, 0, 0]))
    jax.block_until_ready(jax.jit(Rotation.as_quat)(r))


def test_from_quat_wrong_shape(xp):
    # Wrong shape 1d array
    with pytest.raises(ValueError, match="Expected `quat` to have shape"):
        Rotation.from_quat(xp.asarray([1, 2, 3]))

    # Wrong shape 2d array
    with pytest.raises(ValueError, match="Expected `quat` to have shape"):
        Rotation.from_quat(xp.asarray([[1, 2, 3, 4, 5], [4, 5, 6, 7, 8]]))

    # 3d array
    with pytest.raises(ValueError, match="Expected `quat` to have shape"):
        Rotation.from_quat(xp.asarray([[[1, 2, 3, 4]], [[4, 5, 6, 7]]]))

    # 0-length 2d array
    with pytest.raises(ValueError, match="Expected `quat` to have shape"):
        Rotation.from_quat(xp.empty((0, 4)))


def test_zero_norms_from_quat(xp):
    x = xp.asarray([[3, 4, 0, 0], [0, 0, 0, 0], [5, 0, 12, 0]])
    if is_numpy(xp):
        with pytest.raises(ValueError):
            Rotation.from_quat(x)
    else:
        assert xp.all(xp.isnan(Rotation.from_quat(x).as_quat()[1, ...]))


def test_as_matrix_single_1d_quaternion(xp):
    quat = xp.asarray([0, 0, 0, 1])
    mat = Rotation.from_quat(quat).as_matrix()
    # mat.shape == (3,3) due to 1d input
    assert_array_almost_equal(mat, xp.eye(3))


def test_as_matrix_single_2d_quaternion(xp):
    quat = xp.asarray([[0, 0, 1, 1]])
    mat = Rotation.from_quat(quat).as_matrix()
    assert_equal(mat.shape, (1, 3, 3))
    expected_mat = xp.asarray([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    assert_array_almost_equal(mat[0, ...], expected_mat)


def test_as_matrix_from_square_input(xp):
    quats = xp.asarray([[0, 0, 1, 1], [0, 1, 0, 1], [0, 0, 0, 1], [0, 0, 0, -1]])
    mat = Rotation.from_quat(quats).as_matrix()
    assert_equal(mat.shape, (4, 3, 3))

    expected0 = xp.asarray([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    assert_array_almost_equal(mat[0, ...], expected0)

    expected1 = xp.asarray([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
    assert_array_almost_equal(mat[1, ...], expected1)

    assert_array_almost_equal(mat[2, ...], xp.eye(3))
    assert_array_almost_equal(mat[3, ...], xp.eye(3))


def test_as_matrix_from_generic_input(xp):
    quats = xp.asarray([[0, 0, 1, 1], [0, 1, 0, 1], [1, 2, 3, 4]])
    mat = Rotation.from_quat(quats).as_matrix()
    assert_equal(mat.shape, (3, 3, 3))

    expected0 = xp.asarray([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    assert_array_almost_equal(mat[0, ...], expected0)

    expected1 = xp.asarray([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
    assert_array_almost_equal(mat[1, ...], expected1)

    expected2 = xp.asarray([[0.4, -2, 2.2], [2.8, 1, 0.4], [-1, 2, 2]]) / 3
    assert_array_almost_equal(mat[2, ...], expected2)


def test_as_matrix_jax_compile():
    pytest.importorskip("jax")
    import jax

    r = Rotation.from_matrix(jax.numpy.eye(3))
    jax.block_until_ready(jax.jit(Rotation.as_matrix)(r))


def test_from_single_2d_matrix(xp):
    mat = xp.asarray([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    expected_quat = xp.asarray([0.5, 0.5, 0.5, 0.5])
    assert_array_almost_equal(Rotation.from_matrix(mat).as_quat(), expected_quat)


def test_from_single_3d_matrix(xp):
    mat = xp.reshape(xp.asarray([[0, 0, 1], [1, 0, 0], [0, 1, 0]]), (1, 3, 3))
    expected_quat = xp.reshape(xp.asarray([0.5, 0.5, 0.5, 0.5]), (1, 4))
    assert_array_almost_equal(Rotation.from_matrix(mat).as_quat(), expected_quat)


def test_from_matrix_calculation(xp):
    expected_quat = xp.asarray([1.0, 1, 6, 1]) / np.sqrt(39)
    mat = xp.asarray(
        [
            [-0.8974359, -0.2564103, 0.3589744],
            [0.3589744, -0.8974359, 0.2564103],
            [0.2564103, 0.3589744, 0.8974359],
        ]
    )
    assert_array_almost_equal(Rotation.from_matrix(mat).as_quat(), expected_quat)
    assert_array_almost_equal(
        Rotation.from_matrix(xp.reshape(mat, (1, 3, 3))).as_quat(),
        xp.reshape(expected_quat, (1, 4)),
    )


def test_matrix_calculation_pipeline(xp):
    mat = xp.asarray(special_ortho_group.rvs(3, size=10, random_state=0))
    assert_array_almost_equal(Rotation.from_matrix(mat).as_matrix(), mat)


def test_from_matrix_ortho_output(xp):
    rnd = np.random.RandomState(0)
    mat = xp.asarray(rnd.random_sample((100, 3, 3)))
    dets = xp.linalg.det(mat)
    for i in range(dets.shape[0]):
        # Make sure we have a right-handed rotation matrix
        if dets[i] < 0:
            mat = xpx.at(mat)[i, ...].set(-mat[i, ...])
    ortho_mat = Rotation.from_matrix(mat).as_matrix()

    mult_result = xp.matmul(ortho_mat, xp.matrix_transpose(ortho_mat))

    eye3d = xp.zeros((100, 3, 3)) + xp.eye(3)

    assert_array_almost_equal(mult_result, eye3d)


def test_from_matrix_normalize(xp):
    mat = xp.asarray([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
    expected = xp.asarray(
        [[0.894427, 0.447214, 0.0], [-0.447214, 0.894427, 0.0], [0.0, 0.0, 1.0]]
    )
    assert_allclose(Rotation.from_matrix(mat).as_matrix(), expected, atol=1e-6)

    mat = xp.asarray([[0, -0.5, 0], [0.5, 0, 0], [0, 0, 0.5]])
    expected = xp.asarray([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    assert_allclose(Rotation.from_matrix(mat).as_matrix(), expected, atol=1e-6)


def test_from_matrix_non_positive_determinant(xp):
    mat = xp.eye(3)
    mat = xpx.at(mat)[0, 0].set(0)
    # TODO: Unify the error response for all backends
    if is_numpy(xp):
        with pytest.raises(ValueError, match="Non-positive determinant"):
            Rotation.from_matrix(mat)
    elif is_array_api_strict(xp):
        with pytest.raises(ValueError, match="SVD did not converge"):
            Rotation.from_matrix(mat)
    elif is_torch(xp):
        with pytest.raises(Exception, match="linalg.svd: "):
            Rotation.from_matrix(mat)
    else:
        assert xp.all(xp.isnan(Rotation.from_matrix(mat).as_matrix()))

    mat = xpx.at(mat)[0, 0].set(-1)
    if is_numpy(xp):
        with pytest.raises(ValueError, match="Non-positive determinant"):
            Rotation.from_matrix(mat)
    elif is_array_api_strict(xp):
        with pytest.raises(ValueError, match="SVD did not converge"):
            Rotation.from_matrix(mat)
    elif is_torch(xp):
        with pytest.raises(Exception, match="linalg.svd: "):
            Rotation.from_matrix(mat)
    else:
        assert xp.all(xp.isnan(Rotation.from_matrix(mat).as_matrix()))


def test_from_matrix_jax_compile():
    pytest.importorskip("jax")
    import jax

    m = jax.numpy.eye(3)
    r = jax.block_until_ready(jax.jit(Rotation.from_matrix)(m))
    assert isinstance(r, Rotation)


def test_from_1d_single_rotvec(xp):
    rotvec = xp.asarray([1, 0, 0])
    expected_quat = xp.asarray([0.4794255, 0, 0, 0.8775826])
    result = Rotation.from_rotvec(rotvec)
    assert_array_almost_equal(result.as_quat(), expected_quat)


def test_from_2d_single_rotvec(xp):
    rotvec = xp.asarray([[1, 0, 0]])
    expected_quat = xp.asarray([[0.4794255, 0, 0, 0.8775826]])
    result = Rotation.from_rotvec(rotvec)
    assert_array_almost_equal(result.as_quat(), expected_quat)


def test_from_generic_rotvec(xp):
    rotvec = xp.asarray([[1, 2, 2], [1, -1, 0.5], [0, 0, 0]])
    expected_quat = xp.asarray(
        [
            [0.3324983, 0.6649967, 0.6649967, 0.0707372],
            [0.4544258, -0.4544258, 0.2272129, 0.7316889],
            [0, 0, 0, 1],
        ]
    )
    assert_array_almost_equal(Rotation.from_rotvec(rotvec).as_quat(), expected_quat)


def test_from_rotvec_small_angle(xp):
    rotvec = xp.asarray(
        [
            [5e-4 / np.sqrt(3), -5e-4 / np.sqrt(3), 5e-4 / np.sqrt(3)],
            [0.2, 0.3, 0.4],
            [0, 0, 0],
        ]
    )

    quat = Rotation.from_rotvec(rotvec).as_quat()
    # cos(theta/2) ~~ 1 for small theta
    assert_allclose(quat[0, 3], 1)
    # sin(theta/2) / theta ~~ 0.5 for small theta
    assert_allclose(quat[0, :3], rotvec[0, ...] * 0.5)

    assert_allclose(quat[1, 3], 0.9639685)
    assert_allclose(
        quat[1, :3],
        xp.asarray([0.09879603932153465, 0.14819405898230198, 0.19759207864306931]),
    )

    xp_assert_equal(quat[2, ...], xp.asarray([0.0, 0, 0, 1]))


def test_degrees_from_rotvec(xp):
    rotvec1 = xp.asarray([1.0 / np.cbrt(3), 1.0 / np.cbrt(3), 1.0 / np.cbrt(3)])
    rot1 = Rotation.from_rotvec(rotvec1, degrees=True)
    quat1 = rot1.as_quat()

    # deg2rad is not implemented in Array API -> / 180 * np.pi
    rotvec2 = xp.asarray(rotvec1 / 180 * np.pi)
    rot2 = Rotation.from_rotvec(rotvec2)
    quat2 = rot2.as_quat()

    assert_allclose(quat1, quat2)


def test_malformed_1d_from_rotvec(xp):
    with pytest.raises(ValueError, match="Expected `rot_vec` to have shape"):
        Rotation.from_rotvec(xp.asarray([1, 2]))


def test_malformed_2d_from_rotvec(xp):
    with pytest.raises(ValueError, match="Expected `rot_vec` to have shape"):
        Rotation.from_rotvec(xp.asarray([[1, 2, 3, 4], [5, 6, 7, 8]]))


def test_from_rotvec_jax_compile():
    pytest.importorskip("jax")
    import jax

    rot_vec = jax.numpy.array([1, 2, 3])
    r = jax.block_until_ready(jax.jit(Rotation.from_rotvec)(rot_vec))
    assert isinstance(r, Rotation)


def test_as_generic_rotvec(xp):
    quat = xp.asarray([[1, 2, -1, 0.5], [1, -1, 1, 0.0003], [0, 0, 0, 1]])
    quat /= xp.linalg.vector_norm(quat, axis=-1, keepdims=True)

    rotvec = Rotation.from_quat(quat).as_rotvec()
    angle = xp.linalg.vector_norm(rotvec, axis=-1)

    assert_allclose(quat[:, 3], xp.cos(angle / 2))
    assert_allclose(xp.linalg.cross(rotvec, quat[:, :3]), xp.zeros((3, 3)), atol=1e-15)


def test_as_rotvec_single_1d_input(xp):
    quat = xp.asarray([1, 2, -3, 2])
    expected_rotvec = xp.asarray([0.5772381, 1.1544763, -1.7317144])

    actual_rotvec = Rotation.from_quat(quat).as_rotvec()

    assert_equal(actual_rotvec.shape, (3,))
    assert_allclose(actual_rotvec, expected_rotvec)


def test_as_rotvec_single_2d_input(xp):
    quat = xp.asarray([[1, 2, -3, 2]])
    expected_rotvec = xp.asarray([[0.5772381, 1.1544763, -1.7317144]])

    actual_rotvec = Rotation.from_quat(quat).as_rotvec()

    assert_equal(actual_rotvec.shape, (1, 3))
    assert_allclose(actual_rotvec, expected_rotvec)


def test_as_rotvec_degrees(xp):
    # x->y, y->z, z->x
    mat = xp.asarray([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    rot = Rotation.from_matrix(mat)
    rotvec = rot.as_rotvec(degrees=True)
    angle = xp.linalg.vector_norm(rotvec, axis=-1)
    assert_allclose(angle, 120.0)
    assert_allclose(rotvec[0], rotvec[1])
    assert_allclose(rotvec[1], rotvec[2])


def test_rotvec_calc_pipeline(xp):
    # Include small angles
    rotvec = xp.asarray([[0, 0, 0], [1, -1, 2], [-3e-4, 3.5e-4, 7.5e-5]])
    assert_allclose(Rotation.from_rotvec(rotvec).as_rotvec(), rotvec)
    assert_allclose(
        Rotation.from_rotvec(rotvec, degrees=True).as_rotvec(degrees=True), rotvec
    )


def test_as_rotvec_jax_compile():
    pytest.importorskip("jax")
    import jax

    r = Rotation.from_matrix(jax.numpy.eye(3))
    jax.block_until_ready(jax.jit(Rotation.as_rotvec)(r))


def test_from_1d_single_mrp(xp):
    mrp = xp.asarray([0, 0, 1.0])
    expected_quat = xp.asarray([0, 0, 1, 0])
    result = Rotation.from_mrp(mrp)
    assert_array_almost_equal(result.as_quat(), expected_quat)


def test_from_2d_single_mrp(xp):
    mrp = xp.asarray([[0, 0, 1.0]])
    expected_quat = xp.asarray([[0, 0, 1, 0]])
    result = Rotation.from_mrp(mrp)
    assert_array_almost_equal(result.as_quat(), expected_quat)


def test_from_generic_mrp(xp):
    mrp = xp.asarray([[1, 2, 2], [1, -1, 0.5], [0, 0, 0]])
    expected_quat = xp.asarray(
        [
            [0.2, 0.4, 0.4, -0.8],
            [0.61538462, -0.61538462, 0.30769231, -0.38461538],
            [0, 0, 0, 1],
        ]
    )
    assert_array_almost_equal(Rotation.from_mrp(mrp).as_quat(), expected_quat)


def test_malformed_1d_from_mrp(xp):
    with pytest.raises(ValueError, match="Expected `mrp` to have shape"):
        Rotation.from_mrp(xp.asarray([1, 2]))


def test_malformed_2d_from_mrp(xp):
    with pytest.raises(ValueError, match="Expected `mrp` to have shape"):
        Rotation.from_mrp(xp.asarray([[1, 2, 3, 4], [5, 6, 7, 8]]))


def test_from_mrp_jax_compile():
    pytest.importorskip("jax")
    import jax

    mrp = jax.numpy.array([0, 0, 1.0])
    r = jax.block_until_ready(jax.jit(Rotation.from_mrp)(mrp))
    assert isinstance(r, Rotation)


def test_as_generic_mrp(xp):
    quat = xp.asarray([[1, 2, -1, 0.5], [1, -1, 1, 0.0003], [0, 0, 0, 1]])
    quat /= xp.linalg.vector_norm(quat, axis=1)[:, None]

    expected_mrp = xp.asarray(
        [
            [0.33333333, 0.66666667, -0.33333333],
            [0.57725028, -0.57725028, 0.57725028],
            [0, 0, 0],
        ]
    )
    assert_array_almost_equal(Rotation.from_quat(quat).as_mrp(), expected_mrp)


def test_past_180_degree_rotation(xp):
    # ensure that a > 180 degree rotation is returned as a <180 rotation in MRPs
    # in this case 270 should be returned as -90
    expected_mrp = xp.asarray([-np.tan(np.pi / 2 / 4), 0.0, 0])
    assert_array_almost_equal(
        Rotation.from_euler("xyz", [270, 0, 0], degrees=True).as_mrp(), expected_mrp
    )


def test_as_mrp_single_1d_input(xp):
    quat = xp.asarray([1, 2, -3, 2])
    expected_mrp = xp.asarray([0.16018862, 0.32037724, -0.48056586])

    actual_mrp = Rotation.from_quat(quat).as_mrp()

    assert_equal(actual_mrp.shape, (3,))
    assert_allclose(actual_mrp, expected_mrp)


def test_as_mrp_single_2d_input(xp):
    quat = xp.asarray([[1, 2, -3, 2]])
    expected_mrp = xp.asarray([[0.16018862, 0.32037724, -0.48056586]])

    actual_mrp = Rotation.from_quat(quat).as_mrp()

    assert_equal(actual_mrp.shape, (1, 3))
    assert_allclose(actual_mrp, expected_mrp)


def test_mrp_calc_pipeline(xp):
    actual_mrp = xp.asarray(
        [[0, 0, 0], [1, -1, 2], [0.41421356, 0, 0], [0.1, 0.2, 0.1]]
    )
    expected_mrp = xp.asarray(
        [
            [0, 0, 0],
            [-0.16666667, 0.16666667, -0.33333333],
            [0.41421356, 0, 0],
            [0.1, 0.2, 0.1],
        ]
    )
    assert_allclose(Rotation.from_mrp(actual_mrp).as_mrp(), expected_mrp)


def test_as_mrp_jax_compile():
    pytest.importorskip("jax")
    import jax

    r = Rotation.from_matrix(jax.numpy.eye(3))
    jax.block_until_ready(jax.jit(Rotation.as_mrp)(r))


def test_from_euler_single_rotation(xp):
    angles = xp.asarray(90)
    quat = Rotation.from_euler("z", angles, degrees=True).as_quat()
    expected_quat = xp.asarray([0.0, 0, 1, 1]) / np.sqrt(2)
    assert_allclose(quat, expected_quat)


def test_single_intrinsic_extrinsic_rotation(xp):
    angles = xp.asarray(90)
    extrinsic = Rotation.from_euler("z", angles, degrees=True).as_matrix()
    intrinsic = Rotation.from_euler("Z", angles, degrees=True).as_matrix()
    assert_allclose(extrinsic, intrinsic)


def test_from_euler_rotation_order(xp):
    # Intrinsic rotation is same as extrinsic with order reversed
    rnd = np.random.RandomState(0)
    a = xp.asarray(rnd.randint(low=0, high=180, size=(6, 3)))
    b = xp.flip(a, axis=-1)
    x = Rotation.from_euler("xyz", a, degrees=True).as_quat()
    y = Rotation.from_euler("ZYX", b, degrees=True).as_quat()
    assert_allclose(x, y)


def test_from_euler_elementary_extrinsic_rotation(xp):
    # Simple test to check if extrinsic rotations are implemented correctly
    angles = xp.asarray([90, 90])
    mat = Rotation.from_euler("zx", angles, degrees=True).as_matrix()
    expected_mat = np.array([[0, -1, 0], [0, 0, -1], [1, 0, 0]])
    assert_array_almost_equal(mat, expected_mat)


def test_from_euler_intrinsic_rotation_312(xp):
    angles = xp.asarray([[30, 60, 45], [30, 60, 30], [45, 30, 60]])
    mat = Rotation.from_euler("ZXY", angles, degrees=True).as_matrix()

    assert_array_almost_equal(
        mat[0, ...],
        xp.asarray(
            [
                [0.3061862, -0.2500000, 0.9185587],
                [0.8838835, 0.4330127, -0.1767767],
                [-0.3535534, 0.8660254, 0.3535534],
            ]
        ),
    )

    assert_array_almost_equal(
        mat[1, ...],
        xp.asarray(
            [
                [0.5334936, -0.2500000, 0.8080127],
                [0.8080127, 0.4330127, -0.3995191],
                [-0.2500000, 0.8660254, 0.4330127],
            ]
        ),
    )

    assert_array_almost_equal(
        mat[2, ...],
        xp.asarray(
            [
                [0.0473672, -0.6123725, 0.7891491],
                [0.6597396, 0.6123725, 0.4355958],
                [-0.7500000, 0.5000000, 0.4330127],
            ]
        ),
    )


def test_from_euler_intrinsic_rotation_313(xp):
    angles = xp.asarray([[30, 60, 45], [30, 60, 30], [45, 30, 60]])
    mat = Rotation.from_euler("ZXZ", angles, degrees=True).as_matrix()

    assert_array_almost_equal(
        mat[0, ...],
        xp.asarray(
            [
                [0.43559574, -0.78914913, 0.4330127],
                [0.65973961, -0.04736717, -0.750000],
                [0.61237244, 0.61237244, 0.500000],
            ]
        ),
    )

    assert_array_almost_equal(
        mat[1, ...],
        xp.asarray(
            [
                [0.6250000, -0.64951905, 0.4330127],
                [0.64951905, 0.1250000, -0.750000],
                [0.4330127, 0.750000, 0.500000],
            ]
        ),
    )

    assert_array_almost_equal(
        mat[2, ...],
        xp.asarray(
            [
                [-0.1767767, -0.91855865, 0.35355339],
                [0.88388348, -0.30618622, -0.35355339],
                [0.4330127, 0.25000000, 0.8660254],
            ]
        ),
    )


def test_from_euler_extrinsic_rotation_312(xp):
    angles = xp.asarray([[30, 60, 45], [30, 60, 30], [45, 30, 60]])
    mat = Rotation.from_euler("zxy", angles, degrees=True).as_matrix()

    assert_array_almost_equal(
        mat[0, ...],
        xp.asarray(
            [
                [0.91855865, 0.1767767, 0.35355339],
                [0.25000000, 0.4330127, -0.8660254],
                [-0.30618622, 0.88388348, 0.35355339],
            ]
        ),
    )

    assert_array_almost_equal(
        mat[1, ...],
        xp.asarray(
            [
                [0.96650635, -0.0580127, 0.2500000],
                [0.25000000, 0.4330127, -0.8660254],
                [-0.0580127, 0.89951905, 0.4330127],
            ]
        ),
    )

    assert_array_almost_equal(
        mat[2, ...],
        xp.asarray(
            [
                [0.65973961, -0.04736717, 0.7500000],
                [0.61237244, 0.61237244, -0.5000000],
                [-0.43559574, 0.78914913, 0.4330127],
            ]
        ),
    )


def test_from_euler_extrinsic_rotation_313(xp):
    angles = xp.asarray([[30, 60, 45], [30, 60, 30], [45, 30, 60]])
    mat = Rotation.from_euler("zxz", angles, degrees=True).as_matrix()

    assert_array_almost_equal(
        mat[0, ...],
        xp.asarray(
            [
                [0.43559574, -0.65973961, 0.61237244],
                [0.78914913, -0.04736717, -0.61237244],
                [0.4330127, 0.75000000, 0.500000],
            ]
        ),
    )

    assert_array_almost_equal(
        mat[1, ...],
        xp.asarray(
            [
                [0.62500000, -0.64951905, 0.4330127],
                [0.64951905, 0.12500000, -0.750000],
                [0.4330127, 0.75000000, 0.500000],
            ]
        ),
    )

    assert_array_almost_equal(
        mat[2, ...],
        xp.asarray(
            [
                [-0.1767767, -0.88388348, 0.4330127],
                [0.91855865, -0.30618622, -0.250000],
                [0.35355339, 0.35355339, 0.8660254],
            ]
        ),
    )


def test_from_euler_jax_compile():
    pytest.importorskip("jax")
    import jax

    angle = jax.numpy.array([0, 0, 0])
    from_euler = jax.jit(Rotation.from_euler, static_argnums=0)
    r = jax.block_until_ready(from_euler("zxz", angle, degrees=True))
    assert isinstance(r, Rotation)


@pytest.mark.parametrize("seq_tuple", permutations("xyz"))
@pytest.mark.parametrize("intrinsic", (False, True))
def test_as_euler_asymmetric_axes(xp, seq_tuple, intrinsic):
    # helper function for mean error tests
    def test_stats(error, mean_max, rms_max):
        mean = xp.mean(error, axis=0)
        std = xp.std(error, axis=0)
        rms = xp.hypot(mean, std)
        assert xp.all(xp.abs(mean) < mean_max)
        assert xp.all(rms < rms_max)

    rnd = np.random.RandomState(0)
    n = 1000
    angles = np.empty((n, 3))
    angles[:, 0] = rnd.uniform(low=-np.pi, high=np.pi, size=(n,))
    angles[:, 1] = rnd.uniform(low=-np.pi / 2, high=np.pi / 2, size=(n,))
    angles[:, 2] = rnd.uniform(low=-np.pi, high=np.pi, size=(n,))
    angles = xp.asarray(angles)

    seq = "".join(seq_tuple)
    if intrinsic:
        # Extrinsic rotation (wrt to global world) at lower case
        # intrinsic (WRT the object itself) lower case.
        seq = seq.upper()
    rotation = Rotation.from_euler(seq, angles)
    angles_quat = rotation.as_euler(seq)
    assert_allclose(angles, angles_quat, atol=0, rtol=1e-12)
    test_stats(angles_quat - angles, 1e-15, 1e-14)


@pytest.mark.parametrize("seq_tuple", permutations("xyz"))
@pytest.mark.parametrize("intrinsic", (False, True))
def test_as_euler_symmetric_axes(xp, seq_tuple, intrinsic):
    # helper function for mean error tests
    def test_stats(error, mean_max, rms_max):
        mean = xp.mean(error, axis=0)
        std = xp.std(error, axis=0)
        rms = xp.hypot(mean, std)
        assert xp.all(xp.abs(mean) < mean_max)
        assert xp.all(rms < rms_max)

    rnd = np.random.RandomState(0)
    n = 1000
    angles = np.empty((n, 3))
    angles[:, 0] = rnd.uniform(low=-np.pi, high=np.pi, size=(n,))
    angles[:, 1] = rnd.uniform(low=0, high=np.pi, size=(n,))
    angles[:, 2] = rnd.uniform(low=-np.pi, high=np.pi, size=(n,))
    angles = xp.asarray(angles)

    # Rotation of the form A/B/A are rotation around symmetric axes
    seq = "".join([seq_tuple[0], seq_tuple[1], seq_tuple[0]])
    if intrinsic:
        seq = seq.upper()
    rotation = Rotation.from_euler(seq, angles)
    angles_quat = rotation.as_euler(seq)
    assert_allclose(angles, angles_quat, atol=0, rtol=1e-13)
    test_stats(angles_quat - angles, 1e-16, 1e-14)


@pytest.mark.thread_unsafe
@pytest.mark.parametrize("seq_tuple", permutations("xyz"))
@pytest.mark.parametrize("intrinsic", (False, True))
def test_as_euler_degenerate_asymmetric_axes(xp, seq_tuple, intrinsic):
    # Since we cannot check for angle equality, we check for rotation matrix
    # equality
    angles = xp.asarray([[45, 90, 35], [35, -90, 20], [35, 90, 25], [25, -90, 15]])

    seq = "".join(seq_tuple)
    if intrinsic:
        # Extrinsic rotation (wrt to global world) at lower case
        # Intrinsic (WRT the object itself) upper case.
        seq = seq.upper()
    rotation = Rotation.from_euler(seq, angles, degrees=True)
    mat_expected = rotation.as_matrix()

    # We cannot warn on non-NumPy backends because we'd need to condition on traced booleans
    if is_numpy(xp):
        with pytest.warns(UserWarning, match="Gimbal lock"):
            angle_estimates = rotation.as_euler(seq, degrees=True)
    else:
        angle_estimates = rotation.as_euler(seq, degrees=True)

    mat_estimated = Rotation.from_euler(seq, angle_estimates, degrees=True).as_matrix()

    assert_array_almost_equal(mat_expected, mat_estimated)


@pytest.mark.thread_unsafe
@pytest.mark.parametrize("seq_tuple", permutations("xyz"))
@pytest.mark.parametrize("intrinsic", (False, True))
def test_as_euler_degenerate_symmetric_axes(xp, seq_tuple, intrinsic):
    # Since we cannot check for angle equality, we check for rotation matrix
    # equality
    angles = xp.asarray([[15, 0, 60], [35, 0, 75], [60, 180, 35], [15, -180, 25]])

    # Rotation of the form A/B/A are rotation around symmetric axes
    seq = "".join([seq_tuple[0], seq_tuple[1], seq_tuple[0]])
    if intrinsic:
        # Extrinsic rotation (wrt to global world) at lower case
        # Intrinsic (WRT the object itself) upper case.
        seq = seq.upper()
    rotation = Rotation.from_euler(seq, angles, degrees=True)
    mat_expected = rotation.as_matrix()

    # We cannot warn on non-NumPy backends because we'd need to condition on traced booleans
    if is_numpy(xp):
        with pytest.warns(UserWarning, match="Gimbal lock"):
            angle_estimates = rotation.as_euler(seq, degrees=True)
    else:
        angle_estimates = rotation.as_euler(seq, degrees=True)
    mat_estimated = Rotation.from_euler(seq, angle_estimates, degrees=True).as_matrix()

    assert_array_almost_equal(mat_expected, mat_estimated)


def test_as_euler_jax_compile():
    pytest.importorskip("jax")
    import jax

    angle = jax.numpy.array([0, 0, 0])
    from_euler = jax.jit(Rotation.from_euler, static_argnums=0)
    jax.block_until_ready(from_euler("zxz", angle, degrees=True))


def test_inv(xp):
    rnd = np.random.RandomState(0)
    n = 10
    # preserve use of old random_state during SPEC 7 transition
    p = Rotation.random(num=n, random_state=rnd)
    # Transform to xp Rotation. Random is only available for NumPy
    p = Rotation.from_quat(xp.asarray(p.as_quat()))
    q = p.inv()

    p_mat = p.as_matrix()
    q_mat = q.as_matrix()
    result1 = np.einsum("...ij,...jk->...ik", p_mat, q_mat)
    result2 = np.einsum("...ij,...jk->...ik", q_mat, p_mat)

    eye3d = np.empty((n, 3, 3))
    eye3d[:] = np.eye(3)

    assert_array_almost_equal(result1, eye3d)
    assert_array_almost_equal(result2, eye3d)


def test_inv_single_rotation(xp):
    rng = np.random.default_rng(146972845698875399755764481408308808739)
    p = Rotation.random(rng=rng)
    # Transform to xp Rotation. Random is only available for NumPy
    p = Rotation.from_quat(xp.asarray(p.as_quat()))
    q = p.inv()

    p_mat = p.as_matrix()
    q_mat = q.as_matrix()
    res1 = xp.linalg.matmul(p_mat, q_mat)
    res2 = xp.linalg.matmul(q_mat, p_mat)

    eye = xp.eye(3)

    assert_array_almost_equal(res1, eye)
    assert_array_almost_equal(res2, eye)

    x = Rotation.random(num=1, rng=rng)
    # Transform to xp Rotation. Random is only available for NumPy
    x = Rotation.from_quat(xp.asarray(x.as_quat()))
    y = x.inv()

    x_matrix = x.as_matrix()
    y_matrix = y.as_matrix()
    result1 = xp.linalg.matmul(x_matrix, y_matrix)
    result2 = xp.linalg.matmul(y_matrix, x_matrix)

    eye3d = xp.empty((1, 3, 3))
    eye3d = xpx.at(eye3d)[..., :3, :3].set(xp.eye(3))

    assert_array_almost_equal(result1, eye3d)
    assert_array_almost_equal(result2, eye3d)


def test_inv_jax_compile():
    pytest.importorskip("jax")
    import jax

    r = Rotation.from_matrix(jax.numpy.eye(3))
    inv = jax.jit(Rotation.inv)
    jax.block_until_ready(inv(r))


def test_identity_magnitude(xp):
    n = 10
    r = Rotation.identity(n)
    r = Rotation.from_quat(xp.asarray(r.as_quat()))
    assert_allclose(r.magnitude(), 0)
    assert_allclose(r.inv().magnitude(), 0)


def test_single_identity_magnitude(xp):
    r = Rotation.from_quat(xp.asarray(Rotation.identity().as_quat()))
    assert r.magnitude() == 0
    assert r.inv().magnitude() == 0


def test_identity_invariance(xp):
    n = 10
    p = Rotation.random(n, rng=0)
    p = Rotation.from_quat(xp.asarray(p.as_quat()))

    q = Rotation.from_quat(xp.asarray(Rotation.identity(n).as_quat()))
    result = p * q
    assert_array_almost_equal(p.as_quat(), result.as_quat())

    result = result * p.inv()
    assert_array_almost_equal(result.magnitude(), xp.zeros(n))


def test_single_identity_invariance(xp):
    n = 10
    p = Rotation.random(n, rng=0)
    p = Rotation.from_quat(xp.asarray(p.as_quat()))

    q = Rotation.from_quat(xp.asarray(Rotation.identity().as_quat()))
    result = p * q
    assert_array_almost_equal(p.as_quat(), result.as_quat())

    result = result * p.inv()
    assert_array_almost_equal(result.magnitude(), xp.zeros(n))


def test_identity_jax_compile():
    pytest.importorskip("jax")
    import jax

    identity = jax.jit(Rotation.identity)
    jax.block_until_ready(identity())


def test_magnitude(xp):
    r = Rotation.from_quat(xp.eye(4))
    result = r.magnitude()
    assert_array_almost_equal(result, [np.pi, np.pi, np.pi, 0])

    r = Rotation.from_quat(-xp.eye(4))
    result = r.magnitude()
    assert_array_almost_equal(result, [np.pi, np.pi, np.pi, 0])


def test_magnitude_single_rotation(xp):
    r = Rotation.from_quat(xp.eye(4))
    result1 = r[0].magnitude()
    assert_allclose(result1, np.pi)

    result2 = r[3].magnitude()
    assert_allclose(result2, 0)


def test_magnitude_jax_compile():
    pytest.importorskip("jax")
    import jax

    magnitude = jax.jit(Rotation.magnitude)
    r = Rotation.from_matrix(jax.numpy.eye(3))
    jax.block_until_ready(magnitude(r))


def test_approx_equal(xp):
    rng = np.random.default_rng(146972845698875399755764481408308808739)
    p = Rotation.random(10, rng=rng)
    q = Rotation.random(10, rng=rng)
    # Convert random Rotation from numpy to xp
    p = Rotation.from_quat(xp.asarray(p.as_quat()))
    q = Rotation.from_quat(xp.asarray(q.as_quat()))
    r = p * q.inv()
    r_mag = r.magnitude()
    atol = xp.asarray(np.median(r_mag))  # ensure we get mix of Trues and Falses
    xp_assert_equal(p.approx_equal(q, atol), (r_mag < atol))


@pytest.mark.thread_unsafe
def test_approx_equal_single_rotation(xp):
    # also tests passing single argument to approx_equal

    p = Rotation.from_rotvec(xp.asarray([0, 0, 1e-9]))  # less than default atol of 1e-8
    q = Rotation.from_quat(xp.eye(4))
    print(q[3])
    assert p.approx_equal(q[3])
    assert not p.approx_equal(q[0])

    # test passing atol and using degrees
    assert not p.approx_equal(q[3], atol=1e-10)
    assert not p.approx_equal(q[3], atol=1e-8, degrees=True)

    # DECISION: Remove warning for degrees=True and atol=None
    if is_numpy(xp):
        with pytest.warns(UserWarning, match="atol must be set"):
            assert p.approx_equal(q[3], degrees=True)


def test_approx_equal_jax_compile():
    pytest.importorskip("jax")
    import jax

    approx_equal = jax.jit(Rotation.approx_equal)
    r = Rotation.from_matrix(jax.numpy.eye(3))
    jax.block_until_ready(approx_equal(r, r))


def test_mean(xp):
    axes = xp.concat((-xp.eye(3), xp.eye(3)))
    thetas = xp.linspace(0, np.pi / 2, 100)
    for t in thetas:
        r = Rotation.from_rotvec(t * axes)
        assert_allclose(r.mean().magnitude(), 0, atol=1e-10)


def test_weighted_mean(xp):
    # test that doubling a weight is equivalent to including a rotation twice.
    axes = xp.asarray([[0.0, 0, 0], [1, 0, 0], [1, 0, 0]])
    thetas = xp.linspace(0, np.pi / 2, 100)
    for t in thetas:
        rw = Rotation.from_rotvec(t * axes[:2, ...])
        mw = rw.mean(weights=[1, 2])

        r = Rotation.from_rotvec(t * axes)
        m = r.mean()
        assert_allclose((m * mw.inv()).magnitude(), 0, atol=1e-10)


def test_mean_invalid_weights(xp):
    # TOOD: Unify error messages and behavior
    if is_numpy(xp):
        with pytest.raises(ValueError, match="non-negative"):
            r = Rotation.from_quat(xp.eye(4))
            r.mean(weights=-xp.ones(4))
    else:
        r = Rotation.from_quat(xp.eye(4))
        if is_array_api_strict(xp):
            with pytest.raises(ValueError, match="Eigenvalues did not converge"):
                m = r.mean(weights=-xp.ones(4))
        elif is_torch(xp):
            with pytest.raises(RuntimeError, match="The algorithm failed to converge"):
                m = r.mean(weights=-xp.ones(4))
        else:
            m = r.mean(weights=-xp.ones(4))
            assert all(xp.isnan(m._quat))


def test_mean_jax_compile():
    pytest.importorskip("jax")
    import jax

    mean = jax.jit(Rotation.mean)
    r = Rotation.from_matrix(jax.numpy.eye(3))
    jax.block_until_ready(mean(r))


def test_reduction_no_indices(xp):
    r = Rotation.from_quat(xp.asarray([0.0, 0.0, 0.0, 1.0]))
    result = r.reduce(return_indices=False)
    assert isinstance(result, Rotation)


def test_reduction_none_indices(xp):
    r = Rotation.from_quat(xp.asarray([0.0, 0.0, 0.0, 1.0]))
    result = r.reduce(return_indices=True)
    assert type(result) is tuple
    assert len(result) == 3

    _, left_best, right_best = result
    assert left_best is None
    assert right_best is None


def test_reduction_scalar_calculation(xp):
    rng = np.random.default_rng(146972845698875399755764481408308808739)
    l = Rotation.random(5, rng=rng)
    r = Rotation.random(10, rng=rng)
    p = Rotation.random(7, rng=rng)
    reduced, left_best, right_best = p.reduce(l, r, return_indices=True)

    # Loop implementation of the vectorized calculation in Rotation.reduce
    scalars = np.zeros((len(l), len(p), len(r)))
    for i, li in enumerate(l):
        for j, pj in enumerate(p):
            for k, rk in enumerate(r):
                scalars[i, j, k] = np.abs((li * pj * rk).as_quat()[3])
    scalars = np.reshape(np.moveaxis(scalars, 1, 0), (scalars.shape[1], -1))

    max_ind = np.argmax(np.reshape(scalars, (len(p), -1)), axis=1)
    left_best_check = max_ind // len(r)
    right_best_check = max_ind % len(r)
    assert (left_best == left_best_check).all()
    assert (right_best == right_best_check).all()

    reduced_check = l[left_best_check] * p * r[right_best_check]
    mag = (reduced.inv() * reduced_check).magnitude()
    assert_array_almost_equal(mag, np.zeros(len(p)))


def test_reduce_jax_compile():
    pytest.importorskip("jax")
    import jax

    reduce = jax.jit(Rotation.reduce, static_argnums=(3,))
    r = Rotation.from_matrix(jax.numpy.eye(3))
    jax.block_until_ready(reduce(r, return_indices=True))


def test_apply_single_rotation_single_point(xp):
    mat = xp.asarray([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    r_1d = Rotation.from_matrix(mat)
    r_2d = Rotation.from_matrix(xp.expand_dims(mat, axis=0))

    v_1d = xp.asarray([1, 2, 3])
    v_2d = xp.expand_dims(v_1d, axis=0)
    v1d_rotated = xp.asarray([-2, 1, 3])
    v2d_rotated = xp.expand_dims(v1d_rotated, axis=0)

    assert_allclose(r_1d.apply(v_1d), v1d_rotated)
    assert_allclose(r_1d.apply(v_2d), v2d_rotated)
    assert_allclose(r_2d.apply(v_1d), v2d_rotated)
    assert_allclose(r_2d.apply(v_2d), v2d_rotated)

    v1d_inverse = xp.asarray([2, -1, 3])
    v2d_inverse = xp.expand_dims(v1d_inverse, axis=0)

    assert_allclose(r_1d.apply(v_1d, inverse=True), v1d_inverse)
    assert_allclose(r_1d.apply(v_2d, inverse=True), v2d_inverse)
    assert_allclose(r_2d.apply(v_1d, inverse=True), v2d_inverse)
    assert_allclose(r_2d.apply(v_2d, inverse=True), v2d_inverse)


def test_apply_single_rotation_multiple_points(xp):
    mat = xp.asarray([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    r1 = Rotation.from_matrix(mat)
    r2 = Rotation.from_matrix(xp.expand_dims(mat, axis=0))

    v = xp.asarray([[1, 2, 3], [4, 5, 6]])
    v_rotated = xp.asarray([[-2, 1, 3], [-5, 4, 6]])

    assert_allclose(r1.apply(v), v_rotated)
    assert_allclose(r2.apply(v), v_rotated)

    v_inverse = np.array([[2, -1, 3], [5, -4, 6]])

    assert_allclose(r1.apply(v, inverse=True), v_inverse)
    assert_allclose(r2.apply(v, inverse=True), v_inverse)


def test_apply_multiple_rotations_single_point(xp):
    mat = xp.empty((2, 3, 3))
    mat = xpx.at(mat)[0, ...].set(xp.asarray([[0, -1, 0], [1, 0, 0], [0, 0, 1]]))
    mat = xpx.at(mat)[1, ...].set(xp.asarray([[1, 0, 0], [0, 0, -1], [0, 1, 0]]))
    r = Rotation.from_matrix(mat)

    v1 = xp.asarray([1, 2, 3])
    v2 = xp.expand_dims(v1, axis=0)

    v_rotated = xp.asarray([[-2, 1, 3], [1, -3, 2]])

    assert_allclose(r.apply(v1), v_rotated)
    assert_allclose(r.apply(v2), v_rotated)

    v_inverse = xp.asarray([[2, -1, 3], [1, 3, -2]])

    assert_allclose(r.apply(v1, inverse=True), v_inverse)
    assert_allclose(r.apply(v2, inverse=True), v_inverse)


def test_apply_multiple_rotations_multiple_points(xp):
    mat = xp.empty((2, 3, 3))
    mat = xpx.at(mat)[0, ...].set(xp.asarray([[0, -1, 0], [1, 0, 0], [0, 0, 1]]))
    mat = xpx.at(mat)[1, ...].set(xp.asarray([[1, 0, 0], [0, 0, -1], [0, 1, 0]]))
    r = Rotation.from_matrix(mat)

    v = xp.asarray([[1, 2, 3], [4, 5, 6]])
    v_rotated = xp.asarray([[-2, 1, 3], [4, -6, 5]])
    assert_allclose(r.apply(v), v_rotated)

    v_inverse = xp.asarray([[2, -1, 3], [4, 6, -5]])
    assert_allclose(r.apply(v, inverse=True), v_inverse)


def test_apply_jax_compile():
    pytest.importorskip("jax")
    import jax

    apply = jax.jit(Rotation.apply, static_argnums=(2,))
    r = Rotation.from_matrix(jax.numpy.eye(3))
    jax.block_until_ready(apply(r, jax.numpy.eye(3), inverse=False))


def test_getitem(xp):
    mat = xp.empty((2, 3, 3))
    mat = xpx.at(mat)[0, ...].set(xp.asarray([[0, -1, 0], [1, 0, 0], [0, 0, 1]]))
    mat = xpx.at(mat)[1, ...].set(xp.asarray([[1, 0, 0], [0, 0, -1], [0, 1, 0]]))
    r = Rotation.from_matrix(mat)

    assert_allclose(r[0].as_matrix(), mat[0, ...], atol=1e-15)
    assert_allclose(r[1].as_matrix(), mat[1, ...], atol=1e-15)
    assert_allclose(r[:-1].as_matrix(), xp.expand_dims(mat[0, ...], axis=0), atol=1e-15)


def test_getitem_single(xp):
    r = Rotation.from_quat(xp.asarray([0, 0, 0, 1]))
    with pytest.raises(TypeError, match="not subscriptable"):
        r[0]


def test_getitem_jax_compile():
    pytest.importorskip("jax")
    import jax

    getitem = jax.jit(Rotation.__getitem__, static_argnums=(1,))
    r = Rotation.from_matrix(jax.numpy.eye(3).reshape(1, 3, 3))
    jax.block_until_ready(getitem(r, 0))


def test_setitem_single(xp):
    r = Rotation.from_quat(xp.asarray([0, 0, 0, 1]))
    with pytest.raises(TypeError, match="not subscriptable"):
        r[0] = Rotation.from_quat(xp.asarray([0, 0, 0, 1]))


def test_setitem_slice(xp):
    rng = np.random.default_rng(146972845698875399755764481408308808739)
    r1 = Rotation.from_quat(xp.asarray(Rotation.random(10, rng=rng).as_quat()))
    r2 = Rotation.from_quat(xp.asarray(Rotation.random(5, rng=rng).as_quat()))
    r1[1:6] = r2
    xp_assert_equal(r1[1:6].as_quat(), r2.as_quat())


def test_setitem_integer(xp):
    rng = np.random.default_rng(146972845698875399755764481408308808739)
    r1 = Rotation.from_quat(xp.asarray(Rotation.random(10, rng=rng).as_quat()))
    r2 = Rotation.from_quat(xp.asarray(Rotation.random(rng=rng).as_quat()))
    r1[1] = r2
    xp_assert_equal(r1[1].as_quat(), r2.as_quat())


def test_setitem_wrong_type(xp):
    rng = np.random.default_rng(0)
    r = Rotation.from_quat(xp.asarray(Rotation.random(10, rng=rng).as_quat()))
    with pytest.raises(TypeError, match="Rotation object"):
        r[0] = 1


def test_setitem_jax_compile():
    pytest.importorskip("jax")
    import jax

    setitem = jax.jit(Rotation.__setitem__, static_argnums=(1,))
    r = Rotation.from_matrix(jax.numpy.eye(3).reshape(1, 3, 3))
    jax.block_until_ready(setitem(r, 0, r[0]))


def test_n_rotations(xp):
    mat = xp.empty((2, 3, 3))
    mat = xpx.at(mat)[0, ...].set(xp.asarray([[0, -1, 0], [1, 0, 0], [0, 0, 1]]))
    mat = xpx.at(mat)[1, ...].set(xp.asarray([[1, 0, 0], [0, 0, -1], [0, 1, 0]]))
    r = Rotation.from_matrix(mat)

    assert_equal(len(r), 2)
    assert_equal(len(r[:-1]), 1)


def test_len_jax_compile():
    pytest.importorskip("jax")
    import jax

    r = Rotation.from_matrix(jax.numpy.eye(3).reshape(1, 3, 3))
    jax.block_until_ready(jax.jit(lambda rot: len(rot))(r))


def test_random_rotation_shape():
    # No xp testing since random rotations are always using NumPy
    rng = np.random.default_rng(146972845698875399755764481408308808739)
    assert_equal(Rotation.random(rng=rng).as_quat().shape, (4,))
    assert_equal(Rotation.random(None, rng=rng).as_quat().shape, (4,))

    assert_equal(Rotation.random(1, rng=rng).as_quat().shape, (1, 4))
    assert_equal(Rotation.random(5, rng=rng).as_quat().shape, (5, 4))


def test_align_vectors_no_rotation(xp):
    x = xp.asarray([[1, 2, 3], [4, 5, 6]])
    y = xp.asarray(x, copy=True)

    r, rssd = Rotation.align_vectors(x, y)
    assert_array_almost_equal(r.as_matrix(), np.eye(3))
    assert_allclose(rssd, 0, atol=1e-6)


def test_align_vectors_no_noise(xp):
    rng = np.random.default_rng(14697284569885399755764481408308808739)
    c = Rotation.random(rng=rng)
    c = c.from_quat(xp.asarray(c.as_quat()))
    b = xp.asarray(rng.normal(size=(5, 3)))
    a = c.apply(b)

    est, rssd = Rotation.align_vectors(a, b)
    assert_allclose(c.as_quat(), est.as_quat())
    assert_allclose(rssd, 0, atol=1e-7)


def test_align_vectors_improper_rotation(xp):
    # Tests correct logic for issue #10444
    x = xp.asarray(
        [[0.89299824, -0.44372674, 0.0752378], [0.60221789, -0.47564102, -0.6411702]]
    )
    y = xp.asarray(
        [[0.02386536, -0.82176463, 0.5693271], [-0.27654929, -0.95191427, -0.1318321]]
    )

    est, rssd = Rotation.align_vectors(x, y)
    assert_allclose(x, est.apply(y), atol=1e-6)
    assert_allclose(rssd, 0, atol=1e-7)


def test_align_vectors_rssd_sensitivity(xp):
    rssd_expected = xp.asarray(0.141421356237308)
    sens_expected = xp.asarray([[0.2, 0.0, 0.0], [0.0, 1.5, 1.0], [0.0, 1.0, 1.0]])
    atol = 1e-6
    a = xp.asarray([[0, 1, 0], [0, 1, 1], [0, 1, 1]])
    b = xp.asarray([[1, 0, 0], [1, 1.1, 0], [1, 0.9, 0]])
    rot, rssd, sens = Rotation.align_vectors(a, b, return_sensitivity=True)
    assert xp.all(xpx.isclose(rssd, rssd_expected, atol=atol, xp=xp))
    assert_allclose(sens, sens_expected, atol=atol)


def test_align_vectors_scaled_weights(xp):
    n = 10
    a = xp.asarray(Rotation.random(n, rng=0).apply([1, 0, 0]))
    b = xp.asarray(Rotation.random(n, rng=1).apply([1, 0, 0]))
    scale = xp.asarray(2.0)

    est1, rssd1, cov1 = Rotation.align_vectors(a, b, xp.ones(n), True)
    est2, rssd2, cov2 = Rotation.align_vectors(a, b, scale * xp.ones(n), True)

    assert_allclose(est1.as_matrix(), est2.as_matrix())
    assert_allclose(xp.sqrt(scale) * rssd1, rssd2, atol=1e-6)
    assert_allclose(cov1, cov2)


def test_align_vectors_noise(xp):
    rng = np.random.default_rng(146972845698875399755764481408308808739)
    n_vectors = 100
    rot = Rotation.from_quat(xp.asarray(Rotation.random(rng=rng).as_quat()))
    vectors = xp.asarray(rng.normal(size=(n_vectors, 3)))
    result = rot.apply(vectors)

    # The paper adds noise as independently distributed angular errors
    sigma = np.deg2rad(1)
    tolerance = 1.5 * sigma
    noise = Rotation.from_rotvec(
        xp.asarray(rng.normal(size=(n_vectors, 3), scale=sigma))
    )

    # Attitude errors must preserve norm. Hence apply individual random
    # rotations to each vector.
    noisy_result = noise.apply(result)

    est, rssd, cov = Rotation.align_vectors(
        noisy_result, vectors, return_sensitivity=True
    )

    # Use rotation compositions to find out closeness
    error_vector = (rot * est.inv()).as_rotvec()
    assert_allclose(error_vector[0], 0, atol=tolerance)
    assert_allclose(error_vector[1], 0, atol=tolerance)
    assert_allclose(error_vector[2], 0, atol=tolerance)

    # Check error bounds using covariance matrix
    cov *= sigma
    assert_allclose(cov[0, 0], 0, atol=tolerance)
    assert_allclose(cov[1, 1], 0, atol=tolerance)
    assert_allclose(cov[2, 2], 0, atol=tolerance)

    assert_allclose(rssd, xp.sum((noisy_result - est.apply(vectors)) ** 2) ** 0.5)


def test_align_vectors_invalid_input(xp):
    if is_numpy(xp):
        with pytest.raises(ValueError, match="Expected input `a` to have shape"):
            a, b = xp.asarray([1, 2, 3, 4]), xp.asarray([1, 2, 3])
            Rotation.align_vectors(a, b)

    with pytest.raises(ValueError, match="Expected input `b` to have shape"):
        a, b = xp.asarray([1, 2, 3]), xp.asarray([1, 2, 3, 4])
        Rotation.align_vectors(a, b)

    with pytest.raises(
        ValueError, match="Expected inputs `a` and `b` to have same shapes"
    ):
        a, b = xp.asarray([[1, 2, 3], [4, 5, 6]]), xp.asarray([[1, 2, 3]])
        Rotation.align_vectors(a, b)

    with pytest.raises(ValueError, match="Expected `weights` to be 1 dimensional"):
        a, b = xp.asarray([[1, 2, 3]]), xp.asarray([[1, 2, 3]])
        weights = xp.asarray([[1]])
        Rotation.align_vectors(a, b, weights)

    a, b = xp.asarray([[1, 2, 3], [4, 5, 6]]), xp.asarray([[1, 2, 3], [4, 5, 6]])
    weights = xp.asarray([1, 2, 3])
    with pytest.raises(ValueError, match="Expected `weights` to have number of values"):
        Rotation.align_vectors(a, b, weights)

    a, b = xp.asarray([[1, 2, 3]]), xp.asarray([[1, 2, 3]])
    weights = xp.asarray([-1])
    if is_numpy(xp):
        with pytest.raises(
            ValueError, match="`weights` may not contain negative values"
        ):
            Rotation.align_vectors(a, b, weights)
    else:
        r, rssd = Rotation.align_vectors(a, b, weights)
        assert xp.all(xp.isnan(r.as_quat())), "Quaternion should be nan"
        assert xp.isnan(rssd), "RSSD should be nan"

    a, b = xp.asarray([[1, 2, 3], [4, 5, 6]]), xp.asarray([[1, 2, 3], [4, 5, 6]])
    weights = xp.asarray([np.inf, np.inf])
    if is_numpy(xp):
        with pytest.raises(ValueError, match="Only one infinite weight is allowed"):
            Rotation.align_vectors(a, b, weights)
    elif is_array_api_strict(xp):
        with pytest.raises(ValueError, match="SVD did not converge"):
            Rotation.align_vectors(a, b, weights)
    elif is_torch(xp):
        with pytest.raises(Exception, match="linalg.svd: "):
            Rotation.align_vectors(a, b, weights)
    else:
        r, rssd = Rotation.align_vectors(a, b, weights)
        assert xp.all(xp.isnan(r.as_quat())), "Quaternion should be nan"
        assert xp.isnan(rssd), "RSSD should be nan"

    a, b = xp.asarray([[0, 0, 0]]), xp.asarray([[1, 2, 3]])
    if is_numpy(xp):
        with pytest.raises(
            ValueError, match="Cannot align zero length primary vectors"
        ):
            Rotation.align_vectors(a, b)
    else:
        r, rssd = Rotation.align_vectors(a, b)
        assert xp.all(xp.isnan(r.as_quat())), "Quaternion should be nan"
        assert xp.isnan(rssd), "RSSD should be nan"

    a, b = (xp.asarray([[1, 2, 3], [4, 5, 6]]), xp.asarray([[1, 2, 3], [4, 5, 6]]))
    weights = xp.asarray([np.inf, 1])
    if is_numpy(xp):
        with pytest.raises(ValueError, match="Cannot return sensitivity matrix"):
            Rotation.align_vectors(a, b, weights, return_sensitivity=True)
    else:
        r, rssd, sens = Rotation.align_vectors(a, b, weights, return_sensitivity=True)
        assert sens is None, "Sensitivity matrix should be None"

    a, b = xp.asarray([[1, 2, 3]]), xp.asarray([[1, 2, 3]])
    if is_numpy(xp):
        with pytest.raises(ValueError, match="Cannot return sensitivity matrix"):
            Rotation.align_vectors(a, b, return_sensitivity=True)
    else:
        r, rssd, sens = Rotation.align_vectors(a, b, return_sensitivity=True)
        assert sens is None, "Sensitivity matrix should be None"


def test_align_vectors_align_constrain(xp):
    # Align the primary +X B axis with the primary +Y A axis, and rotate about
    # it such that the +Y B axis (residual of the [1, 1, 0] secondary b vector)
    # is aligned with the +Z A axis (residual of the [0, 1, 1] secondary a
    # vector)
    atol = 1e-12
    b = xp.asarray([[1, 0, 0], [1, 1, 0]])
    a = xp.asarray([[0, 1, 0], [0, 1, 1]])
    m_expected = xp.asarray([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    R, rssd = Rotation.align_vectors(a, b, weights=xp.asarray([np.inf, 1]))
    assert_allclose(R.as_matrix(), m_expected, atol=atol)
    assert_allclose(R.apply(b), a, atol=atol)  # Pri and sec align exactly
    assert xpx.isclose(rssd, xp.asarray(0.0), atol=atol, xp=xp)

    # Do the same but with an inexact secondary rotation
    b = xp.asarray([[1, 0, 0], [1, 2, 0]])
    rssd_expected = xp.asarray(1.0)
    R, rssd = Rotation.align_vectors(a, b, weights=xp.asarray([np.inf, 1]))
    assert_allclose(R.as_matrix(), m_expected, atol=atol)
    assert_allclose(R.apply(b)[0, ...], a[0, ...], atol=atol)  # Only pri aligns exactly
    assert xpx.isclose(rssd, rssd_expected, atol=atol, xp=xp)
    a_expected = xp.asarray([[0, 1, 0], [0, 1, 2]])
    assert_allclose(R.apply(b), a_expected, atol=atol)

    # Check random vectors
    b = xp.asarray([[1, 2, 3], [-2, 3, -1]])
    a = xp.asarray([[-1, 3, 2], [1, -1, 2]])
    rssd_expected = xp.asarray(1.3101595297515016)
    R, rssd = Rotation.align_vectors(a, b, weights=xp.asarray([np.inf, 1]))
    assert_allclose(R.apply(b)[0, ...], a[0, ...], atol=atol)  # Only pri aligns exactly
    assert xpx.isclose(rssd, rssd_expected, atol=atol, xp=xp)


def test_align_vectors_near_inf(xp):
    # align_vectors should return near the same result for high weights as for
    # infinite weights. rssd will be different with floating point error on the
    # exactly aligned vector being multiplied by a large non-infinite weight

    # TODO: Consider jitting for JAX. Non-jitted version is quite slow
    n = 100
    mats = []
    for i in range(6):
        mats.append(Rotation.random(n, rng=10 + i).as_matrix())

    for i in range(n):
        # Get random pairs of 3-element vectors
        # Creating tensors from list of numpy arrays fails in PyTorch
        a = xp.asarray(np.array([1 * mats[0][i][0], 2 * mats[1][i][0]]))
        b = xp.asarray(np.array([3 * mats[2][i][0], 4 * mats[3][i][0]]))

        R, _ = Rotation.align_vectors(a, b, weights=[1e10, 1])
        R2, _ = Rotation.align_vectors(a, b, weights=[np.inf, 1])
        assert_allclose(R.as_matrix(), R2.as_matrix(), atol=1e-4)

    for i in range(n):
        # Get random triplets of 3-element vectors
        a = np.array([1 * mats[0][i][0], 2 * mats[1][i][0], 3 * mats[2][i][0]])
        b = np.array([4 * mats[3][i][0], 5 * mats[4][i][0], 6 * mats[5][i][0]])
        a = xp.asarray(a)
        b = xp.asarray(b)

        R, _ = Rotation.align_vectors(a, b, weights=[1e10, 2, 1])
        R2, _ = Rotation.align_vectors(a, b, weights=[np.inf, 2, 1])
        assert_allclose(R.as_matrix(), R2.as_matrix(), atol=1e-4)


def test_align_vectors_parallel(xp):
    atol = 1e-12
    a = xp.asarray([[1, 0, 0], [0, 1, 0]])
    b = xp.asarray([[0, 1, 0], [0, 1, 0]])
    m_expected = xp.asarray([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
    R, _ = Rotation.align_vectors(a, b, weights=[np.inf, 1])
    assert_allclose(R.as_matrix(), m_expected, atol=atol)
    R, _ = Rotation.align_vectors(a[0, ...], b[0, ...])
    assert_allclose(R.as_matrix(), m_expected, atol=atol)
    assert_allclose(R.apply(b[0, ...]), a[0, ...], atol=atol)

    b = xp.asarray([[1, 0, 0], [1, 0, 0]])
    m_expected = xp.asarray([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    R, _ = Rotation.align_vectors(a, b, weights=[np.inf, 1])
    assert_allclose(R.as_matrix(), m_expected, atol=atol)
    R, _ = Rotation.align_vectors(a[0, ...], b[0, ...])
    assert_allclose(R.as_matrix(), m_expected, atol=atol)
    assert_allclose(R.apply(b[0, ...]), a[0, ...], atol=atol)


def test_align_vectors_antiparallel(xp):
    # Test exact 180 deg rotation
    atol = 1e-12
    as_to_test = np.array(
        [[[1, 0, 0], [0, 1, 0]], [[0, 1, 0], [1, 0, 0]], [[0, 0, 1], [0, 1, 0]]]
    )

    bs_to_test = np.array([[-a[0], a[1]] for a in as_to_test])
    for a, b in zip(as_to_test, bs_to_test):
        a, b = xp.asarray(a), xp.asarray(b)
        R, _ = Rotation.align_vectors(a, b, weights=[np.inf, 1])
        assert_allclose(R.magnitude(), np.pi, atol=atol)
        assert_allclose(R.apply(b[0, ...]), a[0, ...], atol=atol)

    # Test exact rotations near 180 deg
    Rs = Rotation.random(100, rng=0)
    dRs = Rotation.from_rotvec(Rs.as_rotvec() * 1e-4)  # scale down to small angle
    a = [[1, 0, 0], [0, 1, 0]]
    b = [[-1, 0, 0], [0, 1, 0]]
    as_to_test = []
    for dR in dRs:
        as_to_test.append(np.array([dR.apply(a[0]), a[1]]))
    for a in as_to_test:
        a, b = xp.asarray(a), xp.asarray(b)
        R, _ = Rotation.align_vectors(a, b, weights=[np.inf, 1])
        R2, _ = Rotation.align_vectors(a, b, weights=[1e10, 1])
        assert_allclose(R.as_matrix(), R2.as_matrix(), atol=atol)


def test_align_vectors_primary_only(xp):
    atol = 1e-12
    mats_a = Rotation.random(100, rng=0).as_matrix()
    mats_b = Rotation.random(100, rng=1).as_matrix()
    for mat_a, mat_b in zip(mats_a, mats_b):
        # Get random 3-element unit vectors
        a = xp.asarray(mat_a[0])
        b = xp.asarray(mat_b[0])

        # Compare to align_vectors with primary only
        R, rssd = Rotation.align_vectors(a, b)
        assert_allclose(R.apply(b), a, atol=atol)
        assert np.isclose(rssd, 0, atol=atol)


def test_align_vectors_jax_compile():
    pytest.importorskip("jax")
    import jax

    rng = np.random.default_rng(14697284569885399755764481408308808739)
    c = Rotation.random(rng=rng)
    c = c.from_quat(jax.numpy.asarray(c.as_quat()))
    b = jax.numpy.asarray(rng.normal(size=(5, 3)))
    a = c.apply(b)

    est, rssd = jax.block_until_ready(jax.jit(Rotation.align_vectors)(a, b))
    assert_allclose(c.as_quat(), est.as_quat())
    assert_allclose(rssd, 0, atol=1e-7)


def test_repr_single_rotation(xp):
    q = xp.asarray([0, 0, 0, 1])
    actual = repr(Rotation.from_quat(q))
    expected = """\
Rotation.from_matrix(array([[1., 0., 0.],
                            [0., 1., 0.],
                            [0., 0., 1.]]))"""
    assert actual == expected


def test_repr_rotation_sequence(xp):
    q = xp.asarray([[0.0, 1, 0, 1], [0, 0, 1, 1]]) / np.sqrt(2)
    actual = f"{Rotation.from_quat(q)!r}"
    expected = """\
Rotation.from_matrix(array([[[ 0.,  0.,  1.],
                             [ 0.,  1.,  0.],
                             [-1.,  0.,  0.]],
                     
                            [[ 0., -1.,  0.],
                             [ 1.,  0.,  0.],
                             [ 0.,  0.,  1.]]]))"""
    assert actual == expected


def test_slerp():
    rnd = np.random.RandomState(0)

    key_rots = Rotation.from_quat(rnd.uniform(size=(5, 4)))
    key_quats = key_rots.as_quat()

    key_times = [0, 1, 2, 3, 4]
    interpolator = Slerp(key_times, key_rots)

    times = [0, 0.5, 0.25, 1, 1.5, 2, 2.75, 3, 3.25, 3.60, 4]
    interp_rots = interpolator(times)
    interp_quats = interp_rots.as_quat()

    # Dot products are affected by sign of quaternions
    interp_quats[interp_quats[:, -1] < 0] *= -1
    # Checking for quaternion equality, perform same operation
    key_quats[key_quats[:, -1] < 0] *= -1

    # Equality at keyframes, including both endpoints
    assert_allclose(interp_quats[0], key_quats[0])
    assert_allclose(interp_quats[3], key_quats[1])
    assert_allclose(interp_quats[5], key_quats[2])
    assert_allclose(interp_quats[7], key_quats[3])
    assert_allclose(interp_quats[10], key_quats[4])

    # Constant angular velocity between keyframes. Check by equating
    # cos(theta) between quaternion pairs with equal time difference.
    cos_theta1 = np.sum(interp_quats[0] * interp_quats[2])
    cos_theta2 = np.sum(interp_quats[2] * interp_quats[1])
    assert_allclose(cos_theta1, cos_theta2)

    cos_theta4 = np.sum(interp_quats[3] * interp_quats[4])
    cos_theta5 = np.sum(interp_quats[4] * interp_quats[5])
    assert_allclose(cos_theta4, cos_theta5)

    # theta1: 0 -> 0.25, theta3 : 0.5 -> 1
    # Use double angle formula for double the time difference
    cos_theta3 = np.sum(interp_quats[1] * interp_quats[3])
    assert_allclose(cos_theta3, 2 * (cos_theta1**2) - 1)

    # Miscellaneous checks
    assert_equal(len(interp_rots), len(times))


def test_slerp_rot_is_rotation():
    with pytest.raises(TypeError, match="must be a `Rotation` instance"):
        r = np.array([[1, 2, 3, 4], [0, 0, 0, 1]])
        t = np.array([0, 1])
        Slerp(t, r)


def test_slerp_single_rot():
    msg = "must be a sequence of at least 2 rotations"
    with pytest.raises(ValueError, match=msg):
        r = Rotation.from_quat([1, 2, 3, 4])
        Slerp([1], r)


def test_slerp_rot_len1():
    msg = "must be a sequence of at least 2 rotations"
    with pytest.raises(ValueError, match=msg):
        r = Rotation.from_quat([[1, 2, 3, 4]])
        Slerp([1], r)


def test_slerp_time_dim_mismatch():
    with pytest.raises(
        ValueError, match="times to be specified in a 1 dimensional array"
    ):
        rnd = np.random.RandomState(0)
        r = Rotation.from_quat(rnd.uniform(size=(2, 4)))
        t = np.array([[1], [2]])
        Slerp(t, r)


def test_slerp_num_rotations_mismatch():
    with pytest.raises(
        ValueError, match="number of rotations to be equal to number of timestamps"
    ):
        rnd = np.random.RandomState(0)
        r = Rotation.from_quat(rnd.uniform(size=(5, 4)))
        t = np.arange(7)
        Slerp(t, r)


def test_slerp_equal_times():
    with pytest.raises(ValueError, match="strictly increasing order"):
        rnd = np.random.RandomState(0)
        r = Rotation.from_quat(rnd.uniform(size=(5, 4)))
        t = [0, 1, 2, 2, 4]
        Slerp(t, r)


def test_slerp_decreasing_times():
    with pytest.raises(ValueError, match="strictly increasing order"):
        rnd = np.random.RandomState(0)
        r = Rotation.from_quat(rnd.uniform(size=(5, 4)))
        t = [0, 1, 3, 2, 4]
        Slerp(t, r)


def test_slerp_call_time_dim_mismatch():
    rnd = np.random.RandomState(0)
    r = Rotation.from_quat(rnd.uniform(size=(5, 4)))
    t = np.arange(5)
    s = Slerp(t, r)

    with pytest.raises(ValueError, match="`times` must be at most 1-dimensional."):
        interp_times = np.array([[3.5], [4.2]])
        s(interp_times)


def test_slerp_call_time_out_of_range():
    rnd = np.random.RandomState(0)
    r = Rotation.from_quat(rnd.uniform(size=(5, 4)))
    t = np.arange(5) + 1
    s = Slerp(t, r)

    with pytest.raises(ValueError, match="times must be within the range"):
        s([0, 1, 2])
    with pytest.raises(ValueError, match="times must be within the range"):
        s([1, 2, 6])


def test_slerp_call_scalar_time():
    r = Rotation.from_euler("X", [0, 80], degrees=True)
    s = Slerp([0, 1], r)

    r_interpolated = s(0.25)
    r_interpolated_expected = Rotation.from_euler("X", 20, degrees=True)

    delta = r_interpolated * r_interpolated_expected.inv()

    assert_allclose(delta.magnitude(), 0, atol=1e-16)


def test_multiplication_stability():
    qs = Rotation.random(50, rng=0)
    rs = Rotation.random(1000, rng=1)
    for q in qs:
        rs *= q * rs
        assert_allclose(np.linalg.norm(rs.as_quat(), axis=1), 1)


def test_pow():
    atol = 1e-14
    p = Rotation.random(10, rng=0)
    p_inv = p.inv()
    # Test the short-cuts and other integers
    for n in [-5, -2, -1, 0, 1, 2, 5]:
        # Test accuracy
        q = p**n
        r = Rotation.identity(10)
        for _ in range(abs(n)):
            if n > 0:
                r = r * p
            else:
                r = r * p_inv
        ang = (q * r.inv()).magnitude()
        assert np.all(ang < atol)

        # Test shape preservation
        r = Rotation.from_quat([0, 0, 0, 1])
        assert (r**n).as_quat().shape == (4,)
        r = Rotation.from_quat([[0, 0, 0, 1]])
        assert (r**n).as_quat().shape == (1, 4)

    # Large angle fractional
    for n in [-1.5, -0.5, -0.0, 0.0, 0.5, 1.5]:
        q = p**n
        r = Rotation.from_rotvec(n * p.as_rotvec())
        assert_allclose(q.as_quat(), r.as_quat(), atol=atol)

    # Small angle
    p = Rotation.from_rotvec([1e-12, 0, 0])
    n = 3
    q = p**n
    r = Rotation.from_rotvec(n * p.as_rotvec())
    assert_allclose(q.as_quat(), r.as_quat(), atol=atol)


def test_pow_errors():
    p = Rotation.random(rng=0)
    with pytest.raises(NotImplementedError, match="modulus not supported"):
        pow(p, 1, 1)


def test_rotation_within_numpy_array():
    single = Rotation.random(rng=0)
    multiple = Rotation.random(2, rng=1)

    array = np.array(single)
    assert_equal(array.shape, ())

    array = np.array(multiple)
    assert_equal(array.shape, (2,))
    assert_allclose(array[0].as_matrix(), multiple[0].as_matrix())
    assert_allclose(array[1].as_matrix(), multiple[1].as_matrix())

    array = np.array([single])
    assert_equal(array.shape, (1,))
    assert_equal(array[0], single)

    array = np.array([multiple])
    assert_equal(array.shape, (1, 2))
    assert_allclose(array[0, 0].as_matrix(), multiple[0].as_matrix())
    assert_allclose(array[0, 1].as_matrix(), multiple[1].as_matrix())

    array = np.array([single, multiple], dtype=object)
    assert_equal(array.shape, (2,))
    assert_equal(array[0], single)
    assert_equal(array[1], multiple)

    array = np.array([multiple, multiple, multiple])
    assert_equal(array.shape, (3, 2))


def test_pickling():
    r = Rotation.from_quat([0, 0, np.sin(np.pi / 4), np.cos(np.pi / 4)])
    pkl = pickle.dumps(r)
    unpickled = pickle.loads(pkl)
    assert_allclose(r.as_matrix(), unpickled.as_matrix(), atol=1e-15)


def test_deepcopy():
    r = Rotation.from_quat([0, 0, np.sin(np.pi / 4), np.cos(np.pi / 4)])
    r1 = copy.deepcopy(r)
    assert_allclose(r.as_matrix(), r1.as_matrix(), atol=1e-15)


def test_as_euler_contiguous():
    # The Array API does not specify contiguous arrays, so we can only check for NumPy
    r = Rotation.from_quat([0, 0, 0, 1])
    e1 = r.as_euler("xyz")  # extrinsic euler rotation
    e2 = r.as_euler("XYZ")  # intrinsic
    assert e1.flags["C_CONTIGUOUS"] is True
    assert e2.flags["C_CONTIGUOUS"] is True
    assert all(i >= 0 for i in e1.strides)
    assert all(i >= 0 for i in e2.strides)


def test_concatenate():
    rotation = Rotation.random(10, rng=0)
    sizes = [1, 2, 3, 1, 3]
    starts = [0] + list(np.cumsum(sizes))
    split = [rotation[i : i + n] for i, n in zip(starts, sizes)]
    result = Rotation.concatenate(split)
    assert_equal(rotation.as_quat(), result.as_quat())

    # Test Rotation input for multiple rotations
    result = Rotation.concatenate(rotation)
    assert_equal(rotation.as_quat(), result.as_quat())

    # Test that a copy is returned
    assert rotation is not result

    # Test Rotation input for single rotations
    result = Rotation.concatenate(Rotation.identity())
    assert_equal(Rotation.identity().as_quat(), result.as_quat())


def test_concatenate_wrong_type():
    with pytest.raises(TypeError, match="Rotation objects only"):
        Rotation.concatenate([Rotation.identity(), 1, None])


# Regression test for gh-16663
def test_len_and_bool():
    rotation_multi_one = Rotation([[0, 0, 0, 1]])
    rotation_multi = Rotation([[0, 0, 0, 1], [0, 0, 0, 1]])
    rotation_single = Rotation([0, 0, 0, 1])

    assert len(rotation_multi_one) == 1
    assert len(rotation_multi) == 2
    with pytest.raises(TypeError, match="Single rotation has no len()."):
        len(rotation_single)

    # Rotation should always be truthy. See gh-16663
    assert rotation_multi_one
    assert rotation_multi
    assert rotation_single


def test_from_davenport_single_rotation():
    axis = [0, 0, 1]
    quat = Rotation.from_davenport(axis, "extrinsic", 90, degrees=True).as_quat()
    expected_quat = np.array([0, 0, 1, 1]) / np.sqrt(2)
    assert_allclose(quat, expected_quat)


def test_from_davenport_one_or_two_axes():
    ez = [0, 0, 1]
    ey = [0, 1, 0]

    # Single rotation, single axis, axes.shape == (3, )
    rot = Rotation.from_rotvec(np.array(ez) * np.pi / 4)
    rot_dav = Rotation.from_davenport(ez, "e", np.pi / 4)
    assert_allclose(rot.as_quat(canonical=True), rot_dav.as_quat(canonical=True))

    # Single rotation, single axis, axes.shape == (1, 3)
    rot = Rotation.from_rotvec([np.array(ez) * np.pi / 4])
    rot_dav = Rotation.from_davenport([ez], "e", [np.pi / 4])
    assert_allclose(rot.as_quat(canonical=True), rot_dav.as_quat(canonical=True))

    # Single rotation, two axes, axes.shape == (2, 3)
    rot = Rotation.from_rotvec([np.array(ez) * np.pi / 4, np.array(ey) * np.pi / 6])
    rot = rot[0] * rot[1]
    rot_dav = Rotation.from_davenport([ey, ez], "e", [np.pi / 6, np.pi / 4])
    assert_allclose(rot.as_quat(canonical=True), rot_dav.as_quat(canonical=True))

    # Two rotations, single axis, axes.shape == (3, )
    rot = Rotation.from_rotvec([np.array(ez) * np.pi / 6, np.array(ez) * np.pi / 4])
    rot_dav = Rotation.from_davenport([ez], "e", [np.pi / 6, np.pi / 4])
    assert_allclose(rot.as_quat(canonical=True), rot_dav.as_quat(canonical=True))


def test_from_davenport_invalid_input():
    ez = [0, 0, 1]
    ey = [0, 1, 0]
    ezy = [0, 1, 1]
    with pytest.raises(ValueError, match="must be orthogonal"):
        Rotation.from_davenport([ez, ezy], "e", [0, 0])
    with pytest.raises(ValueError, match="must be orthogonal"):
        Rotation.from_davenport([ez, ey, ezy], "e", [0, 0, 0])
    with pytest.raises(ValueError, match="order should be"):
        Rotation.from_davenport([ez], "xyz", [0])
    with pytest.raises(ValueError, match="Expected `angles`"):
        Rotation.from_davenport([ez, ey, ez], "e", [0, 1, 2, 3])


def test_as_davenport():
    rnd = np.random.RandomState(0)
    n = 100
    angles = np.empty((n, 3))
    angles[:, 0] = rnd.uniform(low=-np.pi, high=np.pi, size=(n,))
    angles_middle = rnd.uniform(low=0, high=np.pi, size=(n,))
    angles[:, 2] = rnd.uniform(low=-np.pi, high=np.pi, size=(n,))
    lambdas = rnd.uniform(low=0, high=np.pi, size=(20,))

    e1 = np.array([1, 0, 0])
    e2 = np.array([0, 1, 0])

    for lamb in lambdas:
        ax_lamb = [e1, e2, Rotation.from_rotvec(lamb * e2).apply(e1)]
        angles[:, 1] = angles_middle - lamb
        for order in ["extrinsic", "intrinsic"]:
            ax = ax_lamb if order == "intrinsic" else ax_lamb[::-1]
            rot = Rotation.from_davenport(ax, order, angles)
            angles_dav = rot.as_davenport(ax, order)
            assert_allclose(angles_dav, angles)


@pytest.mark.thread_unsafe
def test_as_davenport_degenerate():
    # Since we cannot check for angle equality, we check for rotation matrix
    # equality
    rnd = np.random.RandomState(0)
    n = 5
    angles = np.empty((n, 3))

    # symmetric sequences
    angles[:, 0] = rnd.uniform(low=-np.pi, high=np.pi, size=(n,))
    angles_middle = [rnd.choice([0, np.pi]) for i in range(n)]
    angles[:, 2] = rnd.uniform(low=-np.pi, high=np.pi, size=(n,))
    lambdas = rnd.uniform(low=0, high=np.pi, size=(5,))

    e1 = np.array([1, 0, 0])
    e2 = np.array([0, 1, 0])

    for lamb in lambdas:
        ax_lamb = [e1, e2, Rotation.from_rotvec(lamb * e2).apply(e1)]
        angles[:, 1] = angles_middle - lamb
        for order in ["extrinsic", "intrinsic"]:
            ax = ax_lamb if order == "intrinsic" else ax_lamb[::-1]
            rot = Rotation.from_davenport(ax, order, angles)
            with pytest.warns(UserWarning, match="Gimbal lock"):
                angles_dav = rot.as_davenport(ax, order)
            mat_expected = rot.as_matrix()
            mat_estimated = Rotation.from_davenport(ax, order, angles_dav).as_matrix()
            assert_array_almost_equal(mat_expected, mat_estimated)


def test_compare_from_davenport_from_euler():
    rnd = np.random.RandomState(0)
    n = 100
    angles = np.empty((n, 3))

    # symmetric sequences
    angles[:, 0] = rnd.uniform(low=-np.pi, high=np.pi, size=(n,))
    angles[:, 1] = rnd.uniform(low=0, high=np.pi, size=(n,))
    angles[:, 2] = rnd.uniform(low=-np.pi, high=np.pi, size=(n,))
    for order in ["extrinsic", "intrinsic"]:
        for seq_tuple in permutations("xyz"):
            seq = "".join([seq_tuple[0], seq_tuple[1], seq_tuple[0]])
            ax = [basis_vec(i) for i in seq]
            if order == "intrinsic":
                seq = seq.upper()
            eul = Rotation.from_euler(seq, angles)
            dav = Rotation.from_davenport(ax, order, angles)
            assert_allclose(
                eul.as_quat(canonical=True), dav.as_quat(canonical=True), rtol=1e-12
            )

    # asymmetric sequences
    angles[:, 1] -= np.pi / 2
    for order in ["extrinsic", "intrinsic"]:
        for seq_tuple in permutations("xyz"):
            seq = "".join(seq_tuple)
            ax = [basis_vec(i) for i in seq]
            if order == "intrinsic":
                seq = seq.upper()
            eul = Rotation.from_euler(seq, angles)
            dav = Rotation.from_davenport(ax, order, angles)
            assert_allclose(eul.as_quat(), dav.as_quat(), rtol=1e-12)


def test_compare_as_davenport_as_euler():
    rnd = np.random.RandomState(0)
    n = 100
    angles = np.empty((n, 3))

    # symmetric sequences
    angles[:, 0] = rnd.uniform(low=-np.pi, high=np.pi, size=(n,))
    angles[:, 1] = rnd.uniform(low=0, high=np.pi, size=(n,))
    angles[:, 2] = rnd.uniform(low=-np.pi, high=np.pi, size=(n,))
    for order in ["extrinsic", "intrinsic"]:
        for seq_tuple in permutations("xyz"):
            seq = "".join([seq_tuple[0], seq_tuple[1], seq_tuple[0]])
            ax = [basis_vec(i) for i in seq]
            if order == "intrinsic":
                seq = seq.upper()
            rot = Rotation.from_euler(seq, angles)
            eul = rot.as_euler(seq)
            dav = rot.as_davenport(ax, order)
            assert_allclose(eul, dav, rtol=1e-12)

    # asymmetric sequences
    angles[:, 1] -= np.pi / 2
    for order in ["extrinsic", "intrinsic"]:
        for seq_tuple in permutations("xyz"):
            seq = "".join(seq_tuple)
            ax = [basis_vec(i) for i in seq]
            if order == "intrinsic":
                seq = seq.upper()
            rot = Rotation.from_euler(seq, angles)
            eul = rot.as_euler(seq)
            dav = rot.as_davenport(ax, order)
            assert_allclose(eul, dav, rtol=1e-12)
