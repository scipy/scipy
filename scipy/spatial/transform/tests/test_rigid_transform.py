import pytest

import numpy as np
from numpy.testing import assert_allclose
from scipy.spatial.transform import Rotation, RigidTransform
from scipy.spatial.transform._rigid_transform import normalize_dual_quaternion
from scipy._lib._array_api import is_lazy_array, xp_vector_norm
import scipy._lib.array_api_extra as xpx


def test_repr(xp):
    actual = repr(RigidTransform.identity())
    expected = """\
RigidTransform.from_matrix(array([[1., 0., 0., 0.],
                                  [0., 1., 0., 0.],
                                  [0., 0., 1., 0.],
                                  [0., 0., 0., 1.]]))"""
    assert actual == expected

    actual = repr(RigidTransform.identity(2))
    expected = """\
RigidTransform.from_matrix(array([[[1., 0., 0., 0.],
                                   [0., 1., 0., 0.],
                                   [0., 0., 1., 0.],
                                   [0., 0., 0., 1.]],
                           
                                  [[1., 0., 0., 0.],
                                   [0., 1., 0., 0.],
                                   [0., 0., 1., 0.],
                                   [0., 0., 0., 1.]]]))"""
    assert actual == expected


def test_from_rotation(xp):
    atol = 1e-12

    # Test single rotation
    r = Rotation.from_matrix(xp.eye(3))
    tf = RigidTransform.from_rotation(r)
    assert_allclose(tf.as_matrix(), xp.eye(4), atol=atol)
    assert tf.single

    r = Rotation.from_euler("z", 90, degrees=True)
    r = Rotation.from_quat(xp.asarray(r.as_quat()))
    tf = RigidTransform.from_rotation(r)
    assert_allclose(tf.as_matrix()[:3, :3], r.as_matrix(), atol=atol)
    assert_allclose(tf.as_matrix()[:3, 3], xp.asarray([0, 0, 0]), atol=atol)
    assert_allclose(tf.as_matrix()[3, :], xp.asarray([0, 0, 0, 1]), atol=atol)
    assert tf.single

    # Test multiple rotations
    r = Rotation.from_euler("zyx", xp.asarray([[90, 0, 0], [0, 90, 0]]), degrees=True)
    tf = RigidTransform.from_rotation(r)
    assert_allclose(tf.as_matrix()[:, :3, :3], r.as_matrix(), atol=atol)
    assert_allclose(
        tf.as_matrix()[:, :3, 3], xp.asarray([[0, 0, 0], [0, 0, 0]]), atol=atol
    )
    assert_allclose(
        tf.as_matrix()[:, 3, :], xp.asarray([[0, 0, 0, 1], [0, 0, 0, 1]]), atol=atol
    )
    assert not tf.single


def test_from_rotation_jax_compile():
    pytest.importorskip("jax")
    import jax
    import jax.numpy as jp

    r = Rotation.from_matrix(jp.eye(3))
    from_rotation = jax.jit(RigidTransform.from_rotation)
    jax.block_until_ready(from_rotation(r))


def test_from_translation(xp):
    # Test single translation
    t = xp.asarray([1, 2, 3])
    tf = RigidTransform.from_translation(t)
    expected = xp.eye(4)
    expected = xpx.at(expected)[..., :3, 3].set(t)
    assert_allclose(tf.as_matrix(), expected)
    assert tf.single

    # Test multiple translations
    t = xp.asarray([[1, 2, 3], [4, 5, 6]])
    tf = RigidTransform.from_translation(t)
    for i in range(t.shape[0]):
        expected = xp.eye(4)
        expected = xpx.at(expected)[..., :3, 3].set(t[i, ...])
        assert_allclose(tf.as_matrix()[i, ...], expected)
    assert not tf.single


def test_from_translation_jax_compile():
    pytest.importorskip("jax")
    import jax
    import jax.numpy as jp

    t = jp.asarray([1, 2, 3])
    from_translation = jax.jit(RigidTransform.from_translation)
    jax.block_until_ready(from_translation(t))


def test_from_matrix(xp):
    atol = 1e-12

    # Test single transform matrix
    matrix = xp.eye(4)
    matrix = xpx.at(matrix)[..., :3, 3].set(xp.asarray([1, 2, 3]))
    tf = RigidTransform.from_matrix(matrix)
    assert_allclose(tf.as_matrix(), matrix, atol=atol)
    assert tf.single

    # Test multiple transform matrices
    # torch compat had issues with repeat, so we avoid using it here for testing. See
    # https://github.com/data-apis/array-api-compat/issues/292
    matrices = xp.zeros((2, 4, 4))
    matrices = xpx.at(matrices)[..., :4, :4].set(xp.eye(4))
    matrices = xpx.at(matrices)[..., :3, 3].set(xp.asarray([[1, 2, 3], [4, 5, 6]]))
    tf = RigidTransform.from_matrix(matrices)
    assert_allclose(tf.as_matrix(), matrices, atol=atol)
    assert not tf.single

    # Test non-1 determinant
    matrix = xp.eye(4)
    matrix = xpx.at(matrix)[..., :3, :3].set(xp.eye(3) * 2)
    tf = RigidTransform.from_matrix(matrix)
    assert_allclose(tf.as_matrix(), xp.eye(4), atol=atol)

    # Test non-orthogonal rotation matrix
    matrix = xp.asarray([[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    tf = RigidTransform.from_matrix(matrix)
    expected = xp.asarray(
        [
            [0.894427, 0.447214, 0, 0],
            [-0.447214, 0.894427, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ]
    )
    assert_allclose(tf.as_matrix(), expected, atol=1e-6)

    # Test invalid matrix
    invalid = xp.eye(4)
    invalid = xpx.at(invalid)[..., 3, 3].set(2)  # Invalid last row
    if is_lazy_array(invalid):
        tf = RigidTransform.from_matrix(invalid)
        assert xp.all(xp.isnan(tf.as_matrix()))
    else:
        with pytest.raises(ValueError):
            RigidTransform.from_matrix(invalid)


def test_from_matrix_jax_compile():
    pytest.importorskip("jax")
    import jax
    import jax.numpy as jp

    matrix = jp.eye(4)
    from_matrix = jax.jit(RigidTransform.from_matrix)
    jax.block_until_ready(from_matrix(matrix))


def test_from_components(xp):
    atol = 1e-12

    # Test single rotation and translation
    t = xp.asarray([1, 2, 3])
    r = Rotation.from_euler("zyx", xp.asarray([90, 0, 0]), degrees=True)
    tf = RigidTransform.from_components(t, r)

    expected = xp.zeros((4, 4))
    expected = xpx.at(expected)[..., :3, :3].set(r.as_matrix())
    expected = xpx.at(expected)[..., :3, 3].set(t)
    expected = xpx.at(expected)[..., 3, 3].set(1)
    assert_allclose(tf.as_matrix(), expected, atol=atol)
    assert tf.single

    # Test single rotation and multiple translations
    t = xp.asarray([[1, 2, 3], [4, 5, 6]])
    r = Rotation.from_euler("z", 90, degrees=True)
    r = Rotation.from_quat(xp.asarray(r.as_quat()))
    tf = RigidTransform.from_components(t, r)
    assert not tf.single

    for i in range(t.shape[0]):
        expected = xp.zeros((4, 4))
        expected = xpx.at(expected)[..., :3, :3].set(r.as_matrix())
        expected = xpx.at(expected)[..., :3, 3].set(t[i, ...])
        expected = xpx.at(expected)[..., 3, 3].set(1)
        assert_allclose(tf.as_matrix()[i, ...], expected, atol=atol)

    # Test multiple rotations and translations
    t = xp.asarray([[1, 2, 3], [4, 5, 6]])
    r = Rotation.from_euler("zyx", xp.asarray([[90, 0, 0], [0, 90, 0]]), degrees=True)
    tf = RigidTransform.from_components(t, r)
    assert not tf.single

    for i in range(t.shape[0]):
        expected = xp.zeros((4, 4))
        expected = xpx.at(expected)[..., :3, :3].set(r.as_matrix()[i, ...])
        expected = xpx.at(expected)[..., :3, 3].set(t[i, ...])
        expected = xpx.at(expected)[..., 3, 3].set(1)
        assert_allclose(tf.as_matrix()[i, ...], expected, atol=atol)


def test_from_components_jax_compile():
    pytest.importorskip("jax")
    import jax
    import jax.numpy as jp

    t = jp.asarray([1, 2, 3])
    r = Rotation.from_euler("zyx", jp.asarray([90, 0, 0]), degrees=True)
    from_components = jax.jit(RigidTransform.from_components)
    jax.block_until_ready(from_components(t, r))


def test_as_components(xp):
    atol = 1e-12
    n = 10
    rng = np.random.default_rng(123)
    t = xp.asarray(rng.normal(size=(n, 3)))
    r = Rotation.random(n, rng=rng)
    r = Rotation.from_quat(xp.asarray(r.as_quat()))
    tf = RigidTransform.from_components(t, r)
    new_t, new_r = tf.as_components()
    assert all(new_r.approx_equal(r, atol=atol))
    assert_allclose(new_t, t, atol=atol)


def test_as_components_jax_compile():
    pytest.importorskip("jax")
    import jax
    import jax.numpy as jp

    t = jp.asarray([1, 2, 3])
    r = Rotation.from_euler("zyx", jp.asarray([90, 0, 0]), degrees=True)
    tf = RigidTransform.from_components(t, r)
    as_components = jax.jit(RigidTransform.as_components)
    jax.block_until_ready(as_components(tf))


def test_from_exp_coords(xp):
    # example from 3.3 of
    # https://hades.mech.northwestern.edu/images/2/25/MR-v2.pdf
    angle1 = xp.asarray(30.0 / 180 * np.pi)  # deg2rad is not implemented in Array API
    tf1 = RigidTransform.from_matrix(
        xp.asarray(
            [
                [xp.cos(angle1), -xp.sin(angle1), 0.0, 1.0],
                [xp.sin(angle1), xp.cos(angle1), 0.0, 2.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )
    )
    angle2 = xp.asarray(60.0 / 180 * np.pi)
    tf2 = RigidTransform.from_matrix(
        xp.asarray(
            [
                [xp.cos(angle2), -xp.sin(angle2), 0.0, 2.0],
                [xp.sin(angle2), xp.cos(angle2), 0.0, 1.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )
    )
    expected = tf2 * tf1.inv()
    angle = xp.asarray(30.0 / 180 * np.pi)
    actual = RigidTransform.from_exp_coords(
        angle * xp.asarray([0.0, 0.0, 1.0, 3.37, -3.37, 0.0])
    )
    assert_allclose(actual.as_matrix(), expected.as_matrix(), atol=1e-2)

    # test cases generated by comparison to pytransform3d
    exp_coords = xp.asarray(
        [
            [-2.01041204, -0.52983629, 0.65773501, 0.10386614, 0.05855009, 0.54959179],
            [
                -0.22537438,
                -0.24132627,
                -2.4747121,
                -0.09158594,
                1.88075832,
                -0.03197204,
            ],
        ]
    )
    expected_matrix = [
        [
            [0.76406621, 0.10504613, -0.63652819, -0.10209961],
            [0.59956454, -0.47987325, 0.64050295, 0.40158789],
            [-0.2381705, -0.87102639, -0.42963687, 0.19637636],
            [0.0, 0.0, 0.0, 1.0],
        ],
        [
            [-0.78446989, 0.61157488, 0.10287448, 1.33330055],
            [-0.58017785, -0.78232107, 0.22664378, 0.52660831],
            [0.21909052, 0.11810973, 0.96852952, -0.02968529],
            [0.0, 0.0, 0.0, 1.0],
        ],
    ]
    assert_allclose(
        RigidTransform.from_exp_coords(exp_coords).as_matrix(),
        expected_matrix,
        atol=1e-8,
    )

    # identity
    assert_allclose(
        RigidTransform.from_exp_coords(xp.zeros(6)).as_matrix(), xp.eye(4), atol=1e-12
    )

    # only translation
    expected_matrix = xp.asarray(
        [
            [
                [1.0, 0.0, 0.0, 3.0],
                [0.0, 1.0, 0.0, -5.4],
                [0.0, 0.0, 1.0, 100.2],
                [0.0, 0.0, 0.0, 1.0],
            ],
            [
                [1.0, 0.0, 0.0, -3.0],
                [0.0, 1.0, 0.0, 13.3],
                [0.0, 0.0, 1.0, 1.3],
                [0.0, 0.0, 0.0, 1.0],
            ],
        ]
    )
    actual = RigidTransform.from_exp_coords(
        xp.asarray(
            [
                [0.0, 0.0, 0.0, 3.0, -5.4, 100.2],
                [0.0, 0.0, 0.0, -3.0, 13.3, 1.3],
            ]
        )
    )
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)

    # only rotation
    rot = Rotation.from_euler(
        "zyx", xp.asarray([[34, -12, 0.5], [-102, -55, 30]]), degrees=True
    )
    rotvec = rot.as_rotvec()
    expected_matrix = xp.zeros((2, 4, 4))
    expected_matrix = xpx.at(expected_matrix)[..., :3, :3].set(rot.as_matrix())
    expected_matrix = xpx.at(expected_matrix)[..., 3, 3].set(1)
    exp_coords = xp.concat((rotvec, xp.zeros((2, 3))), axis=-1)
    actual = RigidTransform.from_exp_coords(exp_coords)
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)


def test_from_exp_coords_jax_compile():
    pytest.importorskip("jax")
    import jax
    import jax.numpy as jp

    exp_coords = jp.asarray([0.0, 0.0, 0.0, 3.0, -5.4, 100.2])
    from_exp_coords = jax.jit(RigidTransform.from_exp_coords)
    jax.block_until_ready(from_exp_coords(exp_coords))


def test_as_exp_coords(xp):
    # identity
    expected = xp.zeros(6)
    actual = RigidTransform.from_exp_coords(expected).as_exp_coords()
    assert_allclose(actual, expected, atol=1e-12)

    rng = np.random.default_rng(10)

    # pure rotation
    rot_vec = xp.asarray(rng.normal(scale=0.1, size=(1000, 3)))
    tf = RigidTransform.from_rotation(Rotation.from_rotvec(rot_vec))
    exp_coords = tf.as_exp_coords()
    assert_allclose(exp_coords[:, :3], rot_vec, rtol=1e-13)
    assert_allclose(exp_coords[:, 3:], 0.0, atol=1e-16)

    # pure translation
    translation = xp.asarray(rng.normal(scale=100.0, size=(1000, 3)))
    tf = RigidTransform.from_translation(translation)
    exp_coords = tf.as_exp_coords()
    assert_allclose(exp_coords[:, :3], 0.0, atol=1e-16)
    assert_allclose(exp_coords[:, 3:], translation, rtol=1e-15)


def test_as_exp_coords_jax_compile():
    pytest.importorskip("jax")
    import jax
    import jax.numpy as jp

    tf = RigidTransform.from_translation(jp.asarray([1, 2, 3]))
    as_exp_coords = jax.jit(RigidTransform.as_exp_coords)
    jax.block_until_ready(as_exp_coords(tf))


def test_from_dual_quat(xp):
    # identity
    assert_allclose(
        RigidTransform.from_dual_quat(
            xp.asarray([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0])
        ).as_matrix(),
        xp.eye(4),
        atol=1e-12,
    )
    assert_allclose(
        RigidTransform.from_dual_quat(
            xp.asarray([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]), scalar_first=True
        ).as_matrix(),
        xp.eye(4),
        atol=1e-12,
    )

    # only translation
    actual = RigidTransform.from_dual_quat(
        xp.asarray([0, 0, 0, 1, 0.25, 0.15, -0.7, 0])
    )
    expected_matrix = xp.asarray(
        [[1, 0, 0, 0.5], [0, 1, 0, 0.3], [0, 0, 1, -1.4], [0, 0, 0, 1]]
    )
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)
    actual = RigidTransform.from_dual_quat(
        xp.asarray([1, 0, 0, 0, 0, 0.25, 0.15, -0.7]), scalar_first=True
    )
    expected_matrix = xp.asarray(
        [[1, 0, 0, 0.5], [0, 1, 0, 0.3], [0, 0, 1, -1.4], [0, 0, 0, 1]]
    )
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)

    # only rotation
    actual_rot = Rotation.from_euler("xyz", xp.asarray([65, -13, 90]), degrees=True)
    actual = RigidTransform.from_dual_quat(
        xp.concat((actual_rot.as_quat(), xp.zeros(4)), axis=-1)
    )
    expected_matrix = xp.eye(4)
    expected_matrix = xpx.at(expected_matrix)[..., :3, :3].set(actual_rot.as_matrix())
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)

    actual = RigidTransform.from_dual_quat(
        xp.concat((actual_rot.as_quat(scalar_first=True), xp.zeros(4)), axis=-1),
        scalar_first=True,
    )
    expected_matrix = xp.eye(4)
    expected_matrix = xpx.at(expected_matrix)[..., :3, :3].set(actual_rot.as_matrix())
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)

    # rotation and translation
    actual = RigidTransform.from_dual_quat(
        xp.asarray(
            [
                [
                    0.0617101,
                    -0.06483886,
                    0.31432811,
                    0.94508498,
                    0.04985168,
                    -0.26119618,
                    0.1691491,
                    -0.07743254,
                ],
                [
                    0.19507259,
                    0.49404931,
                    -0.06091285,
                    0.8450749,
                    0.65049656,
                    -0.30782513,
                    0.16566752,
                    0.04174544,
                ],
            ]
        )
    )
    expected_matrix = xp.asarray(
        [
            [
                [0.79398752, -0.60213598, -0.08376202, 0.24605262],
                [0.58613113, 0.79477941, -0.15740392, -0.4932833],
                [0.16135089, 0.07588122, 0.98397557, 0.34262676],
                [0.0, 0.0, 0.0, 1.0],
            ],
            [
                [0.50440981, 0.2957028, 0.81125249, 1.20934468],
                [0.08979911, 0.91647262, -0.3898898, -0.70540077],
                [-0.8587822, 0.26951399, 0.43572393, -0.47776265],
                [0.0, 0.0, 0.0, 1.0],
            ],
        ]
    )
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)

    actual = RigidTransform.from_dual_quat(
        xp.asarray(
            [
                [
                    0.94508498,
                    0.0617101,
                    -0.06483886,
                    0.31432811,
                    -0.07743254,
                    0.04985168,
                    -0.26119618,
                    0.1691491,
                ],
                [
                    0.8450749,
                    0.19507259,
                    0.49404931,
                    -0.06091285,
                    0.04174544,
                    0.65049656,
                    -0.30782513,
                    0.16566752,
                ],
            ]
        ),
        scalar_first=True,
    )
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)

    # unnormalized dual quaternions

    # invalid real quaternion with norm 0
    actual = RigidTransform.from_dual_quat(xp.zeros(8))
    assert_allclose(actual.as_matrix(), xp.eye(4), atol=1e-12)

    # real quaternion with norm != 1
    unnormalized_dual_quat = xp.asarray(
        [
            -0.2547655,
            1.23506123,
            0.20230088,
            0.24247194,  # norm 1.3
            0.38559628,
            0.08184063,
            0.1755943,
            -0.1582222,
        ]  # orthogonal
    )
    assert_allclose(xp_vector_norm(unnormalized_dual_quat[:4]), 1.3, atol=1e-12)
    assert_allclose(
        xp.vecdot(unnormalized_dual_quat[:4], unnormalized_dual_quat[4:]),
        0.0,
        atol=1e-8,
    )
    dual_quat = RigidTransform.from_dual_quat(unnormalized_dual_quat).as_dual_quat()
    assert_allclose(xp_vector_norm(dual_quat[:4]), 1.0, atol=1e-12)
    assert_allclose(xp.vecdot(dual_quat[:4], dual_quat[4:]), 0.0, atol=1e-12)

    # real and dual quaternion are not orthogonal
    unnormalized_dual_quat = xp.asarray(
        [
            0.20824458,
            0.75098079,
            0.54542913,
            -0.30849493,  # unit norm
            -0.16051025,
            0.10742978,
            0.21277201,
            0.20596935,
        ]  # not orthogonal
    )
    assert_allclose(xp_vector_norm(unnormalized_dual_quat[:4]), 1.0, atol=1e-12)
    assert xp.vecdot(unnormalized_dual_quat[:4], unnormalized_dual_quat[4:]) != 0.0

    dual_quat = RigidTransform.from_dual_quat(unnormalized_dual_quat).as_dual_quat()
    assert_allclose(xp_vector_norm(dual_quat[:4]), 1.0, atol=1e-12)
    assert_allclose(xp.vecdot(dual_quat[:4], dual_quat[4:]), 0.0, atol=1e-12)

    # invalid real quaternion with norm 0, non-orthogonal dual quaternion
    unnormalized_dual_quat = xp.asarray(
        [0.0, 0.0, 0.0, 0.0, -0.16051025, 0.10742978, 0.21277201, 0.20596935]
    )
    assert (
        xp.vecdot(xp.asarray([0.0, 0.0, 0.0, 1.0]), unnormalized_dual_quat[4:]) != 0.0
    )
    dual_quat = RigidTransform.from_dual_quat(unnormalized_dual_quat).as_dual_quat()
    assert_allclose(dual_quat[:4], xp.asarray([0, 0, 0, 1]), atol=1e-12)
    assert_allclose(xp.vecdot(dual_quat[:4], dual_quat[4:]), 0.0, atol=1e-12)

    # compensation for precision loss in real quaternion
    rng = np.random.default_rng(1000)
    t = xp.asarray(rng.normal(size=(3,)))
    r = Rotation.random(10, rng=rng)
    r = Rotation.from_quat(xp.asarray(r.as_quat()))
    random_dual_quats = RigidTransform.from_components(t, r).as_dual_quat()

    # ensure that random quaternions are not normalized
    random_dual_quats = xpx.at(random_dual_quats)[:, :4].set(
        random_dual_quats[:, :4] + 0.01
    )
    q_norm = xp_vector_norm(random_dual_quats[:, :4], axis=1)
    assert not xp.any(xp.abs(q_norm - 1.0) < 0.0001)
    dual_quat_norm = RigidTransform.from_dual_quat(random_dual_quats).as_dual_quat()
    assert_allclose(xp_vector_norm(dual_quat_norm[:, :4], axis=1), 1.0, atol=1e-12)

    # compensation for precision loss in dual quaternion, results in violation
    # of orthogonality constraint
    t = xp.asarray(rng.normal(size=(10, 3)))
    r = Rotation.random(10, rng=rng)
    r = Rotation.from_quat(xp.asarray(r.as_quat()))
    random_dual_quats = RigidTransform.from_components(t, r).as_dual_quat()

    # ensure that random quaternions are not normalized
    random_dual_quats = xpx.at(random_dual_quats)[:, 4:].set(
        random_dual_quats[:, 4:] + 0.01
    )
    q_norm = xp.vecdot(random_dual_quats[:, :4], random_dual_quats[:, 4:])
    assert not xp.any(xp.abs(q_norm) < 0.0001)
    dual_quat_norm = RigidTransform.from_dual_quat(random_dual_quats).as_dual_quat()
    assert_allclose(
        xp.vecdot(dual_quat_norm[:, :4], dual_quat_norm[:, 4:]),
        0.0,
        atol=1e-12,
    )
    assert_allclose(random_dual_quats[:, :4], dual_quat_norm[:, :4], atol=1e-12)


def test_from_dual_quat_jax_compile():
    pytest.importorskip("jax")
    import jax
    import jax.numpy as jp

    q = jp.asarray([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0])
    from_dual_quat = jax.jit(RigidTransform.from_dual_quat)
    jax.block_until_ready(from_dual_quat(q))


def test_as_dual_quat(xp):
    # identity
    expected = xp.asarray([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0])
    actual = RigidTransform.identity().as_dual_quat()
    actual = xp.asarray(actual)
    assert_allclose(actual, expected, atol=1e-12)

    expected = xp.asarray([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    actual = RigidTransform.identity().as_dual_quat(scalar_first=True)
    actual = xp.asarray(actual)
    assert_allclose(actual, expected, atol=1e-12)

    rng = np.random.default_rng(10)

    # only rotation
    for _ in range(10):
        real_part = xp.asarray(Rotation.random(rng=rng).as_quat())
        dual_part = xp.zeros(4)
        expected = xp.concat((real_part, dual_part), axis=-1)
        actual = xp.asarray(RigidTransform.from_dual_quat(expected).as_dual_quat())
        # because of double cover:
        if xp.sign(expected[0]) != xp.sign(actual[0]):
            actual = -actual
        assert_allclose(actual, expected, atol=1e-12)

    # only translation
    for _ in range(10):
        tf = xp.asarray(rng.normal(size=3))
        expected = xp.asarray([0.0, 0, 0, 1, *(0.5 * tf), 0])
        actual = xp.asarray(RigidTransform.from_dual_quat(expected).as_dual_quat())
        # because of double cover:
        if xp.sign(expected[0]) != xp.sign(actual[0]):
            actual = -actual
        assert_allclose(actual, expected, atol=1e-12)

    # rotation and translation
    for _ in range(10):
        t = xp.asarray(rng.normal(size=3))
        r = Rotation.from_quat(xp.asarray(Rotation.random(rng=rng).as_quat()))
        expected = RigidTransform.from_components(t, r).as_dual_quat()
        actual = xp.asarray(RigidTransform.from_dual_quat(expected).as_dual_quat())
        # because of double cover:
        if xp.sign(expected[0]) != xp.sign(actual[0]):
            actual = -actual
        assert_allclose(actual, expected, atol=1e-12)


def test_as_dual_quat_jax_compile():
    pytest.importorskip("jax")
    import jax
    import jax.numpy as jp

    tf = RigidTransform.from_matrix(jp.eye(4))
    as_dual_quat = jax.jit(RigidTransform.as_dual_quat)
    jax.block_until_ready(as_dual_quat(tf))


def test_from_as_internal_consistency(xp):
    atol = 1e-12
    n = 1000
    rng = np.random.default_rng(10)
    t = xp.asarray(rng.normal(size=(n, 3)))
    r = Rotation.random(n, rng=rng)
    r = Rotation.from_quat(xp.asarray(r.as_quat()))
    tf0 = RigidTransform.from_components(t, r)

    tf1 = RigidTransform.from_components(*tf0.as_components())
    assert_allclose(tf0.as_matrix(), tf1.as_matrix(), atol=atol)

    tf1 = RigidTransform.from_components(tf0.translation, tf0.rotation)
    assert_allclose(tf0.as_matrix(), tf1.as_matrix(), atol=atol)

    tf1 = RigidTransform.from_exp_coords(tf0.as_exp_coords())
    assert_allclose(tf0.as_matrix(), tf1.as_matrix(), atol=atol)

    tf1 = RigidTransform.from_matrix(tf0.as_matrix())
    assert_allclose(tf0.as_matrix(), tf1.as_matrix(), atol=atol)

    tf1 = RigidTransform.from_dual_quat(tf0.as_dual_quat())
    assert_allclose(tf0.as_matrix(), tf1.as_matrix(), atol=atol)

    # exp_coords small rotation
    t = xp.asarray(rng.normal(scale=1000.0, size=(1000, 3)))
    rot_vec = xp.asarray(rng.normal(scale=1e-10, size=(1000, 3)))
    tf0 = RigidTransform.from_components(t, Rotation.from_rotvec(rot_vec))
    tf1 = RigidTransform.from_exp_coords(tf0.as_exp_coords())
    assert_allclose(tf0.as_matrix(), tf1.as_matrix(), atol=atol)


def test_identity():
    # We do not use xp here because identity always returns numpy arrays
    atol = 1e-12

    # Test single identity
    tf = RigidTransform.identity()
    assert_allclose(tf.as_matrix(), np.eye(4), atol=atol)

    # Test multiple identities
    tf = RigidTransform.identity(5)
    assert_allclose(tf.as_matrix(), np.array([np.eye(4)] * 5), atol=atol)


def test_apply(xp):
    atol = 1e-12

    ## Single transform
    r = Rotation.from_euler("z", 90, degrees=True)
    r = Rotation.from_quat(xp.asarray(r.as_quat()))
    t = xp.asarray([2.0, 3, 4])
    tf = RigidTransform.from_components(t, r)

    # Single vector, single transform
    vec = xp.asarray([1.0, 0, 0])
    expected = t + r.apply(vec)
    res = tf.apply(vec)
    assert_allclose(res, expected, atol=atol)

    # Multiple vectors, single transform
    vecs = xp.asarray([[1.0, 0, 0], [0, 1, 0]])
    expected = t + r.apply(vecs)
    assert_allclose(tf.apply(vecs), expected, atol=atol)

    ## Multiple transforms
    r = Rotation.from_euler("z", xp.asarray([90, 0]), degrees=True)
    t = xp.asarray([[2.0, 3, 4], [5, 6, 7]])
    tf = RigidTransform.from_components(t, r)

    # Single vector, multiple transforms
    vec = xp.asarray([1.0, 0, 0])
    expected = t + r.apply(vec)
    assert_allclose(tf.apply(vec), expected, atol=atol)

    # Multiple vectors, multiple transforms
    vecs = xp.asarray([[1.0, 0, 0], [0, 1, 0]])
    expected = t + r.apply(vecs)
    assert_allclose(tf.apply(vecs), expected, atol=atol)


def test_apply_jax_compile():
    pytest.importorskip("jax")
    import jax
    import jax.numpy as jp

    tf = RigidTransform.from_matrix(jp.eye(4))
    apply = jax.jit(RigidTransform.apply)
    jax.block_until_ready(apply(tf, jp.array([1, 0, 0])))


def test_inverse_apply(xp):
    atol = 1e-12

    # Test applying inverse transform
    t = xp.asarray([1, 2, 3])
    r = Rotation.from_euler("z", 90, degrees=True)
    r = Rotation.from_quat(xp.asarray(r.as_quat()))
    tf = RigidTransform.from_components(t, r)

    # Test single vector
    vec = xp.asarray([1, 0, 0])
    expected = tf.inv().apply(vec)
    assert_allclose(tf.apply(vec, inverse=True), expected, atol=atol)

    # Test multiple vectors
    vecs = xp.asarray([[1, 0, 0], [0, 1, 0]])
    expected = tf.inv().apply(vecs)
    assert_allclose(tf.apply(vecs, inverse=True), expected, atol=atol)


def test_rotation_alone(xp):
    atol = 1e-12

    r = Rotation.from_euler("z", 90, degrees=True)
    r = Rotation.from_quat(xp.asarray(r.as_quat()))
    tf = RigidTransform.from_rotation(r)
    vec = xp.asarray([1, 0, 0])
    expected = r.apply(vec)
    assert_allclose(tf.apply(vec), expected, atol=atol)


def test_translation_alone(xp):
    atol = 1e-12
    t = xp.asarray([1.0, 2, 3])
    tf = RigidTransform.from_translation(t)
    vec = xp.asarray([5.0, 6, 7])
    expected = t + vec
    assert_allclose(tf.apply(vec), expected, atol=atol)


def test_composition(xp):
    atol = 1e-12

    # Test composing single transforms
    t1 = xp.asarray([1.0, 0, 0])
    r1 = Rotation.from_euler("z", 90, degrees=True)
    r1 = Rotation.from_quat(xp.asarray(r1.as_quat()))
    tf1 = RigidTransform.from_components(t1, r1)

    t2 = xp.asarray([0.0, 1, 0])
    r2 = Rotation.from_euler("x", 90, degrees=True)
    r2 = Rotation.from_quat(xp.asarray(r2.as_quat()))
    tf2 = RigidTransform.from_components(t2, r2)

    composed = tf2 * tf1
    vec = xp.asarray([1.0, 0, 0])
    expected = tf2.apply(tf1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)
    assert composed.single

    expected = t2 + r2.apply(t1 + r1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)

    # Multiple transforms with single transform
    t2 = xp.asarray([[1.0, 2, 3], [4, 5, 6]])
    tf2 = RigidTransform.from_components(t2, r2)

    composed = tf2 * tf1
    expected = tf2.apply(tf1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)
    assert not composed.single

    expected = t2 + r2.apply(t1 + r1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)

    # Multiple transforms with multiple transforms
    t1 = xp.asarray([[1.0, 0, 0], [0, -1, 1]])
    tf1 = RigidTransform.from_components(t1, r1)

    composed = tf2 * tf1
    expected = tf2.apply(tf1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)
    assert not composed.single

    expected = t2 + r2.apply(t1 + r1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)


def test_pow(xp):
    atol = 1e-12
    num = 10
    rng = np.random.default_rng(100)
    t = xp.asarray(rng.normal(size=(num, 3)))
    r = Rotation.random(num, rng=rng)
    r = Rotation.from_quat(xp.asarray(r.as_quat()))
    p = RigidTransform.from_components(t, r)
    p_inv = p.inv()

    # Test the short-cuts and other integers
    for n in [-5, -2, -1, 0, 1, 2, 5]:
        q = p**n
        r = RigidTransform.identity(num)
        r = RigidTransform.from_matrix(xp.asarray(r.as_matrix()))
        for _ in range(abs(n)):
            if n > 0:
                r = r * p
            else:
                r = r * p_inv
        assert_allclose(q.as_matrix(), r.as_matrix(), atol=atol)

        # Test shape preservation
        r = RigidTransform.from_rotation(Rotation.from_quat(xp.asarray([0, 0, 0, 1])))
        assert (r**n).as_matrix().shape == (4, 4)
        r = RigidTransform.from_rotation(Rotation.from_quat(xp.asarray([[0, 0, 0, 1]])))
        assert (r**n).as_matrix().shape == (1, 4, 4)

    # Test fractional powers
    q = p**0.5
    assert_allclose((q * q).as_matrix(), p.as_matrix(), atol=atol)
    q = p**-0.5
    assert_allclose((q * q).as_matrix(), p.inv().as_matrix(), atol=atol)
    q = p**1.5
    assert_allclose((q * q).as_matrix(), (p**3).as_matrix(), atol=atol)
    q = p**-1.5
    assert_allclose((q * q).as_matrix(), (p**-3).as_matrix(), atol=atol)

    # pow function
    identity = RigidTransform.from_matrix(xp.eye(4))
    tf = pow(identity, 2)
    assert_allclose(tf.as_matrix(), xp.eye(4), atol=atol)


def test_pow_jax_compile():
    pytest.importorskip("jax")
    import jax
    import jax.numpy as jp

    identity = RigidTransform.from_matrix(jp.eye(4))
    pow = jax.jit(RigidTransform.__pow__)
    jax.block_until_ready(pow(identity, 2))


def test_pow_equivalence_with_rotation(xp):
    atol = 1e-12
    num = 10
    rng = np.random.default_rng(100)
    r = Rotation.random(num, rng=rng)
    r = Rotation.from_quat(xp.asarray(r.as_quat()))
    p = RigidTransform.from_rotation(r)
    for n in [-5, -2, -1.5, -1, -0.5, 0.0, 0.5, 1, 1.5, 2, 5]:
        assert_allclose((p**n).rotation.as_matrix(), (r**n).as_matrix(), atol=atol)


def test_inverse(xp):
    atol = 1e-12

    # Test inverse transform
    r = Rotation.from_euler("z", 90, degrees=True)
    r = Rotation.from_quat(xp.asarray(r.as_quat()))
    t = xp.asarray([1, 2, 3])
    tf = RigidTransform.from_components(t, r)

    # Test that tf * tf.inv() equals identity
    tf_inv = tf.inv()
    composed = tf * tf_inv
    assert_allclose(composed.as_matrix(), xp.eye(4), atol=atol)

    n = 10
    rng = np.random.default_rng(1000)
    t = rng.normal(size=(n, 3))
    r = Rotation.random(n, rng=rng)
    tf = RigidTransform.from_components(t, r)
    tf_inv = tf.inv()
    composed = tf * tf_inv
    expected = xp.zeros((n, 4, 4))
    expected = xpx.at(expected)[...].set(xp.eye(4))
    assert_allclose(composed.as_matrix(), expected, atol=atol)

    # Test multiple transforms
    r = Rotation.from_euler("zyx", xp.asarray([[90, 0, 0], [0, 90, 0]]), degrees=True)
    t = xp.asarray([[1, 2, 3], [4, 5, 6]])
    tf = RigidTransform.from_components(t, r)
    tf_inv = tf.inv()
    composed = tf * tf_inv
    expected = xp.zeros((2, 4, 4))
    expected = xpx.at(expected)[...].set(xp.eye(4))
    assert_allclose(composed.as_matrix(), expected, atol=atol)


def test_properties(xp):
    atol = 1e-12

    # Test rotation and translation properties for single transform
    r = Rotation.from_euler("z", 90, degrees=True)
    r = Rotation.from_quat(xp.asarray(r.as_quat()))
    t = xp.asarray([1, 2, 3])
    tf = RigidTransform.from_components(t, r)

    assert_allclose(tf.rotation.as_matrix(), r.as_matrix(), atol=atol)
    assert tf.rotation.approx_equal(r)
    assert_allclose(tf.translation, t, atol=atol)

    # Test rotation and translation properties for multiple transforms
    r = Rotation.from_euler("zyx", xp.asarray([[90, 0, 0], [0, 90, 0]]), degrees=True)
    t = xp.asarray([[1, 2, 3], [4, 5, 6]])
    tf = RigidTransform.from_components(t, r)

    assert_allclose(tf.rotation.as_matrix(), r.as_matrix(), atol=atol)
    assert all(tf.rotation.approx_equal(r))
    assert_allclose(tf.translation, t, atol=atol)


def test_indexing():
    atol = 1e-12

    # Test indexing for multiple transforms
    r = Rotation.from_euler("zyx", [[90, 0, 0], [0, 90, 0]], degrees=True)
    t = np.array([[1, 2, 3], [4, 5, 6]])
    tf = RigidTransform.from_components(t, r)

    # Test single index
    assert_allclose(tf[0].as_matrix()[:3, :3], r[0].as_matrix(), atol=atol)
    assert_allclose(tf[0].as_matrix()[:3, 3], t[0], atol=atol)

    # Test slice
    tf_slice = tf[0:2]
    assert_allclose(tf_slice.as_matrix()[:, :3, :3], r[0:2].as_matrix(), atol=atol)
    assert_allclose(tf_slice.as_matrix()[:, :3, 3], t[0:2], atol=atol)

    # Test boolean indexing
    tf_masked = tf[[True, True]]
    assert_allclose(tf_masked.as_matrix()[:, :3, :3], r.as_matrix(), atol=atol)
    assert_allclose(tf_masked.as_matrix()[:, :3, 3], t, atol=atol)

    tf_masked = tf[[False, True]]
    assert_allclose(
        tf_masked.as_matrix()[:, :3, :3], r[[False, True]].as_matrix(), atol=atol
    )
    assert_allclose(tf_masked.as_matrix()[:, :3, 3], t[[False, True]], atol=atol)

    tf_masked = tf[[False, False]]
    assert len(tf_masked) == 0


def test_concatenate():
    atol = 1e-12

    # Test concatenation of transforms
    t1 = np.array([1, 0, 0])
    r1 = Rotation.from_euler("z", 90, degrees=True)
    tf1 = RigidTransform.from_components(t1, r1)

    t2 = np.array([0, 1, 0])
    r2 = Rotation.from_euler("x", 90, degrees=True)
    tf2 = RigidTransform.from_components(t2, r2)

    # Concatenate single transforms
    concatenated1 = RigidTransform.concatenate([tf1, tf2])
    assert_allclose(concatenated1[0].as_matrix(), tf1.as_matrix(), atol=atol)
    assert_allclose(concatenated1[1].as_matrix(), tf2.as_matrix(), atol=atol)

    # Concatenate multiple transforms
    concatenated2 = RigidTransform.concatenate([tf1, concatenated1])
    assert_allclose(concatenated2[0].as_matrix(), tf1.as_matrix(), atol=atol)
    assert_allclose(concatenated2[1].as_matrix(), tf1.as_matrix(), atol=atol)
    assert_allclose(concatenated2[2].as_matrix(), tf2.as_matrix(), atol=atol)


def test_input_validation():
    # Test invalid matrix shapes
    inputs = [np.eye(3), np.zeros((4, 3)), [], np.zeros((1, 1, 4, 4))]
    for input in inputs:
        with pytest.raises(ValueError, match="Expected `matrix` to have shape"):
            RigidTransform.from_matrix(input)

    # Test invalid last row
    with pytest.raises(ValueError, match="last row of transformation matrix 0"):
        matrix = np.eye(4)
        matrix[3, :] = [1, 0, 0, 1]
        RigidTransform.from_matrix(matrix)

    # Test invalid last row for multiple transforms
    with pytest.raises(ValueError, match="last row of transformation matrix 1"):
        matrix = np.array([np.eye(4)] * 2)
        matrix[1, 3, :] = [1, 0, 0, 1]
        RigidTransform.from_matrix(matrix)

    # Test left handed rotation matrix
    with pytest.raises(ValueError, match="Non-positive determinant"):
        matrix = np.eye(4)
        matrix[0, 0] = -1
        RigidTransform(matrix, normalize=True)

    # Test non-Rotation input
    with pytest.raises(
        ValueError, match="Expected `rotation` to be a `Rotation` instance"
    ):
        RigidTransform.from_rotation(np.eye(3))


def test_translation_validation():
    # Test invalid translation shapes
    with pytest.raises(ValueError, match="Expected `translation` to have shape"):
        RigidTransform.from_translation([1, 2])

    with pytest.raises(ValueError, match="Expected `translation` to have shape"):
        RigidTransform.from_translation(np.zeros((2, 2)))

    with pytest.raises(ValueError, match="Expected `translation` to have shape"):
        RigidTransform.from_translation(np.zeros((1, 1, 3)))


def test_vector_validation():
    tf = RigidTransform.identity(2)

    # Test invalid vector shapes
    with pytest.raises(ValueError, match="Expected vector to have shape"):
        tf.apply([1, 2])

    with pytest.raises(ValueError, match="Expected vector to have shape"):
        tf.apply(np.zeros((2, 2)))

    with pytest.raises(ValueError, match="Expected vector to have shape"):
        tf.apply(np.zeros((1, 1, 3)))


def test_indexing_validation():
    tf = RigidTransform.identity()

    # Test indexing on single transform
    with pytest.raises(TypeError, match="Single transform is not subscriptable"):
        tf[0]

    with pytest.raises(TypeError, match="Single transform is not subscriptable"):
        tf[0:1]

    # Test length on single transform
    with pytest.raises(TypeError, match="Single transform has no len"):
        len(tf)


def test_composition_validation():
    tf2 = RigidTransform.from_translation([[1, 2, 3], [4, 5, 6]])
    tf3 = RigidTransform.from_translation([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    # Test incompatible shapes
    with pytest.raises(ValueError, match="Expected equal number of transforms"):
        tf2 * tf3


def test_concatenate_validation():
    tf = RigidTransform.identity()

    # Test invalid inputs
    with pytest.raises(TypeError, match="input must contain RigidTransform objects"):
        RigidTransform.concatenate([tf, np.eye(4)])


def test_setitem_validation():
    tf = RigidTransform.from_translation([[1, 2, 3], [4, 5, 6]])
    single = RigidTransform.identity()

    # Test setting item on single transform
    with pytest.raises(TypeError, match="Single transform is not subscriptable"):
        single[0] = tf

    # Test invalid value type
    with pytest.raises(TypeError, match="value must be a RigidTransform"):
        tf[0] = np.eye(4)


def test_copy_flag():
    # Test that copy=True creates new memory
    matrix = np.eye(4)
    tf = RigidTransform(matrix, normalize=False, copy=True)
    matrix[0, 0] = 2
    assert tf.as_matrix()[0, 0] == 1

    # Test that copy=False shares memory
    matrix = np.eye(4)
    tf = RigidTransform(matrix, normalize=False, copy=False)
    matrix[0, 0] = 2
    assert tf.as_matrix()[0, 0] == 2


def test_normalize_dual_quaternion():
    dual_quat = normalize_dual_quaternion(np.zeros((1, 8)))
    assert_allclose(np.linalg.norm(dual_quat[0, :4]), 1.0, atol=1e-12)
    assert_allclose(dual_quat[0, :4] @ dual_quat[0, 4:], 0.0, atol=1e-12)

    rng = np.random.default_rng(103213650)
    dual_quat = rng.normal(size=(1000, 8))
    dual_quat = normalize_dual_quaternion(dual_quat)
    assert_allclose(np.linalg.norm(dual_quat[:, :4], axis=1), 1.0, atol=1e-12)
    assert_allclose(
        np.einsum("ij,ij->i", dual_quat[:, :4], dual_quat[:, 4:]), 0.0, atol=1e-12
    )


def test_empty_transform_construction():
    tf = RigidTransform.from_matrix(np.empty((0, 4, 4)))
    assert len(tf) == 0
    assert not tf.single

    tf = RigidTransform.from_rotation(Rotation.random(0))
    assert len(tf) == 0
    assert not tf.single

    tf = RigidTransform.from_translation(np.empty((0, 3)))
    assert len(tf) == 0
    assert not tf.single

    tf = RigidTransform.from_components(np.empty((0, 3)), Rotation.random(0))
    assert len(tf) == 0
    assert not tf.single

    tf = RigidTransform.from_exp_coords(np.empty((0, 6)))
    assert len(tf) == 0
    assert not tf.single

    tf = RigidTransform.from_dual_quat(np.empty((0, 8)))
    assert len(tf) == 0
    assert not tf.single

    tf = RigidTransform.identity(0)
    assert len(tf) == 0
    assert not tf.single


def test_empty_transform_representation():
    tf = RigidTransform.identity(0)

    assert len(tf.rotation) == 0
    assert tf.translation.shape == (0, 3)

    t, r = tf.as_components()
    assert t.shape == (0, 3)
    assert len(r) == 0

    assert tf.as_matrix().shape == (0, 4, 4)
    assert tf.as_exp_coords().shape == (0, 6)
    assert tf.as_dual_quat().shape == (0, 8)


def test_empty_transform_application():
    tf = RigidTransform.identity(0)

    assert tf.apply(np.zeros((3,))).shape == (0, 3)
    assert tf.apply(np.empty((0, 3))).shape == (0, 3)

    with pytest.raises(ValueError, match="operands could not be broadcast together"):
        tf.apply(np.zeros((2, 3)))


def test_empty_transform_composition():
    tf_empty = RigidTransform.identity(0)
    tf_single = RigidTransform.identity()
    tf_many = RigidTransform.identity(3)

    assert len(tf_empty * tf_empty) == 0
    assert len(tf_empty * tf_single) == 0
    assert len(tf_single * tf_empty) == 0

    with pytest.raises(ValueError, match="Expected equal number of transforms"):
        tf_many * tf_empty

    with pytest.raises(ValueError, match="Expected equal number of transforms"):
        tf_empty * tf_many


def test_empty_transform_concatenation():
    tf_empty = RigidTransform.identity(0)
    tf_single = RigidTransform.identity()
    tf_many = RigidTransform.identity(2)

    assert len(RigidTransform.concatenate([tf_empty, tf_empty])) == 0
    assert len(RigidTransform.concatenate([tf_empty, tf_single])) == 1
    assert len(RigidTransform.concatenate([tf_single, tf_empty])) == 1
    assert len(RigidTransform.concatenate([tf_empty, tf_many])) == 2
    assert len(RigidTransform.concatenate([tf_many, tf_empty])) == 2
    assert len(RigidTransform.concatenate([tf_many, tf_empty, tf_single])) == 3


def test_empty_transform_inv_and_pow():
    tf = RigidTransform.identity(0)
    assert len(tf.inv()) == 0
    assert len(tf**0) == 0
    assert len(tf**1) == 0
    assert len(tf**-1) == 0
    assert len(tf**0.5) == 0


def test_empty_transform_indexing():
    tf_many = RigidTransform.identity(3)
    tf_zero = tf_many[[]]
    assert len(tf_zero) == 0

    assert len(tf_zero[[]]) == 0
    assert len(tf_zero[:5]) == 0  # Slices can go out of bounds.

    with pytest.raises(IndexError):
        tf_zero[0]

    with pytest.raises(IndexError):
        tf_zero[[0, 2]]

    with pytest.raises(IndexError):
        tf_zero[[False, True]]
