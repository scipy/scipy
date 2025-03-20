import pytest

import numpy as np
from numpy.testing import assert_allclose
from scipy.spatial.transform import Rotation, RigidTransform
from scipy.spatial.transform._rigid_transform import normalize_dual_quaternion


def test_repr():
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


def test_from_rotation():
    atol = 1e-12

    # Test single rotation
    r = Rotation.identity()
    tf = RigidTransform.from_rotation(r)
    assert_allclose(tf.as_matrix(), np.eye(4), atol=atol)
    assert tf.single

    r = Rotation.from_euler('z', 90, degrees=True)
    tf = RigidTransform.from_rotation(r)
    assert_allclose(tf.as_matrix()[:3, :3], r.as_matrix(), atol=atol)
    assert_allclose(tf.as_matrix()[:3, 3], [0, 0, 0], atol=atol)
    assert_allclose(tf.as_matrix()[3], [0, 0, 0, 1], atol=atol)
    assert tf.single

    # Test multiple rotations
    r = Rotation.from_euler('zyx', [[90, 0, 0], [0, 90, 0]], degrees=True)
    tf = RigidTransform.from_rotation(r)
    assert_allclose(tf.as_matrix()[:, :3, :3], r.as_matrix(), atol=atol)
    assert_allclose(tf.as_matrix()[:, :3, 3], [[0, 0, 0], [0, 0, 0]], atol=atol)
    assert_allclose(tf.as_matrix()[:, 3], [[0, 0, 0, 1], [0, 0, 0, 1]], atol=atol)
    assert not tf.single


def test_from_translation():
    # Test single translation
    t = np.array([1, 2, 3])
    tf = RigidTransform.from_translation(t)
    expected = np.eye(4)
    expected[:3, 3] = t
    assert_allclose(tf.as_matrix(), expected)
    assert tf.single

    # Test multiple translations
    t = np.array([[1, 2, 3], [4, 5, 6]])
    tf = RigidTransform.from_translation(t)
    for i in range(len(t)):
        expected = np.eye(4)
        expected[:3, 3] = t[i]
        assert_allclose(tf.as_matrix()[i], expected)
    assert not tf.single


def test_from_matrix():
    atol = 1e-12

    # Test single transform matrix
    matrix = np.eye(4)
    matrix[:3, 3] = [1, 2, 3]
    tf = RigidTransform.from_matrix(matrix)
    assert_allclose(tf.as_matrix(), matrix, atol=atol)
    assert tf.single

    # Test multiple transform matrices
    matrices = np.array([np.eye(4)]*2)
    matrices[0, :3, 3] = [1, 2, 3]
    matrices[1, :3, 3] = [4, 5, 6]
    tf = RigidTransform.from_matrix(matrices)
    assert_allclose(tf.as_matrix(), matrices, atol=atol)
    assert not tf.single

    # Test non-1 determinant
    matrix = np.diag([2, 2, 2, 1])
    tf = RigidTransform.from_matrix(matrix)
    assert_allclose(tf.as_matrix(), np.eye(4), atol=atol)

    # Test non-orthogonal rotation matrix
    matrix = np.array([[1, 1, 0, 0],
                       [0, 1, 0, 0],
                       [0, 0, 1, 0],
                       [0, 0, 0, 1]])
    tf = RigidTransform.from_matrix(matrix)
    expected = np.array([[0.894427,  0.447214, 0, 0],
                         [-0.447214,  0.894427, 0, 0],
                         [0, 0, 1, 0],
                         [0, 0, 0, 1]])
    assert_allclose(tf.as_matrix(), expected, atol=1e-6)

    # Test invalid matrix
    with pytest.raises(ValueError):
        invalid = np.eye(4)
        invalid[3, 3] = 2  # Invalid last row
        RigidTransform.from_matrix(invalid)


def test_from_components():
    atol = 1e-12

    # Test single rotation and translation
    t = np.array([1, 2, 3])
    r = Rotation.from_euler('zyx', [90, 0, 0], degrees=True)
    tf = RigidTransform.from_components(t, r)

    expected = np.zeros((4, 4))
    expected[:3, :3] = r.as_matrix()
    expected[:3, 3] = t
    expected[3, 3] = 1
    assert_allclose(tf.as_matrix(), expected, atol=atol)
    assert tf.single

    # Test single rotation and multiple translations
    t = np.array([[1, 2, 3], [4, 5, 6]])
    r = Rotation.from_euler('z', 90, degrees=True)
    tf = RigidTransform.from_components(t, r)
    assert not tf.single

    for i in range(len(t)):
        expected = np.zeros((4, 4))
        expected[:3, :3] = r.as_matrix()
        expected[:3, 3] = t[i]
        expected[3, 3] = 1
        assert_allclose(tf.as_matrix()[i], expected, atol=atol)

    # Test multiple rotations and translations
    t = np.array([[1, 2, 3], [4, 5, 6]])
    r = Rotation.from_euler('zyx', [[90, 0, 0], [0, 90, 0]], degrees=True)
    tf = RigidTransform.from_components(t, r)
    assert not tf.single

    for i in range(len(t)):
        expected = np.zeros((4, 4))
        expected[:3, :3] = r[i].as_matrix()
        expected[:3, 3] = t[i]
        expected[3, 3] = 1
        assert_allclose(tf.as_matrix()[i], expected, atol=atol)


def test_as_components():
    atol = 1e-12
    n = 10
    rng = np.random.default_rng(123)
    t = rng.normal(size=(n, 3))
    r = Rotation.random(n, rng=rng)
    tf = RigidTransform.from_components(t, r)
    new_t, new_r = tf.as_components()
    assert all(new_r.approx_equal(r, atol=atol))
    assert_allclose(new_t, t, atol=atol)


def test_from_exp_coords():
    # example from 3.3 of
    # https://hades.mech.northwestern.edu/images/2/25/MR-v2.pdf
    angle1 = np.deg2rad(30.0)
    tf1 = RigidTransform.from_matrix([
        [np.cos(angle1), -np.sin(angle1), 0.0, 1.0],
        [np.sin(angle1), np.cos(angle1), 0.0, 2.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0]
    ])
    angle2 = np.deg2rad(60.0)
    tf2 = RigidTransform.from_matrix([
        [np.cos(angle2), -np.sin(angle2), 0.0, 2.0],
        [np.sin(angle2), np.cos(angle2), 0.0, 1.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0]
    ])
    expected = tf2 * tf1.inv()
    actual = RigidTransform.from_exp_coords(
        np.deg2rad(30.0) * np.array([0.0, 0.0, 1.0, 3.37, -3.37, 0.0]))
    assert_allclose(actual.as_matrix(), expected.as_matrix(), atol=1e-2)

    # test cases generated by comparison to pytransform3d
    exp_coords = [
        [-2.01041204, -0.52983629, 0.65773501,
         0.10386614, 0.05855009, 0.54959179],
        [-0.22537438, -0.24132627, -2.4747121,
         -0.09158594,  1.88075832, -0.03197204]
    ]
    expected_matrix = [
        [[0.76406621, 0.10504613, -0.63652819, -0.10209961],
         [0.59956454, -0.47987325, 0.64050295, 0.40158789],
         [-0.2381705, -0.87102639, -0.42963687, 0.19637636],
         [0., 0., 0., 1.]],
        [[-0.78446989, 0.61157488, 0.10287448, 1.33330055],
         [-0.58017785, -0.78232107, 0.22664378, 0.52660831],
         [0.21909052, 0.11810973, 0.96852952, -0.02968529],
         [0., 0., 0., 1.]]
    ]
    assert_allclose(
        RigidTransform.from_exp_coords(exp_coords).as_matrix(),
        expected_matrix, atol=1e-8)

    # identity
    assert_allclose(
        RigidTransform.from_exp_coords(np.zeros(6)).as_matrix(),
        np.eye(4), atol=1e-12)

    # only translation
    expected_matrix = np.array([
        [[1.0, 0.0, 0.0, 3.0],
         [0.0, 1.0, 0.0, -5.4],
         [0.0, 0.0, 1.0, 100.2],
         [0.0, 0.0, 0.0, 1.0]],
        [[1.0, 0.0, 0.0, -3.0],
         [0.0, 1.0, 0.0, 13.3],
         [0.0, 0.0, 1.0, 1.3],
         [0.0, 0.0, 0.0, 1.0]]
    ])
    actual = RigidTransform.from_exp_coords([
        [0.0, 0.0, 0.0, 3.0, -5.4, 100.2],
        [0.0, 0.0, 0.0, -3.0, 13.3, 1.3],
    ])
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)

    # only rotation
    rot = Rotation.from_euler(
        'zyx',
        [[34, -12, 0.5],
         [-102, -55, 30]],
        degrees=True)
    rotvec = rot.as_rotvec()
    expected_matrix = np.array([np.eye(4), np.eye(4)])
    expected_matrix[:, :3, :3] = rot.as_matrix()
    actual = RigidTransform.from_exp_coords(
        np.hstack((rotvec, np.zeros((2, 3)))))
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)


def test_as_exp_coords():
    # identity
    expected = np.zeros(6)
    actual = RigidTransform.from_exp_coords(expected).as_exp_coords()
    assert_allclose(actual, expected, atol=1e-12)

    rng = np.random.default_rng(10)

    # pure rotation
    rot_vec = rng.normal(scale=0.1, size=(1000, 3))
    tf = RigidTransform.from_rotation(Rotation.from_rotvec(rot_vec))
    exp_coords = tf.as_exp_coords()
    assert_allclose(exp_coords[:, :3], rot_vec, rtol=1e-13)
    assert_allclose(exp_coords[:, 3:], 0.0, atol=1e-16)

    # pure translation
    translation = rng.normal(scale=100.0, size=(1000, 3))
    tf = RigidTransform.from_translation(translation)
    exp_coords = tf.as_exp_coords()
    assert_allclose(exp_coords[:, :3], 0.0, atol=1e-16)
    assert_allclose(exp_coords[:, 3:], translation, rtol=1e-15)


def test_from_dual_quat():
    # identity
    assert_allclose(
        RigidTransform.from_dual_quat(
            np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0])).as_matrix(),
        np.eye(4), atol=1e-12)
    assert_allclose(
        RigidTransform.from_dual_quat(
            np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            scalar_first=True).as_matrix(),
        np.eye(4), atol=1e-12)

    # only translation
    actual = RigidTransform.from_dual_quat(
        [0, 0, 0, 1, 0.25, 0.15, -0.7, 0])
    expected_matrix = np.array([
        [1, 0, 0, 0.5],
        [0, 1, 0, 0.3],
        [0, 0, 1, -1.4],
        [0, 0, 0, 1]
    ])
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)
    actual = RigidTransform.from_dual_quat(
        [1, 0, 0, 0, 0, 0.25, 0.15, -0.7], scalar_first=True)
    expected_matrix = np.array([
        [1, 0, 0, 0.5],
        [0, 1, 0, 0.3],
        [0, 0, 1, -1.4],
        [0, 0, 0, 1]
    ])
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)

    # only rotation
    actual_rot = Rotation.from_euler("xyz", [65, -13, 90], degrees=True)
    actual = RigidTransform.from_dual_quat(
        np.hstack((actual_rot.as_quat(), np.zeros(4))))
    expected_matrix = np.eye(4)
    expected_matrix[:3, :3] = actual_rot.as_matrix()
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)

    actual = RigidTransform.from_dual_quat(
        np.hstack((actual_rot.as_quat(scalar_first=True), np.zeros(4))),
        scalar_first=True)
    expected_matrix = np.eye(4)
    expected_matrix[:3, :3] = actual_rot.as_matrix()
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)

    # rotation and translation
    actual = RigidTransform.from_dual_quat(
        [[0.0617101, -0.06483886, 0.31432811, 0.94508498,
          0.04985168, -0.26119618, 0.1691491, -0.07743254],
         [0.19507259, 0.49404931, -0.06091285, 0.8450749,
          0.65049656, -0.30782513, 0.16566752, 0.04174544]])
    expected_matrix = np.array(
        [[[0.79398752, -0.60213598, -0.08376202, 0.24605262],
          [0.58613113, 0.79477941, -0.15740392, -0.4932833],
          [0.16135089, 0.07588122, 0.98397557, 0.34262676],
          [0., 0., 0., 1.]],
         [[0.50440981, 0.2957028, 0.81125249, 1.20934468],
          [0.08979911, 0.91647262, -0.3898898, -0.70540077],
          [-0.8587822, 0.26951399, 0.43572393, -0.47776265],
          [0., 0., 0., 1.]]])
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)

    actual = RigidTransform.from_dual_quat(
        [[0.94508498, 0.0617101, -0.06483886, 0.31432811,
          -0.07743254, 0.04985168, -0.26119618, 0.1691491],
         [0.8450749, 0.19507259, 0.49404931, -0.06091285,
          0.04174544, 0.65049656, -0.30782513, 0.16566752]],
        scalar_first=True)
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)

    # unnormalized dual quaternions

    # invalid real quaternion with norm 0
    actual = RigidTransform.from_dual_quat(np.zeros(8))
    assert_allclose(actual.as_matrix(), np.eye(4), atol=1e-12)

    # real quaternion with norm != 1
    unnormalized_dual_quat = np.array(
        [-0.2547655, 1.23506123, 0.20230088, 0.24247194,  # norm 1.3
         0.38559628, 0.08184063, 0.1755943, -0.1582222]  # orthogonal
    )
    assert pytest.approx(np.linalg.norm(unnormalized_dual_quat[:4])) == 1.3
    assert pytest.approx(np.dot(unnormalized_dual_quat[:4],
                                unnormalized_dual_quat[4:]), abs=8) == 0.0
    dual_quat = RigidTransform.from_dual_quat(
        unnormalized_dual_quat).as_dual_quat()
    assert pytest.approx(np.linalg.norm(dual_quat[:4])) == 1.0
    assert pytest.approx(np.dot(dual_quat[:4], dual_quat[4:])) == 0.0

    # real and dual quaternion are not orthogonal
    unnormalized_dual_quat = np.array(
        [0.20824458, 0.75098079, 0.54542913, -0.30849493,  # unit norm
         -0.16051025, 0.10742978, 0.21277201, 0.20596935]  # not orthogonal
    )
    assert pytest.approx(np.linalg.norm(unnormalized_dual_quat[:4])) == 1.0
    assert np.dot(unnormalized_dual_quat[:4],
                  unnormalized_dual_quat[4:]) != 0.0
    dual_quat = RigidTransform.from_dual_quat(
        unnormalized_dual_quat).as_dual_quat()
    assert pytest.approx(np.linalg.norm(dual_quat[:4])) == 1.0
    assert pytest.approx(np.dot(dual_quat[:4], dual_quat[4:])) == 0.0

    # invalid real quaternion with norm 0, non-orthogonal dual quaternion
    unnormalized_dual_quat = np.array(
        [0.0, 0.0, 0.0, 0.0, -0.16051025, 0.10742978, 0.21277201, 0.20596935])
    assert np.dot(np.array([0.0, 0.0, 0.0, 1.0]),
                  unnormalized_dual_quat[4:]) != 0.0
    dual_quat = RigidTransform.from_dual_quat(
        unnormalized_dual_quat).as_dual_quat()
    assert_allclose(dual_quat[:4], np.array([0, 0, 0, 1]), atol=1e-12)
    assert pytest.approx(np.dot(dual_quat[:4], dual_quat[4:])) == 0.0

    # compensation for precision loss in real quaternion
    rng = np.random.default_rng(1000)
    t = rng.normal(size=(3,))
    r = Rotation.random(10, rng=rng)
    random_dual_quats = RigidTransform.from_components(t, r).as_dual_quat()

    # ensure that random quaternions are not normalized
    random_dual_quats[:, :4] = random_dual_quats[:, :4].round(2)
    assert not np.any(np.isclose(
        np.linalg.norm(random_dual_quats[:, :4], axis=1), 1.0, atol=0.0001))
    dual_quat_norm = RigidTransform.from_dual_quat(
        random_dual_quats).as_dual_quat()
    assert_allclose(
        np.linalg.norm(dual_quat_norm[:, :4], axis=1), 1.0, atol=1e-12)

    # compensation for precision loss in dual quaternion, results in violation
    # of orthogonality constraint
    t = rng.normal(size=(10, 3))
    r = Rotation.random(10, rng=rng)
    random_dual_quats = RigidTransform.from_components(t, r).as_dual_quat()

    # ensure that random quaternions are not normalized
    random_dual_quats[:, 4:] = random_dual_quats[:, 4:].round(2)
    assert not np.any(np.isclose(
        np.einsum("ij,ij->i",
                  random_dual_quats[:, :4],
                  random_dual_quats[:, 4:]),
        0.0, atol=0.0001))
    dual_quat_norm = RigidTransform.from_dual_quat(
        random_dual_quats).as_dual_quat()
    assert_allclose(
        np.einsum("ij,ij->i", dual_quat_norm[:, :4], dual_quat_norm[:, 4:]),
        0.0, atol=1e-12)
    assert_allclose(
        random_dual_quats[:, :4], dual_quat_norm[:, :4], atol=1e-12)


def test_as_dual_quat():
    # identity
    expected = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0])
    actual = RigidTransform.identity().as_dual_quat()
    assert_allclose(actual, expected, atol=1e-12)

    expected = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    actual = RigidTransform.identity().as_dual_quat(scalar_first=True)
    assert_allclose(actual, expected, atol=1e-12)

    rng = np.random.default_rng(10)

    # only rotation
    for _ in range(10):
        real_part = Rotation.random(rng=rng).as_quat()
        dual_part = np.zeros(4)
        expected = np.hstack((real_part, dual_part))
        actual = RigidTransform.from_dual_quat(expected).as_dual_quat()
        # because of double cover:
        if np.sign(expected[0]) != np.sign(actual[0]):
            actual *= -1.0
        assert_allclose(actual, expected, atol=1e-12)

    # only translation
    for _ in range(10):
        tf = rng.normal(size=3)
        expected = np.hstack(([0, 0, 0, 1], 0.5 * tf, [0]))
        actual = RigidTransform.from_dual_quat(expected).as_dual_quat()
        # because of double cover:
        if np.sign(expected[0]) != np.sign(actual[0]):
            actual *= -1.0
        assert_allclose(actual, expected, atol=1e-12)

    # rotation and translation
    for _ in range(10):
        t = rng.normal(size=3)
        r = Rotation.random(rng=rng)
        expected = RigidTransform.from_components(t, r).as_dual_quat()
        actual = RigidTransform.from_dual_quat(expected).as_dual_quat()
        # because of double cover:
        if np.sign(expected[0]) != np.sign(actual[0]):
            actual *= -1.0
        assert_allclose(actual, expected, atol=1e-12)


def test_from_as_internal_consistency():
    atol = 1e-12
    n = 1000
    rng = np.random.default_rng(10)
    t = rng.normal(size=(n, 3))
    r = Rotation.random(n, rng=rng)
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
    tf0 = RigidTransform.from_components(
        rng.normal(scale=1000.0, size=(1000, 3)),
        Rotation.from_rotvec(rng.normal(scale=1e-10, size=(1000, 3))))
    tf1 = RigidTransform.from_exp_coords(tf0.as_exp_coords())
    assert_allclose(tf0.as_matrix(), tf1.as_matrix(), atol=atol)


def test_identity():
    atol = 1e-12

    # Test single identity
    tf = RigidTransform.identity()
    assert_allclose(tf.as_matrix(), np.eye(4), atol=atol)

    # Test multiple identities
    tf = RigidTransform.identity(5)
    assert_allclose(tf.as_matrix(), np.array([np.eye(4)] * 5), atol=atol)


def test_apply():
    atol = 1e-12

    ## Single transform
    r = Rotation.from_euler('z', 90, degrees=True)
    t = np.array([2, 3, 4])
    tf = RigidTransform.from_components(t, r)

    # Single vector, single transform
    vec = np.array([1, 0, 0])
    expected = t + r.apply(vec)
    res = tf.apply(vec)
    assert_allclose(res, expected, atol=atol)

    # Multiple vectors, single transform
    vecs = np.array([[1, 0, 0], [0, 1, 0]])
    expected = t + r.apply(vecs)
    assert_allclose(tf.apply(vecs), expected, atol=atol)

    ## Multiple transforms
    r = Rotation.from_euler('z', [90, 0], degrees=True)
    t = np.array([[2, 3, 4], [5, 6, 7]])
    tf = RigidTransform.from_components(t, r)

    # Single vector, multiple transforms
    vec = np.array([1, 0, 0])
    expected = t + r.apply(vec)
    assert_allclose(tf.apply(vec), expected, atol=atol)

    # Multiple vectors, multiple transforms
    vecs = np.array([[1, 0, 0], [0, 1, 0]])
    expected = t + r.apply(vecs)
    assert_allclose(tf.apply(vecs), expected, atol=atol)


def test_inverse_apply():
    atol = 1e-12

    # Test applying inverse transform
    t = np.array([1, 2, 3])
    r = Rotation.from_euler('z', 90, degrees=True)
    tf = RigidTransform.from_components(t, r)

    # Test single vector
    vec = np.array([1, 0, 0])
    expected = tf.inv().apply(vec)
    assert_allclose(tf.apply(vec, inverse=True), expected, atol=atol)

    # Test multiple vectors
    vecs = np.array([[1, 0, 0], [0, 1, 0]])
    expected = tf.inv().apply(vecs)
    assert_allclose(tf.apply(vecs, inverse=True), expected, atol=atol)


def test_rotation_alone():
    atol = 1e-12

    r = Rotation.from_euler('z', 90, degrees=True)
    tf = RigidTransform.from_rotation(r)
    vec = np.array([1, 0, 0])
    expected = r.apply(vec)
    assert_allclose(tf.apply(vec), expected, atol=atol)


def test_translation_alone():
    atol = 1e-12
    t = np.array([1, 2, 3])
    tf = RigidTransform.from_translation(t)
    vec = np.array([5, 6, 7])
    expected = t + vec
    assert_allclose(tf.apply(vec), expected, atol=atol)


def test_composition():
    atol = 1e-12

    # Test composing single transforms
    t1 = np.array([1, 0, 0])
    r1 = Rotation.from_euler('z', 90, degrees=True)
    tf1 = RigidTransform.from_components(t1, r1)

    t2 = np.array([0, 1, 0])
    r2 = Rotation.from_euler('x', 90, degrees=True)
    tf2 = RigidTransform.from_components(t2, r2)

    composed = tf2 * tf1
    vec = np.array([1, 0, 0])
    expected = tf2.apply(tf1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)
    assert composed.single

    expected = t2 + r2.apply(t1 + r1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)

    # Multiple transforms with single transform
    t2 = np.array([[1, 2, 3], [4, 5, 6]])
    tf2 = RigidTransform.from_components(t2, r2)

    composed = tf2 * tf1
    expected = tf2.apply(tf1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)
    assert not composed.single

    expected = t2 + r2.apply(t1 + r1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)

    # Multiple transforms with multiple transforms
    t1 = np.array([[1, 0, 0], [0, -1, 1]])
    tf1 = RigidTransform.from_components(t1, r1)

    composed = tf2 * tf1
    expected = tf2.apply(tf1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)
    assert not composed.single

    expected = t2 + r2.apply(t1 + r1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)


def test_pow():
    atol = 1e-12
    num = 10
    rng = np.random.default_rng(100)
    t = rng.normal(size=(num, 3))
    r = Rotation.random(num, rng=rng)
    p = RigidTransform.from_components(t, r)
    p_inv = p.inv()

    # Test the short-cuts and other integers
    for n in [-5, -2, -1, 0, 1, 2, 5]:
        q = p**n
        r = RigidTransform.identity(num)
        for _ in range(abs(n)):
            if n > 0:
                r = r * p
            else:
                r = r * p_inv
        assert_allclose(q.as_matrix(), r.as_matrix(), atol=atol)

        # Test shape preservation
        r = RigidTransform.from_rotation(Rotation.from_quat([0, 0, 0, 1]))
        assert (r**n).as_matrix().shape == (4, 4)
        r = RigidTransform.from_rotation(Rotation.from_quat([[0, 0, 0, 1]]))
        assert (r**n).as_matrix().shape == (1, 4, 4)

    # Test fractional powers
    q = p**0.5
    assert_allclose((q * q).as_matrix(), p.as_matrix(), atol=atol)
    q = p**-0.5
    assert_allclose((q * q).as_matrix(), p.inv().as_matrix(), atol=atol)
    q = p** 1.5
    assert_allclose((q * q).as_matrix(), (p**3).as_matrix(), atol=atol)
    q = p** -1.5
    assert_allclose((q * q).as_matrix(), (p**-3).as_matrix(), atol=atol)

    # pow function
    tf = pow(RigidTransform.identity(), 2)
    assert_allclose(tf.as_matrix(), np.eye(4), atol=atol)


def test_pow_equivalence_with_rotation():
    atol = 1e-12
    num = 10
    rng = np.random.default_rng(100)
    r = Rotation.random(num, rng=rng)
    p = RigidTransform.from_rotation(r)
    for n in [-5, -2, -1.5, -1, -0.5, 0.0, 0.5, 1, 1.5, 2, 5]:
        assert_allclose((p**n).rotation.as_matrix(), (r**n).as_matrix(), atol=atol)


def test_inverse():
    atol = 1e-12

    # Test inverse transform
    r = Rotation.from_euler('z', 90, degrees=True)
    t = np.array([1, 2, 3])
    tf = RigidTransform.from_components(t, r)

    # Test that tf * tf.inv() equals identity
    tf_inv = tf.inv()
    composed = tf * tf_inv
    assert_allclose(composed.as_matrix(), np.eye(4), atol=atol)

    n = 10
    rng = np.random.default_rng(1000)
    t = rng.normal(size=(n, 3))
    r = Rotation.random(n, rng=rng)
    tf = RigidTransform.from_components(t, r)
    tf_inv = tf.inv()
    composed = tf * tf_inv
    assert_allclose(composed.as_matrix(), np.array([np.eye(4)] * n), atol=atol)

    # Test multiple transforms
    r = Rotation.from_euler('zyx', [[90, 0, 0], [0, 90, 0]], degrees=True)
    t = np.array([[1, 2, 3], [4, 5, 6]])
    tf = RigidTransform.from_components(t, r)
    tf_inv = tf.inv()
    composed = tf * tf_inv
    assert_allclose(composed.as_matrix(), np.array([np.eye(4)] * 2), atol=atol)


def test_properties():
    atol = 1e-12

    # Test rotation and translation properties for single transform
    r = Rotation.from_euler('z', 90, degrees=True)
    t = np.array([1, 2, 3])
    tf = RigidTransform.from_components(t, r)

    assert_allclose(tf.rotation.as_matrix(), r.as_matrix(), atol=atol)
    assert tf.rotation.approx_equal(r)
    assert_allclose(tf.translation, t, atol=atol)

    # Test rotation and translation properties for multiple transforms
    r = Rotation.from_euler('zyx', [[90, 0, 0], [0, 90, 0]], degrees=True)
    t = np.array([[1, 2, 3], [4, 5, 6]])
    tf = RigidTransform.from_components(t, r)

    assert_allclose(tf.rotation.as_matrix(), r.as_matrix(), atol=atol)
    assert all(tf.rotation.approx_equal(r))
    assert_allclose(tf.translation, t, atol=atol)


def test_indexing():
    atol = 1e-12

    # Test indexing for multiple transforms
    r = Rotation.from_euler('zyx', [[90, 0, 0], [0, 90, 0]], degrees=True)
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
    assert_allclose(tf_masked.as_matrix()[:, :3, :3], r[[False, True]].as_matrix(),
                    atol=atol)
    assert_allclose(tf_masked.as_matrix()[:, :3, 3], t[[False, True]], atol=atol)

    tf_masked = tf[[False, False]]
    assert len(tf_masked) == 0


def test_concatenate():
    atol = 1e-12

    # Test concatenation of transforms
    t1 = np.array([1, 0, 0])
    r1 = Rotation.from_euler('z', 90, degrees=True)
    tf1 = RigidTransform.from_components(t1, r1)

    t2 = np.array([0, 1, 0])
    r2 = Rotation.from_euler('x', 90, degrees=True)
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
    with pytest.raises(ValueError,
                       match="Expected `rotation` to be a `Rotation` instance"):
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
    with pytest.raises(TypeError,
                       match="input must contain RigidTransform objects"):
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
    assert_allclose(np.einsum("ij,ij->i", dual_quat[:, :4], dual_quat[:, 4:]),
                    0.0, atol=1e-12)


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
    assert len(tf ** 0) == 0
    assert len(tf ** 1) == 0
    assert len(tf ** -1) == 0
    assert len(tf ** 0.5) == 0


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
