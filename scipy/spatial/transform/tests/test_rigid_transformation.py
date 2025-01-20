import pytest

import numpy as np
from numpy.testing import assert_allclose
from scipy.spatial.transform import Rotation, RigidTransformation


def test_repr():
    actual = repr(RigidTransformation.identity())
    expected = """\
RigidTransformation.from_matrix(array([[1., 0., 0., 0.],
                                       [0., 1., 0., 0.],
                                       [0., 0., 1., 0.],
                                       [0., 0., 0., 1.]]))"""
    assert actual == expected

    actual = repr(RigidTransformation.identity(2))
    expected = """\
RigidTransformation.from_matrix(array([[[1., 0., 0., 0.],
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
    tf = RigidTransformation.from_rotation(r)
    assert_allclose(tf.as_matrix(), np.eye(4), atol=atol)
    assert tf.single

    r = Rotation.from_euler('z', 90, degrees=True)
    tf = RigidTransformation.from_rotation(r)
    assert_allclose(tf.as_matrix()[:3, :3], r.as_matrix(), atol=atol)
    assert_allclose(tf.as_matrix()[:3, 3], [0, 0, 0], atol=atol)
    assert_allclose(tf.as_matrix()[3], [0, 0, 0, 1], atol=atol)
    assert tf.single

    # Test multiple rotations
    r = Rotation.from_euler('zyx', [[90, 0, 0], [0, 90, 0]], degrees=True)
    tf = RigidTransformation.from_rotation(r)
    assert_allclose(tf.as_matrix()[:, :3, :3], r.as_matrix(), atol=atol)
    assert_allclose(tf.as_matrix()[:, :3, 3], [[0, 0, 0], [0, 0, 0]], atol=atol)
    assert_allclose(tf.as_matrix()[:, 3], [[0, 0, 0, 1], [0, 0, 0, 1]], atol=atol)
    assert not tf.single


def test_from_translation():
    # Test single translation
    t = np.array([1, 2, 3])
    tf = RigidTransformation.from_translation(t)
    expected = np.eye(4)
    expected[:3, 3] = t
    assert_allclose(tf.as_matrix(), expected)
    assert tf.single

    # Test multiple translations
    t = np.array([[1, 2, 3], [4, 5, 6]])
    tf = RigidTransformation.from_translation(t)
    for i in range(len(t)):
        expected = np.eye(4)
        expected[:3, 3] = t[i]
        assert_allclose(tf.as_matrix()[i], expected)
    assert not tf.single


def test_from_matrix():
    atol = 1e-12

    # Test single transformation matrix
    matrix = np.eye(4)
    matrix[:3, 3] = [1, 2, 3]
    tf = RigidTransformation.from_matrix(matrix)
    assert_allclose(tf.as_matrix(), matrix, atol=atol)
    assert tf.single

    # Test multiple transformation matrices
    matrices = np.array([np.eye(4)]*2)
    matrices[0, :3, 3] = [1, 2, 3]
    matrices[1, :3, 3] = [4, 5, 6]
    tf = RigidTransformation.from_matrix(matrices)
    assert_allclose(tf.as_matrix(), matrices, atol=atol)
    assert not tf.single

    # Test non-1 determinant
    matrix = np.diag([2, 2, 2, 1])
    tf = RigidTransformation.from_matrix(matrix)
    assert_allclose(tf.as_matrix(), np.eye(4), atol=atol)

    # Test non-orthogonal rotation matrix
    matrix = np.array([[1, 1, 0, 0],
                       [0, 1, 0, 0],
                       [0, 0, 1, 0],
                       [0, 0, 0, 1]])
    tf = RigidTransformation.from_matrix(matrix)
    expected = np.array([[0.894427,  0.447214, 0, 0],
                         [-0.447214,  0.894427, 0, 0],
                         [0, 0, 1, 0],
                         [0, 0, 0, 1]])
    assert_allclose(tf.as_matrix(), expected, atol=1e-6)

    # Test invalid matrix
    with pytest.raises(ValueError):
        invalid = np.eye(4)
        invalid[3, 3] = 2  # Invalid last row
        RigidTransformation.from_matrix(invalid)


def test_from_components():
    atol = 1e-12

    # Test single rotation and translation
    r = Rotation.from_euler('zyx', [90, 0, 0], degrees=True)
    t = np.array([1, 2, 3])
    tf = RigidTransformation.from_components(r, t)

    expected = np.zeros((4, 4))
    expected[:3, :3] = r.as_matrix()
    expected[:3, 3] = t
    expected[3, 3] = 1
    assert_allclose(tf.as_matrix(), expected, atol=atol)
    assert tf.single

    # Test single rotation and multiple translations
    r = Rotation.from_euler('z', 90, degrees=True)
    t = np.array([[1, 2, 3], [4, 5, 6]])
    tf = RigidTransformation.from_components(r, t)
    assert not tf.single

    for i in range(len(t)):
        expected = np.zeros((4, 4))
        expected[:3, :3] = r.as_matrix()
        expected[:3, 3] = t[i]
        expected[3, 3] = 1
        assert_allclose(tf.as_matrix()[i], expected, atol=atol)

    # Test multiple rotations and translations
    r = Rotation.from_euler('zyx', [[90, 0, 0], [0, 90, 0]], degrees=True)
    t = np.array([[1, 2, 3], [4, 5, 6]])
    tf = RigidTransformation.from_components(r, t)
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
    r = Rotation.random(n, rng=10)
    t = np.array([1, 2, 3] * n).reshape(n, 3)
    tf = RigidTransformation.from_components(r, t)
    new_r, new_trans = tf.as_components()
    assert all(new_r.approx_equal(r, atol=atol))
    assert_allclose(new_trans, t, atol=atol)


def test_from_exp_coords():
    # example from 3.3 of
    # https://hades.mech.northwestern.edu/images/2/25/MR-v2.pdf
    angle1 = np.deg2rad(30.0)
    tf1 = RigidTransformation.from_matrix([
        [np.cos(angle1), -np.sin(angle1), 0.0, 1.0],
        [np.sin(angle1), np.cos(angle1), 0.0, 2.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0]
    ])
    angle2 = np.deg2rad(60.0)
    tf2 = RigidTransformation.from_matrix([
        [np.cos(angle2), -np.sin(angle2), 0.0, 2.0],
        [np.sin(angle2), np.cos(angle2), 0.0, 1.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0]
    ])
    expected = tf2 * tf1.inv()
    actual = RigidTransformation.from_exp_coords(
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
        RigidTransformation.from_exp_coords(exp_coords).as_matrix(),
        expected_matrix, atol=1e-8)

    # identity
    assert_allclose(
        RigidTransformation.from_exp_coords(np.zeros(6)).as_matrix(),
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
    actual = RigidTransformation.from_exp_coords([
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
    actual = RigidTransformation.from_exp_coords(
        np.hstack((rotvec, np.zeros((2, 3)))))
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)


def test_as_exp_coords():
    # identity
    expected = np.zeros(6)
    actual = RigidTransformation.from_exp_coords(expected).as_exp_coords()
    assert_allclose(actual, expected, atol=1e-12)

    rng = np.random.default_rng(10)

    # pure rotation
    for _ in range(10):
        expected = np.hstack((rng.normal(size=3), np.zeros(3)))
        actual = RigidTransformation.from_exp_coords(expected).as_exp_coords()
        assert_allclose(actual, expected, atol=1e-12)

    # pure translation
    for _ in range(10):
        expected = np.hstack((np.zeros(3), rng.normal(size=3)))
        actual = RigidTransformation.from_exp_coords(expected).as_exp_coords()
        assert_allclose(actual, expected, atol=1e-12)

    # rotation and translation
    for _ in range(10):
        expected = rng.normal(size=6)
        actual = RigidTransformation.from_exp_coords(expected).as_exp_coords()
        assert_allclose(actual, expected, atol=1e-12)


def test_from_dual_quat():
    # identity
    assert_allclose(
        RigidTransformation.from_dual_quat(
            np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0])).as_matrix(),
        np.eye(4), atol=1e-12)
    assert_allclose(
        RigidTransformation.from_dual_quat(
            np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            scalar_first=True).as_matrix(),
        np.eye(4), atol=1e-12)

    # only translation
    actual = RigidTransformation.from_dual_quat(
        [0, 0, 0, 1, 0.25, 0.15, -0.7, 0])
    expected_matrix = np.array([
        [1, 0, 0, 0.5],
        [0, 1, 0, 0.3],
        [0, 0, 1, -1.4],
        [0, 0, 0, 1]
    ])
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)
    actual = RigidTransformation.from_dual_quat(
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
    actual = RigidTransformation.from_dual_quat(
        np.hstack((actual_rot.as_quat(), np.zeros(4))))
    expected_matrix = np.eye(4)
    expected_matrix[:3, :3] = actual_rot.as_matrix()
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)

    actual = RigidTransformation.from_dual_quat(
        np.hstack((actual_rot.as_quat(scalar_first=True), np.zeros(4))),
        scalar_first=True)
    expected_matrix = np.eye(4)
    expected_matrix[:3, :3] = actual_rot.as_matrix()
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)

    # rotation and translation
    actual = RigidTransformation.from_dual_quat(
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

    actual = RigidTransformation.from_dual_quat(
        [[0.94508498, 0.0617101, -0.06483886, 0.31432811,
          -0.07743254, 0.04985168, -0.26119618, 0.1691491],
         [0.8450749, 0.19507259, 0.49404931, -0.06091285,
          0.04174544, 0.65049656, -0.30782513, 0.16566752]],
        scalar_first=True)
    assert_allclose(actual.as_matrix(), expected_matrix, atol=1e-12)


def test_as_dual_quat():
    # identity
    expected = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0])
    actual = RigidTransformation.identity().as_dual_quat()
    assert_allclose(actual, expected, atol=1e-12)

    expected = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    actual = RigidTransformation.identity().as_dual_quat(scalar_first=True)
    assert_allclose(actual, expected, atol=1e-12)

    rng = np.random.default_rng(10)

    # only rotation
    for _ in range(10):
        real_part = Rotation.random(rng=rng).as_quat()
        dual_part = np.zeros(4)
        expected = np.hstack((real_part, dual_part))
        actual = RigidTransformation.from_dual_quat(expected).as_dual_quat()
        # because of double cover:
        if np.sign(expected[0]) != np.sign(actual[0]):
            actual *= -1.0
        assert_allclose(actual, expected, atol=1e-12)

    # only translation
    for _ in range(10):
        tf = rng.normal(size=3)
        expected = np.hstack(([0, 0, 0, 1], 0.5 * tf, [0]))
        actual = RigidTransformation.from_dual_quat(expected).as_dual_quat()
        # because of double cover:
        if np.sign(expected[0]) != np.sign(actual[0]):
            actual *= -1.0
        assert_allclose(actual, expected, atol=1e-12)

    # rotation and translation
    for _ in range(10):
        expected = RigidTransformation.random(rng=rng).as_dual_quat()
        actual = RigidTransformation.from_dual_quat(expected).as_dual_quat()
        # because of double cover:
        if np.sign(expected[0]) != np.sign(actual[0]):
            actual *= -1.0
        assert_allclose(actual, expected, atol=1e-12)


def test_from_as_internal_consistency():
    atol = 1e-12
    n = 100
    tf0 = RigidTransformation.random(n, rng=10)

    tf1 = RigidTransformation.from_components(*tf0.as_components())
    assert_allclose(tf0.as_matrix(), tf1.as_matrix(), atol=atol)

    tf1 = RigidTransformation.from_components(tf0.rotation, tf0.translation)
    assert_allclose(tf0.as_matrix(), tf1.as_matrix(), atol=atol)

    tf1 = RigidTransformation.from_exp_coords(tf0.as_exp_coords())
    assert_allclose(tf0.as_matrix(), tf1.as_matrix(), atol=atol)

    tf1 = RigidTransformation.from_matrix(tf0.as_matrix())
    assert_allclose(tf0.as_matrix(), tf1.as_matrix(), atol=atol)

    tf1 = RigidTransformation.from_dual_quat(tf0.as_dual_quat())
    assert_allclose(tf0.as_matrix(), tf1.as_matrix(), atol=atol)


def test_identity():
    atol = 1e-12

    # Test single identity
    tf = RigidTransformation.identity()
    assert_allclose(tf.as_matrix(), np.eye(4), atol=atol)

    # Test multiple identities
    tf = RigidTransformation.identity(5)
    assert_allclose(tf.as_matrix(), np.array([np.eye(4)] * 5), atol=atol)


def test_random():
    atol = 1e-6
    n = 5
    tf = RigidTransformation.random(n, rng=0)
    assert len(tf) == n

    # Test repeatability
    expected = np.array([[-0.88199 , -0.36879 ,  0.293406, -0.885639],
                         [ 0.22258 , -0.874756, -0.430418, -0.334283],
                         [ 0.415393, -0.314318,  0.853612, -0.322334],
                         [ 0.      ,  0.      ,  0.      ,  1.      ]])
    assert_allclose(tf[0].as_matrix(), expected, atol=atol)

    # Test that the transformations are valid
    atol = 1e-12
    for m in tf.as_matrix():
        assert_allclose(np.linalg.det(m[:3, :3]), 1, atol=atol)
        assert_allclose(np.linalg.norm(m[:3, 3]), 1, atol=atol)
        assert_allclose(m[3], [0, 0, 0, 1], atol=atol)


def test_apply():
    atol = 1e-12

    ## Single transformation
    r = Rotation.from_euler('z', 90, degrees=True)
    t = np.array([2, 3, 4])
    tf = RigidTransformation.from_components(r, t)

    # Single vector, single transformation
    vec = np.array([1, 0, 0])
    expected = t + r.apply(vec)
    res = tf.apply(vec)
    assert_allclose(res, expected, atol=atol)

    # Multiple vectors, single transformation
    vecs = np.array([[1, 0, 0], [0, 1, 0]])
    expected = t + r.apply(vecs)
    assert_allclose(tf.apply(vecs), expected, atol=atol)

    ## Multiple transformations
    r = Rotation.from_euler('z', [90, 0], degrees=True)
    t = np.array([[2, 3, 4], [5, 6, 7]])
    tf = RigidTransformation.from_components(r, t)

    # Single vector, multiple transformations
    vec = np.array([1, 0, 0])
    expected = t + r.apply(vec)
    assert_allclose(tf.apply(vec), expected, atol=atol)

    # Multiple vectors, multiple transformations
    vecs = np.array([[1, 0, 0], [0, 1, 0]])
    expected = t + r.apply(vecs)
    assert_allclose(tf.apply(vecs), expected, atol=atol)


def test_inverse_apply():
    atol = 1e-12

    # Test applying inverse transformation
    r = Rotation.from_euler('z', 90, degrees=True)
    t = np.array([1, 2, 3])
    tf = RigidTransformation.from_components(r, t)

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
    tf = RigidTransformation.from_rotation(r)
    vec = np.array([1, 0, 0])
    expected = r.apply(vec)
    assert_allclose(tf.apply(vec), expected, atol=atol)


def test_translation_alone():
    atol = 1e-12
    t = np.array([1, 2, 3])
    tf = RigidTransformation.from_translation(t)
    vec = np.array([5, 6, 7])
    expected = t + vec
    assert_allclose(tf.apply(vec), expected, atol=atol)


def test_composition():
    atol = 1e-12

    # Test composing single transformations
    r1 = Rotation.from_euler('z', 90, degrees=True)
    t1 = np.array([1, 0, 0])
    tf1 = RigidTransformation.from_components(r1, t1)

    r2 = Rotation.from_euler('x', 90, degrees=True)
    t2 = np.array([0, 1, 0])
    tf2 = RigidTransformation.from_components(r2, t2)

    composed = tf2 * tf1
    vec = np.array([1, 0, 0])
    expected = tf2.apply(tf1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)
    assert composed.single

    expected = t2 + r2.apply(t1 + r1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)

    # Multiple transformations with single transformation
    t2 = np.array([[1, 2, 3], [4, 5, 6]])
    tf2 = RigidTransformation.from_components(r2, t2)

    composed = tf2 * tf1
    expected = tf2.apply(tf1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)
    assert not composed.single

    expected = t2 + r2.apply(t1 + r1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)

    # Multiple transformations with multiple transformations
    t1 = np.array([[1, 0, 0], [0, -1, 1]])
    tf1 = RigidTransformation.from_components(r1, t1)

    composed = tf2 * tf1
    expected = tf2.apply(tf1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)
    assert not composed.single

    expected = t2 + r2.apply(t1 + r1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)


def test_pow():
    atol = 1e-12
    num = 10
    p = RigidTransformation.random(num, rng=0)
    p_inv = p.inv()

    # Test the short-cuts and other integers
    for n in [-5, -2, -1, 0, 1, 2, 5]:
        q = p**n
        r = RigidTransformation.identity(num)
        for _ in range(abs(n)):
            if n > 0:
                r = r * p
            else:
                r = r * p_inv
        assert_allclose(q.as_matrix(), r.as_matrix(), atol=atol)

        # Test shape preservation
        r = RigidTransformation.from_rotation(Rotation.from_quat([0, 0, 0, 1]))
        assert (r**n).as_matrix().shape == (4, 4)
        r = RigidTransformation.from_rotation(Rotation.from_quat([[0, 0, 0, 1]]))
        assert (r**n).as_matrix().shape == (1, 4, 4)

    # pow function
    tf = pow(RigidTransformation.identity(), 2)
    assert_allclose(tf.as_matrix(), np.eye(4), atol=atol)


def test_inverse():
    atol = 1e-12

    # Test inverse transformation
    r = Rotation.from_euler('z', 90, degrees=True)
    t = np.array([1, 2, 3])
    tf = RigidTransformation.from_components(r, t)

    # Test that tf * tf.inv() equals identity
    tf_inv = tf.inv()
    composed = tf * tf_inv
    assert_allclose(composed.as_matrix(), np.eye(4), atol=atol)

    n = 10
    tf = RigidTransformation.random(n, rng=1)
    tf_inv = tf.inv()
    composed = tf * tf_inv
    assert_allclose(composed.as_matrix(), np.array([np.eye(4)] * n), atol=atol)

    # Test multiple transformations
    r = Rotation.from_euler('zyx', [[90, 0, 0], [0, 90, 0]], degrees=True)
    t = np.array([[1, 2, 3], [4, 5, 6]])
    tf = RigidTransformation.from_components(r, t)
    tf_inv = tf.inv()
    composed = tf * tf_inv
    assert_allclose(composed.as_matrix(), np.array([np.eye(4)] * 2), atol=atol)


def test_properties():
    atol = 1e-12

    # Test rotation and translation properties for single transformation
    r = Rotation.from_euler('z', 90, degrees=True)
    t = np.array([1, 2, 3])
    tf = RigidTransformation.from_components(r, t)

    assert_allclose(tf.rotation.as_matrix(), r.as_matrix(), atol=atol)
    assert tf.rotation.approx_equal(r)
    assert_allclose(tf.translation, t, atol=atol)

    # Test rotation and translation properties for multiple transformations
    r = Rotation.from_euler('zyx', [[90, 0, 0], [0, 90, 0]], degrees=True)
    t = np.array([[1, 2, 3], [4, 5, 6]])
    tf = RigidTransformation.from_components(r, t)

    assert_allclose(tf.rotation.as_matrix(), r.as_matrix(), atol=atol)
    assert all(tf.rotation.approx_equal(r))
    assert_allclose(tf.translation, t, atol=atol)


def test_indexing():
    atol = 1e-12

    # Test indexing for multiple transformations
    r = Rotation.from_euler('zyx', [[90, 0, 0], [0, 90, 0]], degrees=True)
    t = np.array([[1, 2, 3], [4, 5, 6]])
    tf = RigidTransformation.from_components(r, t)

    # Test single index
    assert_allclose(tf[0].as_matrix()[:3, :3], r[0].as_matrix(), atol=atol)
    assert_allclose(tf[0].as_matrix()[:3, 3], t[0], atol=atol)

    # Test slice
    tf_slice = tf[0:2]
    assert_allclose(tf_slice.as_matrix()[:, :3, :3], r[0:2].as_matrix(), atol=atol)
    assert_allclose(tf_slice.as_matrix()[:, :3, 3], t[0:2], atol=atol)


def test_concatenate():
    atol = 1e-12

    # Test concatenation of transformations
    r1 = Rotation.from_euler('z', 90, degrees=True)
    t1 = np.array([1, 0, 0])
    tf1 = RigidTransformation.from_components(r1, t1)

    r2 = Rotation.from_euler('x', 90, degrees=True)
    t2 = np.array([0, 1, 0])
    tf2 = RigidTransformation.from_components(r2, t2)

    # Concatenate single transformations
    concatenated1 = RigidTransformation.concatenate([tf1, tf2])
    assert_allclose(concatenated1[0].as_matrix(), tf1.as_matrix(), atol=atol)
    assert_allclose(concatenated1[1].as_matrix(), tf2.as_matrix(), atol=atol)

    # Concatenate multiple transformations
    concatenated2 = RigidTransformation.concatenate([tf1, concatenated1])
    assert_allclose(concatenated2[0].as_matrix(), tf1.as_matrix(), atol=atol)
    assert_allclose(concatenated2[1].as_matrix(), tf1.as_matrix(), atol=atol)
    assert_allclose(concatenated2[2].as_matrix(), tf2.as_matrix(), atol=atol)


def test_input_validation():
    # Test invalid matrix shapes
    with pytest.raises(ValueError, match="Expected `matrix` to have shape"):
        RigidTransformation.from_matrix(np.eye(3))

    with pytest.raises(ValueError, match="Expected `matrix` to have shape"):
        RigidTransformation.from_matrix(np.zeros((4, 3)))

    with pytest.raises(ValueError, match="Expected `matrix` to have shape"):
        RigidTransformation.from_matrix([])

    with pytest.raises(ValueError, match="Expected `matrix` to have shape"):
        RigidTransformation.from_matrix(np.zeros((0, 4, 4)))

    with pytest.raises(ValueError, match="Expected `matrix` to have shape"):
        RigidTransformation.from_matrix(np.zeros((1, 1, 4, 4)))

    # Test invalid last row
    with pytest.raises(ValueError, match="Expected last row.*to be"):
        matrix = np.eye(4)
        matrix[3, :] = [1, 0, 0, 1]
        RigidTransformation.from_matrix(matrix)

    # Test invalid last row for multiple transformations
    with pytest.raises(ValueError, match="Expected last row.*to be"):
        matrix = np.array([np.eye(4)] * 2)
        matrix[1, 3, :] = [1, 0, 0, 1]
        RigidTransformation.from_matrix(matrix)

    # Test non-rotation matrix
    with pytest.raises(ValueError,
                       match="matrix 0 be orthonormal:"):
        matrix = np.eye(4)
        matrix[:3, :3] *= 2
        RigidTransformation(matrix, normalize=False)

    # Test left handed rotation matrix
    with pytest.raises(ValueError,
                       match="Non-positive determinant"):
        matrix = np.eye(4)
        matrix[0, 0] = -1
        RigidTransformation(matrix, normalize=True)

    # Test non-Rotation input
    with pytest.raises(ValueError,
                       match="Expected `rotation` to be a `Rotation` instance"):
        RigidTransformation.from_rotation(np.eye(3))


def test_translation_validation():
    # Test invalid translation shapes
    with pytest.raises(ValueError, match="Expected `translation` to have shape"):
        RigidTransformation.from_translation([1, 2])

    with pytest.raises(ValueError, match="Expected `translation` to have shape"):
        RigidTransformation.from_translation(np.zeros((2, 2)))

    with pytest.raises(ValueError, match="Expected `translation` to have shape"):
        RigidTransformation.from_translation(np.zeros((0, 3)))

    with pytest.raises(ValueError, match="Expected `translation` to have shape"):
        RigidTransformation.from_translation(np.zeros((1, 1, 3)))


def test_vector_validation():
    tf = RigidTransformation.identity(2)

    # Test invalid vector shapes
    with pytest.raises(ValueError, match="Expected vector to have shape"):
        tf.apply([1, 2])

    with pytest.raises(ValueError, match="Expected vector to have shape"):
        tf.apply(np.zeros((2, 2)))

    with pytest.raises(ValueError, match="Expected vector to have shape"):
        tf.apply(np.zeros((0, 3)))

    with pytest.raises(ValueError, match="Expected vector to have shape"):
        tf.apply(np.zeros((1, 1, 3)))


def test_indexing_validation():
    tf = RigidTransformation.identity()

    # Test indexing on single transformation
    with pytest.raises(TypeError, match="Single transformation is not subscriptable"):
        tf[0]

    with pytest.raises(TypeError, match="Single transformation is not subscriptable"):
        tf[0:1]

    # Test length on single transformation
    with pytest.raises(TypeError, match="Single transformation has no len"):
        len(tf)


def test_composition_validation():
    tf2 = RigidTransformation.from_translation([[1, 2, 3], [4, 5, 6]])
    tf3 = RigidTransformation.from_translation([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    # Test incompatible shapes
    with pytest.raises(ValueError, match="Expected equal number of transformations"):
        tf2 * tf3


def test_concatenate_validation():
    tf = RigidTransformation.identity()

    # Test invalid inputs
    with pytest.raises(TypeError,
                       match="input must contain RigidTransformation objects"):
        RigidTransformation.concatenate([tf, np.eye(4)])


def test_setitem_validation():
    tf = RigidTransformation.from_translation([[1, 2, 3], [4, 5, 6]])
    single = RigidTransformation.identity()

    # Test setting item on single transformation
    with pytest.raises(TypeError, match="Single transformation is not subscriptable"):
        single[0] = tf

    # Test invalid value type
    with pytest.raises(TypeError, match="value must be a RigidTransformation"):
        tf[0] = np.eye(4)


def test_copy_flag():
    # Test that copy=True creates new memory
    matrix = np.eye(4)
    tf = RigidTransformation(matrix, normalize=False, copy=True)
    matrix[0, 0] = 2
    assert tf.as_matrix()[0, 0] == 1

    # Test that copy=False shares memory
    matrix = np.eye(4)
    tf = RigidTransformation(matrix, normalize=False, copy=False)
    matrix[0, 0] = 2
    assert tf.as_matrix()[0, 0] == 2
