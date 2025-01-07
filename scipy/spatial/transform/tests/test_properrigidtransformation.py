import pytest

import numpy as np
from numpy.testing import assert_allclose
from scipy.spatial.transform import ProperRigidTransformation as T
from scipy.spatial.transform import Rotation as R


def test_from_rotation():
    atol = 1e-12

    # Test single rotation
    r = R.identity()
    t = T.from_rotation(r)
    assert_allclose(t.as_matrix(), np.eye(4), atol=atol)
    assert t.single

    r = R.from_euler('z', 90, degrees=True)
    t = T.from_rotation(r)
    assert_allclose(t.as_matrix()[:3, :3], r.as_matrix(), atol=atol)
    assert_allclose(t.as_matrix()[:3, 3], [0, 0, 0], atol=atol)
    assert_allclose(t.as_matrix()[3], [0, 0, 0, 1], atol=atol)
    assert t.single

    # Test multiple rotations
    r = R.from_euler('zyx', [[90, 0, 0], [0, 90, 0]], degrees=True)
    t = T.from_rotation(r)
    assert_allclose(t.as_matrix()[:, :3, :3], r.as_matrix(), atol=atol)
    assert_allclose(t.as_matrix()[:, :3, 3], [[0, 0, 0], [0, 0, 0]], atol=atol)
    assert_allclose(t.as_matrix()[:, 3], [[0, 0, 0, 1], [0, 0, 0, 1]], atol=atol)
    assert not t.single


def test_from_translation():
    # Test single translation
    trans = np.array([1, 2, 3])
    t = T.from_translation(trans)
    expected = np.eye(4)
    expected[:3, 3] = trans
    assert_allclose(t.as_matrix(), expected)
    assert t.single

    # Test multiple translations
    trans = np.array([[1, 2, 3], [4, 5, 6]])
    t = T.from_translation(trans)
    for i in range(len(trans)):
        expected = np.eye(4)
        expected[:3, 3] = trans[i]
        assert_allclose(t.as_matrix()[i], expected)
    assert not t.single


def test_from_matrix():
    atol = 1e-12

    # Test single transformation matrix
    matrix = np.eye(4)
    matrix[:3, 3] = [1, 2, 3]
    t = T.from_matrix(matrix)
    assert_allclose(t.as_matrix(), matrix, atol=atol)
    assert t.single

    # Test multiple transformation matrices
    matrices = np.array([np.eye(4)]*2)
    matrices[0, :3, 3] = [1, 2, 3]
    matrices[1, :3, 3] = [4, 5, 6]
    t = T.from_matrix(matrices)
    assert_allclose(t.as_matrix(), matrices, atol=atol)
    assert not t.single

    # Test invalid matrix
    with pytest.raises(ValueError):
        invalid = np.eye(4)
        invalid[3, 3] = 2  # Invalid last row
        T.from_matrix(invalid)


def test_from_transrot():
    atol = 1e-12

    # Test single rotation and translation
    r = R.from_euler('zyx', [90, 0, 0], degrees=True)
    trans = np.array([1, 2, 3])
    t = T.from_transrot(trans, r)

    expected = np.zeros((4, 4))
    expected[:3, :3] = r.as_matrix()
    expected[:3, 3] = trans
    expected[3, 3] = 1
    assert_allclose(t.as_matrix(), expected, atol=atol)
    assert t.single

    # Test single rotation and multiple translations
    r = R.from_euler('z', 90, degrees=True)
    trans = np.array([[1, 2, 3], [4, 5, 6]])
    t = T.from_transrot(trans, r)
    assert not t.single

    for i in range(len(trans)):
        expected = np.zeros((4, 4))
        expected[:3, :3] = r.as_matrix()
        expected[:3, 3] = trans[i]
        expected[3, 3] = 1
        assert_allclose(t.as_matrix()[i], expected, atol=atol)

    # Test multiple rotations and translations
    r = R.from_euler('zyx', [[90, 0, 0], [0, 90, 0]], degrees=True)
    trans = np.array([[1, 2, 3], [4, 5, 6]])
    t = T.from_transrot(trans, r)
    assert not t.single

    for i in range(len(trans)):
        expected = np.zeros((4, 4))
        expected[:3, :3] = r[i].as_matrix()
        expected[:3, 3] = trans[i]
        expected[3, 3] = 1
        assert_allclose(t.as_matrix()[i], expected, atol=atol)


def test_from_expcoords():
    # example from 3.3 of
    # https://hades.mech.northwestern.edu/images/2/25/MR-v2.pdf
    angle1 = np.deg2rad(30.0)
    T1 = T.from_matrix([
        [np.cos(angle1), -np.sin(angle1), 0.0, 1.0],
        [np.sin(angle1), np.cos(angle1), 0.0, 2.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0]
    ])
    angle2 = np.deg2rad(60.0)
    T2 = T.from_matrix([
        [np.cos(angle2), -np.sin(angle2), 0.0, 2.0],
        [np.sin(angle2), np.cos(angle2), 0.0, 1.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0]
    ])
    T_expected = T2 * T1.inv()
    T_actual = T.from_expcoords(
        np.deg2rad(30.0) * np.array([0.0, 0.0, 1.0, 3.37, -3.37, 0.0]))
    assert_allclose(T_actual.as_matrix(), T_expected.as_matrix(), atol=1e-2)


def test_identity():
    atol = 1e-12

    # Test single identity
    t = T.identity()
    assert_allclose(t.as_matrix(), np.eye(4), atol=atol)

    # Test multiple identities
    t = T.identity(5)
    assert_allclose(t.as_matrix(), np.array([np.eye(4)] * 5), atol=atol)


def test_random():
    atol = 1e-6
    n = 5
    t = T.random(n, rng=0)
    assert len(t) == n

    # Test repeatability
    expected = np.array([[-0.88199 , -0.36879 ,  0.293406, -0.885639],
                         [ 0.22258 , -0.874756, -0.430418, -0.334283],
                         [ 0.415393, -0.314318,  0.853612, -0.322334],
                         [ 0.      ,  0.      ,  0.      ,  1.      ]])
    assert_allclose(t[0].as_matrix(), expected, atol=atol)

    # Test that the transformations are valid
    atol = 1e-12
    for m in t.as_matrix():
        assert_allclose(np.linalg.det(m[:3, :3]), 1, atol=atol)
        assert_allclose(np.linalg.norm(m[:3, 3]), 1, atol=atol)
        assert_allclose(m[3], [0, 0, 0, 1], atol=atol)


def test_apply():
    atol = 1e-12

    # Single transformation
    r = R.from_euler('z', 90, degrees=True)
    trans = np.array([2, 3, 4])
    t = T.from_transrot(trans, r)

    # Single vector, single transformation
    vec = np.array([1, 0, 0])
    expected = trans + r.apply(vec)
    res = t.apply(vec)
    assert_allclose(res, expected, atol=atol)

    # Multiple vectors, single transformation
    vecs = np.array([[1, 0, 0], [0, 1, 0]])
    expected = trans + r.apply(vecs)
    assert_allclose(t.apply(vecs), expected, atol=atol)

    # Multiple transformations
    r = R.from_euler('z', [90, 0], degrees=True)
    trans = np.array([[2, 3, 4], [5, 6, 7]])
    t = T.from_transrot(trans, r)

    # Single vector, multiple transformations
    vec = np.array([1, 0, 0])
    expected = trans + r.apply(vec)
    assert_allclose(t.apply(vec), expected, atol=atol)

    # Multiple vectors, multiple transformations
    vecs = np.array([[1, 0, 0], [0, 1, 0]])
    expected = trans + r.apply(vecs)
    assert_allclose(t.apply(vecs), expected, atol=atol)


def test_inverse_apply():
    atol = 1e-12

    # Test applying inverse transformation
    r = R.from_euler('z', 90, degrees=True)
    trans = np.array([1, 2, 3])
    t = T.from_transrot(trans, r)

    # Test single vector
    vec = np.array([1, 0, 0])
    expected = t.inv().apply(vec)
    assert_allclose(t.apply(vec, inverse=True), expected, atol=atol)

    # Test multiple vectors
    vecs = np.array([[1, 0, 0], [0, 1, 0]])
    expected = t.inv().apply(vecs)
    assert_allclose(t.apply(vecs, inverse=True), expected, atol=atol)


def test_rotation_alone():
    atol = 1e-12

    r = R.from_euler('z', 90, degrees=True)
    t = T.from_rotation(r)
    vec = np.array([1, 0, 0])
    expected = r.apply(vec)
    assert_allclose(t.apply(vec), expected, atol=atol)


def test_translation_alone():
    atol = 1e-12
    trans = np.array([1, 2, 3])
    t = T.from_translation(trans)
    vec = np.array([5, 6, 7])
    expected = trans + vec
    assert_allclose(t.apply(vec), expected, atol=atol)


def test_composition():
    atol = 1e-12

    # Test composing single transformations
    r1 = R.from_euler('z', 90, degrees=True)
    t1 = np.array([1, 0, 0])
    tf1 = T.from_transrot(t1, r1)

    r2 = R.from_euler('x', 90, degrees=True)
    t2 = np.array([0, 1, 0])
    tf2 = T.from_transrot(t2, r2)

    composed = tf2 * tf1
    vec = np.array([1, 0, 0])
    expected = tf2.apply(tf1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)
    assert composed.single

    # Multiple transformations with single transformation
    t2 = np.array([t2, t2])
    tf2 = T.from_transrot(t2, r2)

    composed = tf2 * tf1
    expected = tf2.apply(tf1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)
    assert not composed.single

    composed = tf1 * tf2
    expected = tf1.apply(tf2.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)
    assert not composed.single

    # Multiple transformations with multiple transformations
    t1 = np.array([t1, t1])
    tf1 = T.from_transrot(t1, r1)

    composed = tf2 * tf1
    expected = tf2.apply(tf1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)
    assert not composed.single


def test_pow():
    atol = 1e-12
    p = T.from_rotation(R.random(10, rng=0))
    p_inv = p.inv()

    # Test the short-cuts and other integers
    for n in [-5, -2, -1, 0, 1, 2, 5]:
        q = p**n
        r = T.identity(10)
        for _ in range(abs(n)):
            if n > 0:
                r = r * p
            else:
                r = r * p_inv
        assert_allclose(q.as_matrix(), r.as_matrix(), atol=atol)

        # Test shape preservation
        r = T.from_rotation(R.from_quat([0, 0, 0, 1]))
        assert (r**n).as_matrix().shape == (4, 4)
        r = T.from_rotation(R.from_quat([[0, 0, 0, 1]]))
        assert (r**n).as_matrix().shape == (1, 4, 4)


def test_inverse():
    atol = 1e-12

    # Test inverse transformation
    r = R.from_euler('z', 90, degrees=True)
    trans = np.array([1, 2, 3])
    t = T.from_transrot(trans, r)

    # Test that t * t.inv() equals identity
    t_inv = t.inv()
    composed = t * t_inv
    assert_allclose(composed.as_matrix(), np.eye(4), atol=atol)

    # Test multiple transformations
    r = R.from_euler('zyx', [[90, 0, 0], [0, 90, 0]], degrees=True)
    trans = np.array([[1, 2, 3], [4, 5, 6]])
    t = T.from_transrot(trans, r)
    t_inv = t.inv()
    composed = t * t_inv
    assert_allclose(composed.as_matrix(), np.array([np.eye(4)] * 2), atol=atol)


def test_properties():
    atol = 1e-12

    # Test rotation and translation properties
    r = R.from_euler('z', 90, degrees=True)
    trans = np.array([1, 2, 3])
    t = T.from_transrot(trans, r)

    assert_allclose(t.rotation.as_matrix(), r.as_matrix(), atol=atol)
    assert t.rotation.approx_equal(r)
    assert_allclose(t.translation, trans, atol=atol)


def test_indexing():
    atol = 1e-12

    # Test indexing for multiple transformations
    r = R.from_euler('zyx', [[90, 0, 0], [0, 90, 0]], degrees=True)
    trans = np.array([[1, 2, 3], [4, 5, 6]])
    t = T.from_transrot(trans, r)

    # Test single index
    assert_allclose(t[0].as_matrix()[:3, :3], r[0].as_matrix(), atol=atol)
    assert_allclose(t[0].as_matrix()[:3, 3], trans[0], atol=atol)

    # Test slice
    t_slice = t[0:2]
    assert_allclose(t_slice.as_matrix()[:, :3, :3], r[0:2].as_matrix(), atol=atol)
    assert_allclose(t_slice.as_matrix()[:, :3, 3], trans[0:2], atol=atol)


def test_concatenate():
    atol = 1e-12

    # Test concatenation of transformations
    r1 = R.from_euler('z', 90, degrees=True)
    t1 = np.array([1, 0, 0])
    tf1 = T.from_transrot(t1, r1)

    r2 = R.from_euler('x', 90, degrees=True)
    t2 = np.array([0, 1, 0])
    tf2 = T.from_transrot(t2, r2)

    # Concatenate single transformations
    concatenated = T.concatenate([tf1, tf2])
    assert_allclose(concatenated[0].as_matrix(), tf1.as_matrix(), atol=atol)
    assert_allclose(concatenated[1].as_matrix(), tf2.as_matrix(), atol=atol)

    # Concatenate multiple transformations
    concatenated = T.concatenate([tf1, T.concatenate([tf2, tf1])])
    assert_allclose(concatenated[0].as_matrix(), tf1.as_matrix(), atol=atol)
    assert_allclose(concatenated[1].as_matrix(), tf2.as_matrix(), atol=atol)
    assert_allclose(concatenated[2].as_matrix(), tf1.as_matrix(), atol=atol)


def test_input_validation():
    # Test invalid matrix shapes
    with pytest.raises(ValueError, match="Expected `matrix` to have shape"):
        T.from_matrix(np.eye(3))

    with pytest.raises(ValueError, match="Expected `matrix` to have shape"):
        T.from_matrix(np.zeros((4, 3)))

    with pytest.raises(ValueError, match="Expected `matrix` to have shape"):
        T.from_matrix([])

    with pytest.raises(ValueError, match="Expected `matrix` to have shape"):
        T.from_matrix(np.zeros((0, 4, 4)))

    with pytest.raises(ValueError, match="Expected `matrix` to have shape"):
        T.from_matrix(np.zeros((1, 1, 4, 4)))

    # Test invalid last row
    with pytest.raises(ValueError, match="Expected last row.*to be"):
        matrix = np.eye(4)
        matrix[3, :] = [1, 0, 0, 1]
        T.from_matrix(matrix)

    # Test non-rotation matrix
    with pytest.raises(ValueError,
                       match="Expected rotation matrix to have determinant 1"):
        matrix = np.eye(4)
        matrix[:3, :3] *= 2
        T.from_matrix(matrix)


def test_translation_validation():
    # Test invalid translation shapes
    with pytest.raises(ValueError, match="Expected `translation` to have shape"):
        T.from_translation([1, 2])

    with pytest.raises(ValueError, match="Expected `translation` to have shape"):
        T.from_translation(np.zeros((2, 2)))

    with pytest.raises(ValueError, match="Expected `translation` to have shape"):
        T.from_translation(np.zeros((0, 3)))

    with pytest.raises(ValueError, match="Expected `translation` to have shape"):
        T.from_translation(np.zeros((1, 1, 3)))


def test_vector_validation():
    t = T.identity()

    # Test invalid vector shapes
    with pytest.raises(ValueError, match="Expected vector to have shape"):
        t.apply([1, 2])

    with pytest.raises(ValueError, match="Expected vector to have shape"):
        t.apply(np.zeros((2, 2)))

    with pytest.raises(ValueError, match="Expected vector to have shape"):
        t.apply(np.zeros((0, 3)))

    with pytest.raises(ValueError, match="Expected vector to have shape"):
        t.apply(np.zeros((1, 1, 3)))


def test_indexing_validation():
    t = T.identity()

    # Test indexing on single transformation
    with pytest.raises(TypeError, match="Single transformation is not subscriptable"):
        t[0]

    with pytest.raises(TypeError, match="Single transformation is not subscriptable"):
        t[0:1]

    # Test length on single transformation
    with pytest.raises(TypeError, match="Single transformation has no len"):
        len(t)


def test_composition_validation():
    t2 = T.from_translation([[1, 2, 3], [4, 5, 6]])
    t3 = T.from_translation([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    # Test incompatible shapes
    with pytest.raises(ValueError, match="Expected equal number of transformations"):
        t2 * t3


def test_concatenate_validation():
    t = T.identity()

    # Test invalid inputs
    with pytest.raises(TypeError,
                       match="input must contain ProperRigidTransformation objects"):
        T.concatenate([t, np.eye(4)])


def test_setitem_validation():
    t = T.from_translation([[1, 2, 3], [4, 5, 6]])
    single = T.identity()

    # Test setting item on single transformation
    with pytest.raises(TypeError, match="Single transformation is not subscriptable"):
        single[0] = t

    # Test invalid value type
    with pytest.raises(TypeError, match="value must be a ProperRigidTransformation"):
        t[0] = np.eye(4)


def test_copy_flag():
    # Test that copy=True creates new memory
    matrix = np.eye(4)
    t = T.from_matrix(matrix, copy=True)
    matrix[0, 0] = 2
    assert t.as_matrix()[0, 0] == 1

    # Test that copy=False shares memory
    matrix = np.eye(4)
    t = T.from_matrix(matrix, copy=False)
    matrix[0, 0] = 2
    assert t.as_matrix()[0, 0] == 2
