import pytest

import numpy as np
from numpy.testing import assert_allclose
from scipy.spatial.transform import Rotation, ProperRigidTransformation


def test_from_rotation():
    atol = 1e-12

    # Test single rotation
    r = Rotation.identity()
    t = ProperRigidTransformation.from_rotation(r)
    assert_allclose(t.as_matrix(), np.eye(4), atol=atol)
    assert t.single

    r = Rotation.from_euler('z', 90, degrees=True)
    t = ProperRigidTransformation.from_rotation(r)
    assert_allclose(t.as_matrix()[:3, :3], r.as_matrix(), atol=atol)
    assert_allclose(t.as_matrix()[:3, 3], [0, 0, 0], atol=atol)
    assert_allclose(t.as_matrix()[3], [0, 0, 0, 1], atol=atol)
    assert t.single

    # Test multiple rotations
    r = Rotation.from_euler('zyx', [[90, 0, 0], [0, 90, 0]], degrees=True)
    t = ProperRigidTransformation.from_rotation(r)
    assert_allclose(t.as_matrix()[:, :3, :3], r.as_matrix(), atol=atol)
    assert_allclose(t.as_matrix()[:, :3, 3], [[0, 0, 0], [0, 0, 0]], atol=atol)
    assert_allclose(t.as_matrix()[:, 3], [[0, 0, 0, 1], [0, 0, 0, 1]], atol=atol)
    assert not t.single


def test_from_translation():
    # Test single translation
    trans = np.array([1, 2, 3])
    t = ProperRigidTransformation.from_translation(trans)
    expected = np.eye(4)
    expected[:3, 3] = trans
    assert_allclose(t.as_matrix(), expected)
    assert t.single

    # Test multiple translations
    trans = np.array([[1, 2, 3], [4, 5, 6]])
    t = ProperRigidTransformation.from_translation(trans)
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
    t = ProperRigidTransformation.from_matrix(matrix)
    assert_allclose(t.as_matrix(), matrix, atol=atol)
    assert t.single

    # Test multiple transformation matrices
    matrices = np.array([np.eye(4)]*2)
    matrices[0, :3, 3] = [1, 2, 3]
    matrices[1, :3, 3] = [4, 5, 6]
    t = ProperRigidTransformation.from_matrix(matrices)
    assert_allclose(t.as_matrix(), matrices, atol=atol)
    assert not t.single

    # Test invalid matrix
    with pytest.raises(ValueError):
        invalid = np.eye(4)
        invalid[3, 3] = 2  # Invalid last row
        ProperRigidTransformation.from_matrix(invalid)


def test_from_transrot():
    atol = 1e-12

    # Test single rotation and translation
    r = Rotation.from_euler('zyx', [90, 0, 0], degrees=True)
    trans = np.array([1, 2, 3])
    t = ProperRigidTransformation.from_transrot(trans, r)

    expected = np.zeros((4, 4))
    expected[:3, :3] = r.as_matrix()
    expected[:3, 3] = trans
    expected[3, 3] = 1
    assert_allclose(t.as_matrix(), expected, atol=atol)
    assert t.single

    # Test single rotation and multiple translations
    r = Rotation.from_euler('z', 90, degrees=True)
    trans = np.array([[1, 2, 3], [4, 5, 6]])
    t = ProperRigidTransformation.from_transrot(trans, r)
    assert not t.single

    for i in range(len(trans)):
        expected = np.zeros((4, 4))
        expected[:3, :3] = r.as_matrix()
        expected[:3, 3] = trans[i]
        expected[3, 3] = 1
        assert_allclose(t.as_matrix()[i], expected, atol=atol)

    # Test multiple rotations and translations
    r = Rotation.from_euler('zyx', [[90, 0, 0], [0, 90, 0]], degrees=True)
    trans = np.array([[1, 2, 3], [4, 5, 6]])
    t = ProperRigidTransformation.from_transrot(trans, r)
    assert not t.single

    for i in range(len(trans)):
        expected = np.zeros((4, 4))
        expected[:3, :3] = r[i].as_matrix()
        expected[:3, 3] = trans[i]
        expected[3, 3] = 1
        assert_allclose(t.as_matrix()[i], expected, atol=atol)


def test_identity():
    atol = 1e-12

    # Test single identity
    t = ProperRigidTransformation.identity()
    assert_allclose(t.as_matrix(), np.eye(4), atol=atol)

    # Test multiple identities
    t = ProperRigidTransformation.identity(5)
    assert_allclose(t.as_matrix(), np.array([np.eye(4)] * 5), atol=atol)


def test_random():
    atol = 1e-6
    n = 5
    t = ProperRigidTransformation.random(n, rng=0)
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

    ## Single transformation
    r = Rotation.from_euler('z', 90, degrees=True)
    trans = np.array([2, 3, 4])
    t = ProperRigidTransformation.from_transrot(trans, r)

    # Single vector, single transformation
    vec = np.array([1, 0, 0])
    expected = trans + r.apply(vec)
    res = t.apply(vec)
    assert_allclose(res, expected, atol=atol)

    # Multiple vectors, single transformation
    vecs = np.array([[1, 0, 0], [0, 1, 0]])
    expected = trans + r.apply(vecs)
    assert_allclose(t.apply(vecs), expected, atol=atol)

    ## Multiple transformations
    r = Rotation.from_euler('z', [90, 0], degrees=True)
    trans = np.array([[2, 3, 4], [5, 6, 7]])
    t = ProperRigidTransformation.from_transrot(trans, r)

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
    r = Rotation.from_euler('z', 90, degrees=True)
    trans = np.array([1, 2, 3])
    t = ProperRigidTransformation.from_transrot(trans, r)

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

    r = Rotation.from_euler('z', 90, degrees=True)
    t = ProperRigidTransformation.from_rotation(r)
    vec = np.array([1, 0, 0])
    expected = r.apply(vec)
    assert_allclose(t.apply(vec), expected, atol=atol)


def test_translation_alone():
    atol = 1e-12
    trans = np.array([1, 2, 3])
    t = ProperRigidTransformation.from_translation(trans)
    vec = np.array([5, 6, 7])
    expected = trans + vec
    assert_allclose(t.apply(vec), expected, atol=atol)


def test_composition():
    atol = 1e-12

    # Test composing single transformations
    r1 = Rotation.from_euler('z', 90, degrees=True)
    t1 = np.array([1, 0, 0])
    tf1 = ProperRigidTransformation.from_transrot(t1, r1)

    r2 = Rotation.from_euler('x', 90, degrees=True)
    t2 = np.array([0, 1, 0])
    tf2 = ProperRigidTransformation.from_transrot(t2, r2)

    composed = tf2 * tf1
    vec = np.array([1, 0, 0])
    expected = tf2.apply(tf1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)
    assert composed.single

    expected = t2 + r2.apply(t1 + r1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)

    # Multiple transformations with single transformation
    t2 = np.array([[1, 2, 3], [4, 5, 6]])
    tf2 = ProperRigidTransformation.from_transrot(t2, r2)

    composed = tf2 * tf1
    expected = tf2.apply(tf1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)
    assert not composed.single

    expected = t2 + r2.apply(t1 + r1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)

    # Multiple transformations with multiple transformations
    t1 = np.array([[1, 0, 0], [0, -1, 1]])
    tf1 = ProperRigidTransformation.from_transrot(t1, r1)

    composed = tf2 * tf1
    expected = tf2.apply(tf1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)
    assert not composed.single

    expected = t2 + r2.apply(t1 + r1.apply(vec))
    assert_allclose(composed.apply(vec), expected, atol=atol)


def test_pow():
    atol = 1e-12
    num = 10
    p = ProperRigidTransformation.random(num, rng=0)
    p_inv = p.inv()

    # Test the short-cuts and other integers
    for n in [-5, -2, -1, 0, 1, 2, 5]:
        q = p**n
        r = ProperRigidTransformation.identity(num)
        for _ in range(abs(n)):
            if n > 0:
                r = r * p
            else:
                r = r * p_inv
        assert_allclose(q.as_matrix(), r.as_matrix(), atol=atol)

        # Test shape preservation
        r = ProperRigidTransformation.from_rotation(Rotation.from_quat([0, 0, 0, 1]))
        assert (r**n).as_matrix().shape == (4, 4)
        r = ProperRigidTransformation.from_rotation(Rotation.from_quat([[0, 0, 0, 1]]))
        assert (r**n).as_matrix().shape == (1, 4, 4)


def test_inverse():
    atol = 1e-12

    # Test inverse transformation
    r = Rotation.from_euler('z', 90, degrees=True)
    trans = np.array([1, 2, 3])
    t = ProperRigidTransformation.from_transrot(trans, r)

    # Test that t * t.inv() equals identity
    t_inv = t.inv()
    composed = t * t_inv
    assert_allclose(composed.as_matrix(), np.eye(4), atol=atol)

    n = 10
    t = ProperRigidTransformation.random(n, rng=1)
    t_inv = t.inv()
    composed = t * t_inv
    assert_allclose(composed.as_matrix(), np.array([np.eye(4)] * n), atol=atol)

    # Test multiple transformations
    r = Rotation.from_euler('zyx', [[90, 0, 0], [0, 90, 0]], degrees=True)
    trans = np.array([[1, 2, 3], [4, 5, 6]])
    t = ProperRigidTransformation.from_transrot(trans, r)
    t_inv = t.inv()
    composed = t * t_inv
    assert_allclose(composed.as_matrix(), np.array([np.eye(4)] * 2), atol=atol)


def test_properties():
    atol = 1e-12

    # Test rotation and translation properties
    r = Rotation.from_euler('z', 90, degrees=True)
    trans = np.array([1, 2, 3])
    t = ProperRigidTransformation.from_transrot(trans, r)

    assert_allclose(t.rotation.as_matrix(), r.as_matrix(), atol=atol)
    assert t.rotation.approx_equal(r)
    assert_allclose(t.translation, trans, atol=atol)


def test_indexing():
    atol = 1e-12

    # Test indexing for multiple transformations
    r = Rotation.from_euler('zyx', [[90, 0, 0], [0, 90, 0]], degrees=True)
    trans = np.array([[1, 2, 3], [4, 5, 6]])
    t = ProperRigidTransformation.from_transrot(trans, r)

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
    r1 = Rotation.from_euler('z', 90, degrees=True)
    t1 = np.array([1, 0, 0])
    tf1 = ProperRigidTransformation.from_transrot(t1, r1)

    r2 = Rotation.from_euler('x', 90, degrees=True)
    t2 = np.array([0, 1, 0])
    tf2 = ProperRigidTransformation.from_transrot(t2, r2)

    # Concatenate single transformations
    concatenated1 = ProperRigidTransformation.concatenate([tf1, tf2])
    assert_allclose(concatenated1[0].as_matrix(), tf1.as_matrix(), atol=atol)
    assert_allclose(concatenated1[1].as_matrix(), tf2.as_matrix(), atol=atol)

    # Concatenate multiple transformations
    concatenated2 = ProperRigidTransformation.concatenate([tf1, concatenated1])
    assert_allclose(concatenated2[0].as_matrix(), tf1.as_matrix(), atol=atol)
    assert_allclose(concatenated2[1].as_matrix(), tf1.as_matrix(), atol=atol)
    assert_allclose(concatenated2[2].as_matrix(), tf2.as_matrix(), atol=atol)


def test_input_validation():
    # Test invalid matrix shapes
    with pytest.raises(ValueError, match="Expected `matrix` to have shape"):
        ProperRigidTransformation.from_matrix(np.eye(3))

    with pytest.raises(ValueError, match="Expected `matrix` to have shape"):
        ProperRigidTransformation.from_matrix(np.zeros((4, 3)))

    with pytest.raises(ValueError, match="Expected `matrix` to have shape"):
        ProperRigidTransformation.from_matrix([])

    with pytest.raises(ValueError, match="Expected `matrix` to have shape"):
        ProperRigidTransformation.from_matrix(np.zeros((0, 4, 4)))

    with pytest.raises(ValueError, match="Expected `matrix` to have shape"):
        ProperRigidTransformation.from_matrix(np.zeros((1, 1, 4, 4)))

    # Test invalid last row
    with pytest.raises(ValueError, match="Expected last row.*to be"):
        matrix = np.eye(4)
        matrix[3, :] = [1, 0, 0, 1]
        ProperRigidTransformation.from_matrix(matrix)

    # Test non-rotation matrix
    with pytest.raises(ValueError,
                       match="Expected rotation matrix to have determinant 1"):
        matrix = np.eye(4)
        matrix[:3, :3] *= 2
        ProperRigidTransformation.from_matrix(matrix)


def test_translation_validation():
    # Test invalid translation shapes
    with pytest.raises(ValueError, match="Expected `translation` to have shape"):
        ProperRigidTransformation.from_translation([1, 2])

    with pytest.raises(ValueError, match="Expected `translation` to have shape"):
        ProperRigidTransformation.from_translation(np.zeros((2, 2)))

    with pytest.raises(ValueError, match="Expected `translation` to have shape"):
        ProperRigidTransformation.from_translation(np.zeros((0, 3)))

    with pytest.raises(ValueError, match="Expected `translation` to have shape"):
        ProperRigidTransformation.from_translation(np.zeros((1, 1, 3)))


def test_vector_validation():
    t = ProperRigidTransformation.identity()

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
    t = ProperRigidTransformation.identity()

    # Test indexing on single transformation
    with pytest.raises(TypeError, match="Single transformation is not subscriptable"):
        t[0]

    with pytest.raises(TypeError, match="Single transformation is not subscriptable"):
        t[0:1]

    # Test length on single transformation
    with pytest.raises(TypeError, match="Single transformation has no len"):
        len(t)


def test_composition_validation():
    t2 = ProperRigidTransformation.from_translation([[1, 2, 3], [4, 5, 6]])
    t3 = ProperRigidTransformation.from_translation([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    # Test incompatible shapes
    with pytest.raises(ValueError, match="Expected equal number of transformations"):
        t2 * t3


def test_concatenate_validation():
    t = ProperRigidTransformation.identity()

    # Test invalid inputs
    with pytest.raises(TypeError,
                       match="input must contain ProperRigidTransformation objects"):
        ProperRigidTransformation.concatenate([t, np.eye(4)])


def test_setitem_validation():
    t = ProperRigidTransformation.from_translation([[1, 2, 3], [4, 5, 6]])
    single = ProperRigidTransformation.identity()

    # Test setting item on single transformation
    with pytest.raises(TypeError, match="Single transformation is not subscriptable"):
        single[0] = t

    # Test invalid value type
    with pytest.raises(TypeError, match="value must be a ProperRigidTransformation"):
        t[0] = np.eye(4)


def test_copy_flag():
    # Test that copy=True creates new memory
    matrix = np.eye(4)
    t = ProperRigidTransformation.from_matrix(matrix, copy=True)
    matrix[0, 0] = 2
    assert t.as_matrix()[0, 0] == 1

    # Test that copy=False shares memory
    matrix = np.eye(4)
    t = ProperRigidTransformation.from_matrix(matrix, copy=False)
    matrix[0, 0] = 2
    assert t.as_matrix()[0, 0] == 2
