"""Tests for functions in special_matrices.py."""

from numpy import arange, add, array, eye, all, copy
from numpy.testing import *

from scipy.linalg import toeplitz, hankel, circulant, hadamard, leslie, \
                            companion, tri, triu, tril, kron, block_diag


def get_mat(n):
    data = arange(n)
    data = add.outer(data,data)
    return data


class TestTri(TestCase):
    def test_basic(self):
        assert_equal(tri(4),array([[1,0,0,0],
                                   [1,1,0,0],
                                   [1,1,1,0],
                                   [1,1,1,1]]))
        assert_equal(tri(4,dtype='f'),array([[1,0,0,0],
                                                [1,1,0,0],
                                                [1,1,1,0],
                                                [1,1,1,1]],'f'))
    def test_diag(self):
        assert_equal(tri(4,k=1),array([[1,1,0,0],
                                       [1,1,1,0],
                                       [1,1,1,1],
                                       [1,1,1,1]]))
        assert_equal(tri(4,k=-1),array([[0,0,0,0],
                                        [1,0,0,0],
                                        [1,1,0,0],
                                        [1,1,1,0]]))
    def test_2d(self):
        assert_equal(tri(4,3),array([[1,0,0],
                                     [1,1,0],
                                     [1,1,1],
                                     [1,1,1]]))
        assert_equal(tri(3,4),array([[1,0,0,0],
                                     [1,1,0,0],
                                     [1,1,1,0]]))
    def test_diag2d(self):
        assert_equal(tri(3,4,k=2),array([[1,1,1,0],
                                         [1,1,1,1],
                                         [1,1,1,1]]))
        assert_equal(tri(4,3,k=-2),array([[0,0,0],
                                          [0,0,0],
                                          [1,0,0],
                                          [1,1,0]]))

class TestTril(TestCase):
    def test_basic(self):
        a = (100*get_mat(5)).astype('l')
        b = a.copy()
        for k in range(5):
            for l in range(k+1,5):
                b[k,l] = 0
        assert_equal(tril(a),b)

    def test_diag(self):
        a = (100*get_mat(5)).astype('f')
        b = a.copy()
        for k in range(5):
            for l in range(k+3,5):
                b[k,l] = 0
        assert_equal(tril(a,k=2),b)
        b = a.copy()
        for k in range(5):
            for l in range(max((k-1,0)),5):
                b[k,l] = 0
        assert_equal(tril(a,k=-2),b)


class TestTriu(TestCase):
    def test_basic(self):
        a = (100*get_mat(5)).astype('l')
        b = a.copy()
        for k in range(5):
            for l in range(k+1,5):
                b[l,k] = 0
        assert_equal(triu(a),b)

    def test_diag(self):
        a = (100*get_mat(5)).astype('f')
        b = a.copy()
        for k in range(5):
            for l in range(max((k-1,0)),5):
                b[l,k] = 0
        assert_equal(triu(a,k=2),b)
        b = a.copy()
        for k in range(5):
            for l in range(k+3,5):
                b[l,k] = 0
        assert_equal(triu(a,k=-2),b)


class TestToeplitz(TestCase):
    
    def test_basic(self):
        y = toeplitz([1,2,3])
        assert_array_equal(y,[[1,2,3],[2,1,2],[3,2,1]])
        y = toeplitz([1,2,3],[1,4,5])
        assert_array_equal(y,[[1,4,5],[2,1,4],[3,2,1]])
        
    def test_complex_01(self):
        data = (1.0 + arange(3.0)) * (1.0 + 1.0j)
        x = copy(data)
        t = toeplitz(x)
        # Calling toeplitz should not change x.
        assert_array_equal(x, data)
        # According to the docstring, x should be the first column of t.
        col0 = t[:,0]
        assert_array_equal(col0, data)
        assert_array_equal(t[0,1:], data[1:].conj())

    def test_scalar_00(self):
        """Scalar arguments still produce a 2D array."""
        t = toeplitz(10)
        assert_array_equal(t, [[10]])
        t = toeplitz(10, 20)
        assert_array_equal(t, [[10]])
        
    def test_scalar_01(self):
        c = array([1,2,3])
        t = toeplitz(c, 1)
        assert_array_equal(t, [[1],[2],[3]])

    def test_scalar_02(self):
        c = array([1,2,3])
        t = toeplitz(c, array(1))
        assert_array_equal(t, [[1],[2],[3]])

    def test_scalar_03(self):
        c = array([1,2,3])
        t = toeplitz(c, array([1]))
        assert_array_equal(t, [[1],[2],[3]])

    def test_scalar_04(self):
        r = array([10,2,3])
        t = toeplitz(1, r)
        assert_array_equal(t, [[1,2,3]])


class TestHankel(TestCase):
    def test_basic(self):
        y = hankel([1,2,3])
        assert_array_equal(y, [[1,2,3], [2,3,0], [3,0,0]])
        y = hankel([1,2,3], [3,4,5])
        assert_array_equal(y, [[1,2,3], [2,3,4], [3,4,5]])


class TestCirculant(TestCase):
    def test_basic(self):
        y = circulant([1,2,3])
        assert_array_equal(y, [[1,3,2], [2,1,3], [3,2,1]])


class TestHadamard(TestCase):

    def test_basic(self):

        y = hadamard(1)
        assert_array_equal(y, [[1]])

        y = hadamard(2, dtype=float)
        assert_array_equal(y, [[1.0, 1.0], [1.0, -1.0]])

        y = hadamard(4)        
        assert_array_equal(y, [[1,1,1,1], [1,-1,1,-1], [1,1,-1,-1], [1,-1,-1,1]])

        assert_raises(ValueError, hadamard, 0)
        assert_raises(ValueError, hadamard, 5)


class TestLeslie(TestCase):

    def test_bad_shapes(self):
        assert_raises(ValueError, leslie, [[1,1],[2,2]], [3,4,5])        
        assert_raises(ValueError, leslie, [3,4,5], [[1,1],[2,2]])
        assert_raises(ValueError, leslie, [1,2], [1,2])
        assert_raises(ValueError, leslie, [1], [])

    def test_basic(self):
        a = leslie([1, 2, 3], [0.25, 0.5])
        expected = array([
            [1.0,  2.0, 3.0],
            [0.25, 0.0, 0.0],
            [0.0,  0.5, 0.0]])
        assert_array_equal(a, expected)


class TestCompanion(TestCase):

    def test_bad_shapes(self):
        assert_raises(ValueError, companion, [[1,1],[2,2]])        
        assert_raises(ValueError, companion, [0,4,5])
        assert_raises(ValueError, companion, [1])
        assert_raises(ValueError, companion, [])

    def test_basic(self):
        c = companion([1, 2, 3])
        expected = array([
            [-2.0, -3.0],
            [ 1.0,  0.0]])
        assert_array_equal(c, expected)

        c = companion([2.0, 5.0, -10.0])
        expected = array([
            [-2.5, 5.0],
            [ 1.0, 0.0]])
        assert_array_equal(c, expected)


class TestBlockDiag:
    def test_basic(self):
        x = block_diag(eye(2), [[1,2], [3,4], [5,6]], [[1, 2, 3]])
        assert all(x == [[1, 0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0, 0],
                         [0, 0, 1, 2, 0, 0, 0],
                         [0, 0, 3, 4, 0, 0, 0],
                         [0, 0, 5, 6, 0, 0, 0],
                         [0, 0, 0, 0, 1, 2, 3]])

    def test_dtype(self):
        x = block_diag([[1.5]])
        assert_equal(x.dtype, float)

        x = block_diag([[True]])
        assert_equal(x.dtype, bool)
        
    def test_scalar_and_1d_args(self):
        a = block_diag(1)
        assert_equal(a.shape, (1,1))
        assert_array_equal(a, [[1]])
        
        a = block_diag([2,3], 4)
        assert_array_equal(a, [[2, 3, 0], [0, 0, 4]])

    def test_bad_arg(self):
        assert_raises(ValueError, block_diag, [[[1]]])

    def test_no_args(self):
        a = block_diag()
        assert_equal(a.ndim, 2)
        assert_equal(a.nbytes, 0)


class TestKron:

    def test_basic(self):

        a = kron(array([[1, 2], [3, 4]]), array([[1, 1, 1]]))
        assert_array_equal(a, array([[1, 1, 1, 2, 2, 2],
                                     [3, 3, 3, 4, 4, 4]]))

        m1 = array([[1, 2], [3, 4]])
        m2 = array([[10], [11]])
        a = kron(m1, m2)
        expected = array([[ 10, 20 ],
                          [ 11, 22 ],
                          [ 30, 40 ],
                          [ 33, 44 ]])
        assert_array_equal(a, expected)


if __name__ == "__main__":
    run_module_suite()
