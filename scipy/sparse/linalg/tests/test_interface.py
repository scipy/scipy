"""Test functions for the sparse.linalg.interface module
"""

from __future__ import division, print_function, absolute_import

from numpy.testing import (TestCase, assert_, assert_equal,
        assert_raises, assert_allclose)

import numpy as np
import scipy.sparse as sparse

from scipy.sparse.linalg import interface


class TestLinearOperator(TestCase):
    def setUp(self):
        self.matvecs = []

        # these matvecs do not preserve type or shape
        def matvec1(x):
            return np.array([1*x[0] + 2*x[1] + 3*x[2],
                              4*x[0] + 5*x[1] + 6*x[2]])

        def matvec2(x):
            return np.matrix(matvec1(x).reshape(2,1))

        self.matvecs.append(matvec1)
        self.matvecs.append(matvec2)

    def test_matvec(self):

        for matvec in self.matvecs:
            A = interface.LinearOperator((2,3), matvec)

            assert_equal(A.matvec(np.array([1,2,3])), [14,32])
            assert_equal(A.matvec(np.array([[1],[2],[3]])), [[14],[32]])
            assert_equal(A * np.array([1,2,3]), [14,32])
            assert_equal(A * np.array([[1],[2],[3]]), [[14],[32]])
            assert_equal(A.dot(np.array([1,2,3])), [14,32])
            assert_equal(A.dot(np.array([[1],[2],[3]])), [[14],[32]])

            assert_equal(A.matvec(np.matrix([[1],[2],[3]])), [[14],[32]])
            assert_equal(A * np.matrix([[1],[2],[3]]), [[14],[32]])
            assert_equal(A.dot(np.matrix([[1],[2],[3]])), [[14],[32]])

            assert_(isinstance(A.matvec(np.array([1,2,3])), np.ndarray))
            assert_(isinstance(A.matvec(np.array([[1],[2],[3]])), np.ndarray))
            assert_(isinstance(A * np.array([1,2,3]), np.ndarray))
            assert_(isinstance(A * np.array([[1],[2],[3]]), np.ndarray))
            assert_(isinstance(A.dot(np.array([1,2,3])), np.ndarray))
            assert_(isinstance(A.dot(np.array([[1],[2],[3]])), np.ndarray))

            assert_(isinstance(A.matvec(np.matrix([[1],[2],[3]])), np.ndarray))
            assert_(isinstance(A * np.matrix([[1],[2],[3]]), np.ndarray))
            assert_(isinstance(A.dot(np.matrix([[1],[2],[3]])), np.ndarray))

            assert_raises(ValueError, A.matvec, np.array([1,2]))
            assert_raises(ValueError, A.matvec, np.array([1,2,3,4]))
            assert_raises(ValueError, A.matvec, np.array([[1],[2]]))
            assert_raises(ValueError, A.matvec, np.array([[1],[2],[3],[4]]))


class TestAsLinearOperator(TestCase):
    def setUp(self):
        self.cases = []

        def make_cases(dtype):
            self.cases.append(np.matrix([[1,2,3],[4,5,6]], dtype=dtype))
            self.cases.append(np.array([[1,2,3],[4,5,6]], dtype=dtype))
            self.cases.append(sparse.csr_matrix([[1,2,3],[4,5,6]], dtype=dtype))

            class matlike:
                def __init__(self, dtype):
                    self.dtype = np.dtype(dtype)
                    self.shape = (2,3)

                def matvec(self,x):
                    y = np.array([1*x[0] + 2*x[1] + 3*x[2],
                                   4*x[0] + 5*x[1] + 6*x[2]], dtype=self.dtype)
                    if len(x.shape) == 2:
                        y = y.reshape(-1,1)
                    return y

                def rmatvec(self,x):
                    return np.array([1*x[0] + 4*x[1],
                                      2*x[0] + 5*x[1],
                                      3*x[0] + 6*x[1]], dtype=self.dtype)
            self.cases.append(matlike('int'))

        make_cases('int32')
        make_cases('float32')
        make_cases('float64')

    def test_basic(self):

        for M in self.cases:
            A = interface.aslinearoperator(M)
            M,N = A.shape

            assert_equal(A.matvec(np.array([1,2,3])), [14,32])
            assert_equal(A.matvec(np.array([[1],[2],[3]])), [[14],[32]])

            assert_equal(A * np.array([1,2,3]), [14,32])
            assert_equal(A * np.array([[1],[2],[3]]), [[14],[32]])

            assert_equal(A.rmatvec(np.array([1,2])), [9,12,15])
            assert_equal(A.rmatvec(np.array([[1],[2]])), [[9],[12],[15]])

            assert_equal(
                    A.matmat(np.array([[1,4],[2,5],[3,6]])),
                    [[14,32],[32,77]])

            assert_equal(A * np.array([[1,4],[2,5],[3,6]]), [[14,32],[32,77]])

            if hasattr(M,'dtype'):
                assert_equal(A.dtype, M.dtype)

    def test_dot(self):

        for M in self.cases:
            A = interface.aslinearoperator(M)
            M,N = A.shape

            assert_equal(A.dot(np.array([1,2,3])), [14,32])
            assert_equal(A.dot(np.array([[1],[2],[3]])), [[14],[32]])

            assert_equal(
                    A.dot(np.array([[1,4],[2,5],[3,6]])),
                    [[14,32],[32,77]])


class TestProductAndPowerOperators(TestCase):

    def test_matrix_product_operator(self):
        np.random.seed(1234)
        na, nb, nc = 14, 20, 16
        k = 2
        A = np.random.randn(na, nb)
        B = np.random.randn(nb, nc)
        C = np.random.randn(nc, na)
        D = np.random.randn(na, k)
        v = np.random.randn(na)
        M_explicit = A.dot(B).dot(C)
        M_implicit = interface.MatrixProductOperator(A, B, C)
        M_implicit_T = interface.MatrixProductOperator(C.T, B.T, A.T)
        M_op = interface._ProductOp(A, B, C)
        for target in (D, v):
            assert_allclose(M_explicit.dot(target), M_implicit.dot(target))
            assert_allclose(M_implicit.dot(target), M_op.dot(target))
            assert_allclose(M_implicit_T.dot(target), M_op.T.dot(target))
        assert_allclose(M_explicit.dot(D), M_implicit.matmat(D))
        assert_allclose(M_explicit.dot(v), M_implicit.matvec(v))
        assert_allclose(M_explicit.T.dot(D), M_implicit_T.matmat(D))
        assert_allclose(M_explicit.T.dot(v), M_implicit_T.matvec(v))

        # check compatibility with the matrix linear operator
        matrix_operator = interface.MatrixLinearOperator(M_explicit)
        assert_allclose(matrix_operator.matmat(D), M_implicit.matmat(D))
        assert_allclose(matrix_operator.matvec(v), M_implicit.matvec(v))

    def test_matrix_power_operator(self):
        np.random.seed(1234)
        n = 4
        k = 2
        A = np.random.randn(n, n)
        B = np.random.randn(n, k)
        v = np.random.randn(n)
        for p in range(1, 4):
            M_explicit = np.linalg.matrix_power(A, p)
            M_implicit = interface.MatrixPowerOperator(A, p)
            M_implicit_T = interface.MatrixPowerOperator(A.T, p)
            M_op = interface._PowerOp(A, p)
            for target in (B, v):
                assert_allclose(M_explicit.dot(target), M_implicit.dot(target))
                assert_allclose(M_implicit.dot(target), M_op.dot(target))
                assert_allclose(M_implicit_T.dot(target), M_op.T.dot(target))
            assert_allclose(M_explicit.dot(B), M_implicit.matmat(B))
            assert_allclose(M_explicit.dot(v), M_implicit.matvec(v))
            assert_allclose(M_explicit.T.dot(B), M_implicit_T.matmat(B))
            assert_allclose(M_explicit.T.dot(v), M_implicit_T.matvec(v))

            # check compatibility with the matrix linear operator
            matrix_operator = interface.MatrixLinearOperator(M_explicit)
            assert_allclose(matrix_operator.matmat(B), M_implicit.matmat(B))
            assert_allclose(matrix_operator.matvec(v), M_implicit.matvec(v))


