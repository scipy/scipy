from __future__ import division, print_function, absolute_import

import warnings

from numpy import array, finfo, arange, eye, all, unique, ones, dot, matrix
import numpy.random as random
from numpy.testing import TestCase, run_module_suite, assert_array_almost_equal, \
    assert_raises, assert_almost_equal, assert_equal, assert_array_equal, assert_

import scipy.linalg
from scipy.linalg import norm, inv
from scipy.sparse import spdiags, SparseEfficiencyWarning, csc_matrix, csr_matrix
from scipy.sparse.linalg.dsolve import spsolve, use_solver, splu, spilu

warnings.simplefilter('ignore',SparseEfficiencyWarning)

# TODO add more comprehensive tests
use_solver(useUmfpack=False)


class TestLinsolve(TestCase):
    def test_singular(self):
        A = csc_matrix((5,5), dtype='d')
        b = array([1, 2, 3, 4, 5],dtype='d')
        x = spsolve(A, b, use_umfpack=False)

    def test_twodiags(self):
        A = spdiags([[1, 2, 3, 4, 5], [6, 5, 8, 9, 10]], [0, 1], 5, 5)
        b = array([1, 2, 3, 4, 5])

        # condition number of A
        cond_A = norm(A.todense(),2) * norm(inv(A.todense()),2)

        for t in ['f','d','F','D']:
            eps = finfo(t).eps  # floating point epsilon
            b = b.astype(t)

            for format in ['csc','csr']:
                Asp = A.astype(t).asformat(format)

                x = spsolve(Asp,b)

                assert_(norm(b - Asp*x) < 10 * cond_A * eps)

    def test_bvector_smoketest(self):
        Adense = matrix([[0., 1., 1.],
                         [1., 0., 1.],
                         [0., 0., 1.]])
        As = csc_matrix(Adense)
        random.seed(1234)
        x = random.randn(3)
        b = As*x
        x2 = spsolve(As, b)

        assert_array_almost_equal(x, x2)

    def test_bmatrix_smoketest(self):
        Adense = matrix([[0., 1., 1.],
                         [1., 0., 1.],
                         [0., 0., 1.]])
        As = csc_matrix(Adense)
        random.seed(1234)
        x = random.randn(3, 4)
        Bdense = As.dot(x)
        Bs = csc_matrix(Bdense)
        x2 = spsolve(As, Bs)
        assert_array_almost_equal(x, x2.todense())

    def test_non_square(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=SparseEfficiencyWarning)
            # A is not square.
            A = ones((3, 4))
            b = ones((4, 1))
            assert_raises(ValueError, spsolve, A, b)
            # A2 and b2 have incompatible shapes.
            A2 = csc_matrix(eye(3))
            b2 = array([1.0, 2.0])
            assert_raises(ValueError, spsolve, A2, b2)

    def test_example_comparison(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=SparseEfficiencyWarning)
            row = array([0,0,1,2,2,2])
            col = array([0,2,2,0,1,2])
            data = array([1,2,3,-4,5,6])
            sM = csr_matrix((data,(row,col)), shape=(3,3), dtype=float)
            M = sM.todense()

            row = array([0,0,1,1,0,0])
            col = array([0,2,1,1,0,0])
            data = array([1,1,1,1,1,1])
            sN = csr_matrix((data, (row,col)), shape=(3,3), dtype=float)
            N = sN.todense()

            sX = spsolve(sM, sN)
            X = scipy.linalg.solve(M, N)

            assert_array_almost_equal(X, sX.todense())


class TestSplu(object):
    def setUp(self):
        n = 40
        d = arange(n) + 1
        self.n = n
        self.A = spdiags((d, 2*d, d[::-1]), (-3, 0, 5), n, n)
        random.seed(1234)

    def test_splu_smoketest(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=SparseEfficiencyWarning)
            # Check that splu works at all
            x = random.rand(self.n)
            lu = splu(self.A)
            r = self.A*lu.solve(x)
            assert_(abs(x - r).max() < 1e-13)

    def test_spilu_smoketest(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=SparseEfficiencyWarning)
            # Check that spilu works at all
            x = random.rand(self.n)
            lu = spilu(self.A, drop_tol=1e-2, fill_factor=5)
            r = self.A*lu.solve(x)
            assert_(abs(x - r).max() < 1e-2)
            assert_(abs(x - r).max() > 1e-5)

    def test_splu_nnz0(self):
        A = csc_matrix((5,5), dtype='d')
        assert_raises(RuntimeError, splu, A)

    def test_spilu_nnz0(self):
        A = csc_matrix((5,5), dtype='d')
        assert_raises(RuntimeError, spilu, A)

    def test_splu_basic(self):
        # Test basic splu functionality.
        n = 30
        a = random.random((n, n))
        a[a < 0.95] = 0
        # First test with a singular matrix
        a[:, 0] = 0
        a_ = csc_matrix(a)
        # Matrix is exactly singular
        assert_raises(RuntimeError, splu, a_)

        # Make a diagonal dominant, to make sure it is not singular
        a += 4*eye(n)
        a_ = csc_matrix(a)
        lu = splu(a_)
        b = ones(n)
        x = lu.solve(b)
        assert_almost_equal(dot(a, x), b)

    def test_splu_perm(self):
        # Test the permutation vectors exposed by splu.
        n = 30
        a = random.random((n, n))
        a[a < 0.95] = 0
        # Make a diagonal dominant, to make sure it is not singular
        a += 4*eye(n)
        a_ = csc_matrix(a)
        lu = splu(a_)
        # Check that the permutation indices do belong to [0, n-1].
        for perm in (lu.perm_r, lu.perm_c):
            assert_(all(perm > -1))
            assert_(all(perm < n))
            assert_equal(len(unique(perm)), len(perm))

        # Now make a symmetric, and test that the two permutation vectors are
        # the same
        a += a.T
        a_ = csc_matrix(a)
        lu = splu(a_)
        assert_array_equal(lu.perm_r, lu.perm_c)

    def test_lu_refcount(self):
        # Test that we are keeping track of the reference count with splu.
        n = 30
        a = random.random((n, n))
        a[a < 0.95] = 0
        # Make a diagonal dominant, to make sure it is not singular
        a += 4*eye(n)
        a_ = csc_matrix(a)
        lu = splu(a_)

        # And now test that we don't have a refcount bug
        import gc
        import sys
        rc = sys.getrefcount(lu)
        for attr in ('perm_r', 'perm_c'):
            perm = getattr(lu, attr)
            assert_equal(sys.getrefcount(lu), rc + 1)
            del perm
            assert_equal(sys.getrefcount(lu), rc)


if __name__ == "__main__":
    run_module_suite()
