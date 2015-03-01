from __future__ import division, print_function, absolute_import

import warnings

import numpy as np
from numpy import array, finfo, arange, eye, all, unique, ones, dot, matrix
import numpy.random as random
from numpy.testing import (TestCase, run_module_suite,
        assert_array_almost_equal, assert_raises, assert_almost_equal,
        assert_equal, assert_array_equal, assert_, assert_allclose)

import scipy.linalg
from scipy.linalg import norm, inv
from scipy.sparse import (spdiags, SparseEfficiencyWarning, csc_matrix,
        csr_matrix, isspmatrix, dok_matrix, lil_matrix, bsr_matrix)
from scipy.sparse.linalg.dsolve import (spsolve, use_solver, splu, spilu,
        MatrixRankWarning, _superlu)

warnings.simplefilter('ignore',SparseEfficiencyWarning)

# TODO add more comprehensive tests
use_solver(useUmfpack=False)


def toarray(a):
    if isspmatrix(a):
        return a.toarray()
    else:
        return a


class TestLinsolve(TestCase):
    def test_singular(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=MatrixRankWarning)

            A = csc_matrix((5,5), dtype='d')
            b = array([1, 2, 3, 4, 5],dtype='d')
            x = spsolve(A, b, use_umfpack=False)
            assert_(not np.isfinite(x).any())

    def test_singular_gh_3312(self):
        # "Bad" test case that leads SuperLU to call LAPACK with invalid
        # arguments. Check that it fails moderately gracefully.
        ij = np.array([(17, 0), (17, 6), (17, 12), (10, 13)], dtype=np.int32)
        v = np.array([0.284213, 0.94933781, 0.15767017, 0.38797296])
        A = csc_matrix((v, ij.T), shape=(20, 20))
        b = np.arange(20)

        with warnings.catch_warnings():
            try:
                # should either raise a runtimeerror or return value
                # appropriate for singular input
                x = spsolve(A, b, use_umfpack=False)
                assert_(not np.isfinite(x).any())
            except RuntimeError:
                pass

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

    def test_shape_compatibility(self):
        A = csc_matrix([[1., 0], [0, 2]])
        bs = [
            [1, 6],
            array([1, 6]),
            [[1], [6]],
            array([[1], [6]]),
            csc_matrix([[1], [6]]),
            csr_matrix([[1], [6]]),
            dok_matrix([[1], [6]]),
            bsr_matrix([[1], [6]]),
            array([[1., 2., 3.], [6., 8., 10.]]),
            csc_matrix([[1., 2., 3.], [6., 8., 10.]]),
            csr_matrix([[1., 2., 3.], [6., 8., 10.]]),
            dok_matrix([[1., 2., 3.], [6., 8., 10.]]),
            bsr_matrix([[1., 2., 3.], [6., 8., 10.]]),
            ]

        for b in bs:
            x = np.linalg.solve(A.toarray(), toarray(b))
            for spmattype in [csc_matrix, csr_matrix, dok_matrix, lil_matrix]:
                x1 = spsolve(spmattype(A), b, use_umfpack=True)
                x2 = spsolve(spmattype(A), b, use_umfpack=False)

                # check solution
                if x.ndim == 2 and x.shape[1] == 1:
                    # interprets also these as "vectors"
                    x = x.ravel()

                assert_array_almost_equal(toarray(x1), x, err_msg=repr((b, spmattype, 1)))
                assert_array_almost_equal(toarray(x2), x, err_msg=repr((b, spmattype, 2)))

                # dense vs. sparse output  ("vectors" are always dense)
                if isspmatrix(b) and x.ndim > 1:
                    assert_(isspmatrix(x1), repr((b, spmattype, 1)))
                    assert_(isspmatrix(x2), repr((b, spmattype, 2)))
                else:
                    assert_(isinstance(x1, np.ndarray), repr((b, spmattype, 1)))
                    assert_(isinstance(x2, np.ndarray), repr((b, spmattype, 2)))

                # check output shape
                if x.ndim == 1:
                    # "vector"
                    assert_equal(x1.shape, (A.shape[1],))
                    assert_equal(x2.shape, (A.shape[1],))
                else:
                    # "matrix"
                    assert_equal(x1.shape, x.shape)
                    assert_equal(x2.shape, x.shape)

        A = csc_matrix((3, 3))
        b = csc_matrix((1, 3))
        assert_raises(ValueError, spsolve, A, b)

    def test_ndarray_support(self):
        A = array([[1., 2.], [2., 0.]])
        x = array([[1., 1.], [0.5, -0.5]])
        b = array([[2., 0.], [2., 2.]])

        assert_array_almost_equal(x, spsolve(A, b))

    def test_gssv_badinput(self):
        N = 10
        d = arange(N) + 1.0
        A = spdiags((d, 2*d, d[::-1]), (-3, 0, 5), N, N)

        for spmatrix in (csc_matrix, csr_matrix):
            A = spmatrix(A)
            b = np.arange(N)

            def not_c_contig(x):
                return x.repeat(2)[::2]

            def not_1dim(x):
                return x[:,None]

            def bad_type(x):
                return x.astype(bool)

            def too_short(x):
                return x[:-1]

            badops = [not_c_contig, not_1dim, bad_type, too_short]

            for badop in badops:
                msg = "%r %r" % (spmatrix, badop)
                # Not C-contiguous
                assert_raises((ValueError, TypeError), _superlu.gssv,
                              N, A.nnz, badop(A.data), A.indices, A.indptr,
                              b, int(spmatrix == csc_matrix), err_msg=msg)
                assert_raises((ValueError, TypeError), _superlu.gssv,
                              N, A.nnz, A.data, badop(A.indices), A.indptr,
                              b, int(spmatrix == csc_matrix), err_msg=msg)
                assert_raises((ValueError, TypeError), _superlu.gssv,
                              N, A.nnz, A.data, A.indices, badop(A.indptr),
                              b, int(spmatrix == csc_matrix), err_msg=msg)

    def test_sparsity_preservation(self):
        ident = csc_matrix([
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]])
        b = csc_matrix([
            [0, 1],
            [1, 0],
            [0, 0]])
        x = spsolve(ident, b)
        assert_equal(ident.nnz, 3)
        assert_equal(b.nnz, 2)
        assert_equal(x.nnz, 2)
        assert_allclose(x.A, b.A, atol=1e-12, rtol=1e-12)
    
    def test_dense_return(self):
        ident = csc_matrix([
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]])
        b = csc_matrix([
            [0, 1],
            [1, 0],
            [0, 0]])
        x_dense = spsolve(ident, b, return_dense=True)
        x_sparse = spsolve(ident, b)
        self.assertTrue(isinstance(x_dense, np.ndarray))
        self.assertTrue(isinstance(x_sparse, csc_matrix))
        assert_equal(x_dense, x_sparse.todense())

class TestSplu(object):
    def setUp(self):
        n = 40
        d = arange(n) + 1
        self.n = n
        self.A = spdiags((d, 2*d, d[::-1]), (-3, 0, 5), n, n)
        random.seed(1234)

    def _smoketest(self, spxlu, check, dtype):
        if np.issubdtype(dtype, np.complexfloating):
            A = self.A + 1j*self.A.T
        else:
            A = self.A

        A = A.astype(dtype)
        lu = spxlu(A)

        # Input shapes
        for k in [None, 1, 2, self.n, self.n+2]:
            msg = "k=%r" % (k,)

            if k is None:
                b = random.rand(self.n)
            else:
                b = random.rand(self.n, k)

            if np.issubdtype(dtype, np.complexfloating):
                b = b + 1j*random.rand(*b.shape)
            b = b.astype(dtype)

            x = lu.solve(b)
            check(A, b, x, msg)

            x = lu.solve(b, 'T')
            check(A.T, b, x, msg)

            x = lu.solve(b, 'H')
            check(A.T.conj(), b, x, msg)

    def test_splu_smoketest(self):
        # Check that splu works at all
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=SparseEfficiencyWarning)

            def check(A, b, x, msg=""):
                eps = np.finfo(A.dtype).eps
                r = A * x
                assert_(abs(r - b).max() < 1e3*eps, msg)

            self._smoketest(splu, check, np.float32)
            self._smoketest(splu, check, np.float64)
            self._smoketest(splu, check, np.complex64)
            self._smoketest(splu, check, np.complex128)

    def test_spilu_smoketest(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=SparseEfficiencyWarning)

            errors = []

            def check(A, b, x, msg=""):
                r = A * x
                err = abs(r - b).max()
                assert_(err < 1e-2, msg)
                if b.dtype in (np.float64, np.complex128):
                    errors.append(err)

            self._smoketest(spilu, check, np.float32)
            self._smoketest(spilu, check, np.float64)
            self._smoketest(spilu, check, np.complex64)
            self._smoketest(spilu, check, np.complex128)

            assert_(max(errors) > 1e-5)

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
        # Note: a += a.T relies on undefined behavior.
        a = a + a.T
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
        import sys
        rc = sys.getrefcount(lu)
        for attr in ('perm_r', 'perm_c'):
            perm = getattr(lu, attr)
            assert_equal(sys.getrefcount(lu), rc + 1)
            del perm
            assert_equal(sys.getrefcount(lu), rc)

    def test_bad_inputs(self):
        A = self.A.tocsc()

        assert_raises(ValueError, splu, A[:,:4])
        assert_raises(ValueError, spilu, A[:,:4])

        for lu in [splu(A), spilu(A)]:
            b = random.rand(42)
            B = random.rand(42, 3)
            BB = random.rand(self.n, 3, 9)
            assert_raises(ValueError, lu.solve, b)
            assert_raises(ValueError, lu.solve, B)
            assert_raises(ValueError, lu.solve, BB)
            assert_raises(TypeError, lu.solve,
                          b.astype(np.complex64))
            assert_raises(TypeError, lu.solve,
                          b.astype(np.complex128))

    def test_superlu_dlamch_i386_nan(self):
        # SuperLU 4.3 calls some functions returning floats without
        # declaring them. On i386@linux call convention, this fails to
        # clear floating point registers after call. As a result, NaN
        # can appear in the next floating point operation made.
        #
        # Here's a test case that triggered the issue.
        n = 8
        d = np.arange(n) + 1
        A = spdiags((d, 2*d, d[::-1]), (-3, 0, 5), n, n)
        A = A.astype(np.float32)
        spilu(A)
        A = A + 1j*A
        B = A.A
        assert_(not np.isnan(B).any())

    def test_lu_attr(self):

        def check(dtype, complex_2=False):
            A = self.A.astype(dtype)

            if complex_2:
                A = A + 1j*A.T

            n = A.shape[0]
            lu = splu(A)

            # Check that the decomposition is as advertized

            Pc = np.zeros((n, n))
            Pc[np.arange(n), lu.perm_c] = 1

            Pr = np.zeros((n, n))
            Pr[lu.perm_r, np.arange(n)] = 1

            Ad = A.toarray()
            lhs = Pr.dot(Ad).dot(Pc)
            rhs = (lu.L * lu.U).toarray()

            eps = np.finfo(dtype).eps

            assert_allclose(lhs, rhs, atol=100*eps)

        check(np.float32)
        check(np.float64)
        check(np.complex64)
        check(np.complex128)
        check(np.complex64, True)
        check(np.complex128, True)


if __name__ == "__main__":
    run_module_suite()
