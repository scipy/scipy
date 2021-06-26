import re
import numpy as np

from numpy.testing import (assert_allclose, assert_array_almost_equal_nulp,
                           assert_equal, assert_array_equal)
from pytest import raises as assert_raises
import pytest

from scipy.linalg import hilbert, svd
from scipy.sparse import csc_matrix, csr_matrix, isspmatrix
from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import svds
from scipy.sparse.linalg.eigen.arpack import ArpackNoConvergence


# --- Helper Functions / Classes ---


def sorted_svd(m, k, which='LM'):
    # Compute svd of a dense matrix m, and return singular vectors/values
    # sorted.
    if isspmatrix(m):
        m = m.todense()
    u, s, vh = svd(m)
    if which == 'LM':
        ii = np.argsort(s)[-k:]
    elif which == 'SM':
        ii = np.argsort(s)[:k]
    else:
        raise ValueError("unknown which=%r" % (which,))

    return u[:, ii], s[ii], vh[ii]


def svd_estimate(u, s, vh):
    return np.dot(u, np.dot(np.diag(s), vh))


def _check_svds(A, k, u, s, vh, which="LM",
                check_usvh_A=True, check_svd=False):
    n, m = A.shape
    atol = 1e-12

    # Check shapes.
    assert_equal(u.shape, (n, k))
    assert_equal(s.shape, (k,))
    assert_equal(vh.shape, (k, m))

    # Check that the original matrix can be reconstituted.
    A_rebuilt = (u*s).dot(vh)
    assert_equal(A_rebuilt.shape, A.shape)
    if check_usvh_A:
        assert_allclose(A_rebuilt, A)

    # Check that u is a semi-orthogonal matrix.
    uh_u = np.dot(u.T.conj(), u)
    assert_equal(uh_u.shape, (k, k))
    assert_allclose(uh_u, np.identity(k), atol=atol)

    # Check that V is a semi-orthogonal matrix.
    vh_v = np.dot(vh, vh.T.conj())
    assert_equal(vh_v.shape, (k, k))
    assert_allclose(vh_v, np.identity(k), atol=atol)

    # Check that scipy.sparse.linalg.svds ~ scipy.linalg.svd
    if check_svd:
        u2, s2, vh2 = sorted_svd(A, k, which)
        assert_allclose(np.abs(u), np.abs(u2), atol=atol)
        assert_allclose(s, s2, atol=atol)
        assert_allclose(np.abs(vh), np.abs(vh2), atol=atol)


class CheckingLinearOperator(LinearOperator):
    def __init__(self, A):
        self.A = A
        self.dtype = A.dtype
        self.shape = A.shape

    def _matvec(self, x):
        assert_equal(max(x.shape), np.size(x))
        return self.A.dot(x)

    def _rmatvec(self, x):
        assert_equal(max(x.shape), np.size(x))
        return self.A.T.conjugate().dot(x)


# --- Test Input Validation ---
# Tests input validation on parameters `k` and `which`
# Needs better input validation checks for all other parameters

class SVDSCommonTests:

    solver = None

    # some of these IV tests could run only once, say with solver=None

    def test_svds_input_validation_A_1(self):
        message = "invalid shape"
        with pytest.raises(ValueError, match=message):
            svds([[[1., 2.], [3., 4.]]], k=1, solver=self.solver)

        message = "`A` must not be empty."
        with pytest.raises(ValueError, match=message):
            svds([[]], k=1, solver=self.solver)

    @pytest.mark.parametrize("A", ("hi", 1, [[1, 2], [3, 4]]))
    def test_svds_input_validation_A_2(self, A):
        message = "`A` must be of floating or complex floating data type."
        with pytest.raises(ValueError, match=message):
            svds(A, k=1, solver=self.solver)

    @pytest.mark.parametrize("k", list(range(-1, 6)) + [1.5, "1"])
    def test_svds_input_validation_k_1(self, k):
        rng = np.random.default_rng(0)
        A = rng.random((4, 3))

        if k in {1, 2}:
            u, s, vh = svds(A, k=k, solver=self.solver)
            # partial decomposition, so don't check that u@diag(s)@vh=A;
            # do check that scipy.sparse.linalg.svds ~ scipy.linalg.svd
            _check_svds(A, k, u, s, vh, check_usvh_A=False, check_svd=True)
        else:
            message = ("`k` must be an integer satisfying")
            with pytest.raises(ValueError, match=message):
                svds(A, k=k, solver=self.solver)

    def test_svds_input_validation_k_2(self):
        # I think the stack trace is reasonable when `k` can't be converted
        # to an int.
        message = "int() argument must be a"
        with pytest.raises(TypeError, match=re.escape(message)):
            svds(np.eye(10), k=[], solver=self.solver)

        message = "invalid literal for int()"
        with pytest.raises(ValueError, match=message):
            svds(np.eye(10), k="hi", solver=self.solver)

    @pytest.mark.parametrize("tol", (-1, np.inf, np.nan))
    def test_svds_input_validation_tol_1(self, tol):
        message = "`tol` must be a non-negative floating point value."
        with pytest.raises(ValueError, match=message):
            svds(np.eye(10), tol=tol, solver=self.solver)

    @pytest.mark.parametrize("tol", ([], 'hi'))
    def test_svds_input_validation_tol_2(self, tol):
        # I think the stack trace is reasonable here
        message = "'<' not supported between instances"
        with pytest.raises(TypeError, match=message):
            svds(np.eye(10), tol=tol, solver=self.solver)

    @pytest.mark.parametrize("which", ('LM', 'SM', 'LA', 'SA', 'ekki', 0))
    def test_svds_input_validation_which(self, which):
        # Regression test for a github issue.
        # https://github.com/scipy/scipy/issues/4590
        # Function was not checking for eigenvalue type and unintended
        # values could be returned.
        rng = np.random.default_rng(0)
        A = rng.random((4, 3))
        k = 2

        if which in {'LM', 'SM'}:
            u, s, vh = svds(A, k=k, which=which, solver=self.solver)
            # partial decomposition, so don't check that u@diag(s)@vh=A;
            # do check that scipy.sparse.linalg.svds ~ scipy.linalg.svd
            _check_svds(A, k, u, s, vh, which,
                        check_usvh_A=False, check_svd=True)
        else:
            with pytest.raises(ValueError, match="`which` must be in"):
                svds(A, k=k, which=which, solver=self.solver)

    @pytest.mark.parametrize("transpose", (True, False))
    @pytest.mark.parametrize("n", range(4, 9))
    def test_svds_input_validation_v0_1(self, transpose, n):
        rng = np.random.default_rng(0)
        A = rng.random((5, 7))
        v0 = rng.random(n)
        if transpose:
            A = A.T
        k = 2
        message = "`v0` must have shape"
        if n != 5:
            with pytest.raises(ValueError, match=message):
                svds(A, k=k, v0=v0, solver=self.solver)
        else:
            u, s, vh = svds(A, k=k, v0=v0, solver=self.solver)
            # partial decomposition, so don't check that u@diag(s)@vh=A;
            # do check that scipy.sparse.linalg.svds ~ scipy.linalg.svd
            _check_svds(A, k, u, s, vh, which="LM",
                        check_usvh_A=False, check_svd=True)

    def test_svds_input_validation_v0_2(self):
        A = np.ones((10, 10))
        v0 = np.ones((1, 10))
        message = "`v0` must have shape"
        with pytest.raises(ValueError, match=message):
            svds(A, k=1, v0=v0, solver=self.solver)

    @pytest.mark.parametrize("v0", ("hi", 1, np.ones(10, dtype=int)))
    def test_svds_input_validation_v0_3(self, v0):
        A = np.ones((10, 10))
        message = "`v0` must be of floating or complex floating data type."
        with pytest.raises(ValueError, match=message):
            svds(A, k=1, v0=v0, solver=self.solver)

    @pytest.mark.parametrize("maxiter", (-1, 0, 5.5))
    def test_svds_input_validation_maxiter_1(self, maxiter):
        message = ("`maxiter` must be a non-negative integer.")
        with pytest.raises(ValueError, match=message):
            svds(np.eye(10), maxiter=maxiter, solver=self.solver)

    def test_svds_input_validation_maxiter_2(self):
        # I think the stack trace is reasonable when `k` can't be converted
        # to an int.
        message = "int() argument must be a"
        with pytest.raises(TypeError, match=re.escape(message)):
            svds(np.eye(10), maxiter=[], solver=self.solver)

        message = "invalid literal for int()"
        with pytest.raises(ValueError, match=message):
            svds(np.eye(10), maxiter="hi", solver=self.solver)

    @pytest.mark.parametrize("rsv", ('ekki', 10))
    def test_svds_input_validation_return_singular_vectors(self, rsv):
        message = "`return_singular_vectors` must be in"
        with pytest.raises(ValueError, match=message):
            svds(np.eye(10), return_singular_vectors=rsv, solver=self.solver)

    # --- Test Parameters ---
    # tests that `which`, `v0`, `maxiter`, and `return_singular_vectors` do
    # what they are purported to do. Should also check `ncv` and `tol` somehow.

    def test_svd_which(self):
        # check that the which parameter works as expected
        solver = self.solver

        x = hilbert(6)
        for which in ['LM', 'SM']:
            _, s, _ = sorted_svd(x, 2, which=which)

            ss = svds(x, 2, which=which, return_singular_vectors=False,
                      solver=solver)
            ss.sort()
            assert_allclose(s, ss, atol=np.sqrt(1e-15))

    def test_svd_v0(self):
        # check that the v0 parameter works as expected
        solver = self.solver

        x = np.array([[1, 2, 3, 4], [5, 6, 7, 8]], float)

        u, s, vh = svds(x, 1, solver=solver)
        u2, s2, vh2 = svds(x, 1, v0=u[:, 0], solver=solver)

        assert_allclose(s, s2, atol=np.sqrt(1e-15))

    def test_svd_return(self):
        # check that the return_singular_vectors parameter works as expected
        solver = self.solver

        x = hilbert(6)
        _, s, _ = sorted_svd(x, 2)
        ss = svds(x, 2, solver=solver, return_singular_vectors=False)
        assert_allclose(s, ss)

    def test_svds_partial_return(self):
        solver = self.solver

        x = np.array([[1, 2, 3],
                      [3, 4, 3],
                      [1, 0, 2],
                      [0, 0, 1]], float)
        # test vertical matrix
        z = csr_matrix(x)
        vh_full = svds(z, 2, solver=solver)[-1]
        vh_partial = svds(z, 2, return_singular_vectors='vh', solver=solver)[-1]
        dvh = np.linalg.norm(np.abs(vh_full) - np.abs(vh_partial))
        if dvh > 1e-10:
            raise AssertionError('right eigenvector matrices differ when using '
                                 'return_singular_vectors parameter')
        if svds(z, 2, return_singular_vectors='vh', solver=solver)[0] is not None:
            raise AssertionError('left eigenvector matrix was computed when it '
                                 'should not have been')
        # test horizontal matrix
        z = csr_matrix(x.T)
        u_full = svds(z, 2, solver=solver)[0]
        u_partial = svds(z, 2, return_singular_vectors='vh', solver=solver)[0]
        du = np.linalg.norm(np.abs(u_full) - np.abs(u_partial))
        if du > 1e-10:
            raise AssertionError('left eigenvector matrices differ when using '
                                 'return_singular_vectors parameter')
        if svds(z, 2, return_singular_vectors='u', solver=solver)[-1] is not None:
            raise AssertionError('right eigenvector matrix was computed when it '
                                 'should not have been')

    # --- Test Basic Functionality ---
    # Tests the accuracy of each solver for real and complex matrices provided
    # as dense array, sparse matrix, and LinearOperator. Could be written
    # more concisely and use parametrization.

    def test_svd_simple_real(self):
        solver = self.solver
        np.random.seed(0)  # set random seed for generating propack v0

        x = np.array([[1, 2, 3],
                      [3, 4, 3],
                      [1, 0, 2],
                      [0, 0, 1]], float)
        y = np.array([[1, 2, 3, 8],
                      [3, 4, 3, 5],
                      [1, 0, 2, 3],
                      [0, 0, 1, 0]], float)
        z = csc_matrix(x)

        for m in [x.T, x, y, z, z.T]:
            for k in range(1, min(m.shape)):
                u, s, vh = sorted_svd(m, k)
                su, ss, svh = svds(m, k, solver=solver)

                m_hat = svd_estimate(u, s, vh)
                sm_hat = svd_estimate(su, ss, svh)

                assert_array_almost_equal_nulp(
                    m_hat, sm_hat, nulp=1000 if solver != 'propack' else 1436)

    def test_svd_simple_complex(self):
        solver = self.solver

        x = np.array([[1, 2, 3],
                      [3, 4, 3],
                      [1 + 1j, 0, 2],
                      [0, 0, 1]], complex)
        y = np.array([[1, 2, 3, 8 + 5j],
                      [3 - 2j, 4, 3, 5],
                      [1, 0, 2, 3],
                      [0, 0, 1, 0]], complex)
        z = csc_matrix(x)

        for m in [x, x.T.conjugate(), x.T, y, y.conjugate(), z, z.T]:
            for k in range(1, min(m.shape) - 1):
                u, s, vh = sorted_svd(m, k)
                su, ss, svh = svds(m, k, solver=solver)

                m_hat = svd_estimate(u, s, vh)
                sm_hat = svd_estimate(su, ss, svh)

                assert_array_almost_equal_nulp(
                    m_hat, sm_hat, nulp=1000 if solver != 'propack' else 1575)

    def test_svd_linop(self):
        solver = self.solver

        nmks = [(6, 7, 3),
                (9, 5, 4),
                (10, 8, 5)]

        def reorder(args):
            U, s, VH = args
            j = np.argsort(s)
            return U[:, j], s[j], VH[j, :]

        for n, m, k in nmks:
            # Test svds on a LinearOperator.
            A = np.random.RandomState(52).randn(n, m)
            L = CheckingLinearOperator(A)

            if solver == 'propack':
                v0 = np.ones(n)
            else:
                v0 = np.ones(min(A.shape))

            U1, s1, VH1 = reorder(svds(A, k, v0=v0, solver=solver))
            U2, s2, VH2 = reorder(svds(L, k, v0=v0, solver=solver))

            assert_allclose(np.abs(U1), np.abs(U2))
            assert_allclose(s1, s2)
            assert_allclose(np.abs(VH1), np.abs(VH2))
            assert_allclose(np.dot(U1, np.dot(np.diag(s1), VH1)),
                            np.dot(U2, np.dot(np.diag(s2), VH2)))

            # Try again with which="SM".
            A = np.random.RandomState(1909).randn(n, m)
            L = CheckingLinearOperator(A)

            # TODO: arpack crashes when v0=v0, which="SM"
            kwargs = {'v0': v0} if solver not in {None, 'arpack'} else {}
            U1, s1, VH1 = reorder(svds(A, k, which="SM", solver=solver,
                                       **kwargs))
            U2, s2, VH2 = reorder(svds(L, k, which="SM", solver=solver,
                                       **kwargs))

            assert_allclose(np.abs(U1), np.abs(U2))
            assert_allclose(s1, s2)
            assert_allclose(np.abs(VH1), np.abs(VH2))
            assert_allclose(np.dot(U1, np.dot(np.diag(s1), VH1)),
                            np.dot(U2, np.dot(np.diag(s2), VH2)))

            if k < min(n, m) - 1:
                # Complex input and explicit which="LM".
                for (dt, eps) in [(complex, 1e-7), (np.complex64, 1e-3)]:
                    rng = np.random.RandomState(1648)
                    A = (rng.randn(n, m) + 1j * rng.randn(n, m)).astype(dt)
                    L = CheckingLinearOperator(A)

                    U1, s1, VH1 = reorder(svds(A, k, which="LM", solver=solver))
                    U2, s2, VH2 = reorder(svds(L, k, which="LM", solver=solver))

                    assert_allclose(np.abs(U1), np.abs(U2), rtol=eps)
                    assert_allclose(s1, s2, rtol=eps)
                    assert_allclose(np.abs(VH1), np.abs(VH2), rtol=eps)
                    assert_allclose(np.dot(U1, np.dot(np.diag(s1), VH1)),
                                    np.dot(U2, np.dot(np.diag(s2), VH2)),
                                    rtol=eps)

    # --- Test Edge Cases ---
    # Checks a few edge cases. There are obvious ones missing (e.g. empty inpout)
    # but I don't think we need to substantially expand these.

    def test_svd_LM_ones_matrix(self):
        # Check that svds can deal with matrix_rank less than k in LM mode.
        solver = self.solver

        k = 3
        for n, m in (6, 5), (5, 5), (5, 6):
            for t in float, complex:
                A = np.ones((n, m), dtype=t)

                U, s, VH = svds(A, k, solver=solver)

                # Check some generic properties of svd.
                _check_svds(A, k, U, s, VH)

                # Check that the largest singular value is near sqrt(n*m)
                # and the other singular values have been forced to zero.
                assert_allclose(np.max(s), np.sqrt(n*m))
                assert_array_equal(sorted(s)[:-1], 0)

    def test_svd_LM_zeros_matrix(self):
        # Check that svds can deal with matrices containing only zeros.
        solver = self.solver

        k = 1
        for n, m in (3, 4), (4, 4), (4, 3):
            for t in float, complex:
                A = np.zeros((n, m), dtype=t)

                U, s, VH = svds(A, k, solver=solver)

                # Check some generic properties of svd.
                _check_svds(A, k, U, s, VH)

                # Check that the singular values are zero.
                assert_array_equal(s, 0)

    def test_svd_LM_zeros_matrix_gh_3452(self):
        # Regression test for a github issue.
        # https://github.com/scipy/scipy/issues/3452
        # Note that for complex dype the size of this matrix is too small for k=1.
        solver = self.solver

        n, m, k = 4, 2, 1
        A = np.zeros((n, m))

        U, s, VH = svds(A, k, solver=solver)

        # Check some generic properties of svd.
        _check_svds(A, k, U, s, VH)

        # Check that the singular values are zero.
        assert_array_equal(s, 0)


# --- Perform tests with each solver ---

class Test_SVDS_once():
    @pytest.mark.parametrize("solver", ['ekki', object])
    def test_svds_input_validation_solver(self, solver):
        message = "solver must be one of"
        with pytest.raises(ValueError, match=message):
            svds(np.ones((3, 4)), k=2, solver=solver)


class Test_SVDS_ARPACK(SVDSCommonTests):

    def setup_method(self):
        self.solver = 'arpack'

    @pytest.mark.parametrize("ncv", list(range(-1, 8)) + [4.5, "5"])
    def test_svds_input_validation_ncv_1(self, ncv):
        rng = np.random.default_rng(0)
        A = rng.random((6, 7))
        k = 3
        if ncv in {4, 5}:
            u, s, vh = svds(A, k=k, ncv=ncv, solver=self.solver)
        # partial decomposition, so don't check that u@diag(s)@vh=A;
        # do check that scipy.sparse.linalg.svds ~ scipy.linalg.svd
            _check_svds(A, k, u, s, vh, check_usvh_A=False, check_svd=True)
        else:
            message = ("`ncv` must be an integer satisfying")
            with pytest.raises(ValueError, match=message):
                svds(A, k=k, ncv=ncv, solver=self.solver)

    def test_svds_input_validation_ncv_2(self):
        # I think the stack trace is reasonable when `ncv` can't be converted
        # to an int.
        message = "int() argument must be a"
        with pytest.raises(TypeError, match=re.escape(message)):
            svds(np.eye(10), ncv=[], solver=self.solver)

        message = "invalid literal for int()"
        with pytest.raises(ValueError, match=message):
            svds(np.eye(10), ncv="hi", solver=self.solver)

    def test_svd_maxiter(self):
        # check that maxiter works as expected
        x = hilbert(6)
        # ARPACK shouldn't converge on such an ill-conditioned matrix with just
        # one iteration
        assert_raises(ArpackNoConvergence, svds, x, 1, maxiter=1, ncv=3)
        # but 100 iterations should be more than enough
        u, s, vt = svds(x, 1, maxiter=100, ncv=3)
        assert_allclose(s, [1.7], atol=0.5)


class Test_SVDS_LOBPCG(SVDSCommonTests):

    def setup_method(self):
        self.solver = 'lobpcg'


class Test_SVDS_PROPACK(SVDSCommonTests):

    def setup_method(self):
        self.solver = 'propack'
