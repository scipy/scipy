import numpy as np

from numpy.testing import (assert_allclose, assert_array_almost_equal_nulp,
                           assert_equal, assert_array_equal)
from pytest import raises as assert_raises

from scipy.linalg import hilbert, svd
from scipy.sparse import csc_matrix, csr_matrix, isspmatrix
from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import svds
from scipy.sparse.linalg.eigen.arpack import ArpackNoConvergence


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


def svd_test_input_check():
    x = np.array([[1, 2, 3],
                  [3, 4, 3],
                  [1, 0, 2],
                  [0, 0, 1]], float)

    assert_raises(ValueError, svds, x, k=-1)
    assert_raises(ValueError, svds, x, k=0)
    assert_raises(ValueError, svds, x, k=10)
    assert_raises(ValueError, svds, x, k=x.shape[0])
    assert_raises(ValueError, svds, x, k=x.shape[1])
    assert_raises(ValueError, svds, x.T, k=x.shape[0])
    assert_raises(ValueError, svds, x.T, k=x.shape[1])


def test_svd_simple_real():
    x = np.array([[1, 2, 3],
                  [3, 4, 3],
                  [1, 0, 2],
                  [0, 0, 1]], float)
    y = np.array([[1, 2, 3, 8],
                  [3, 4, 3, 5],
                  [1, 0, 2, 3],
                  [0, 0, 1, 0]], float)
    z = csc_matrix(x)

    for solver in [None, 'arpack', 'lobpcg']:
        for m in [x.T, x, y, z, z.T]:
            for k in range(1, min(m.shape)):
                u, s, vh = sorted_svd(m, k)
                su, ss, svh = svds(m, k, solver=solver)

                m_hat = svd_estimate(u, s, vh)
                sm_hat = svd_estimate(su, ss, svh)

                assert_array_almost_equal_nulp(m_hat, sm_hat, nulp=1000)


def test_svd_simple_complex():
    x = np.array([[1, 2, 3],
                  [3, 4, 3],
                  [1 + 1j, 0, 2],
                  [0, 0, 1]], complex)
    y = np.array([[1, 2, 3, 8 + 5j],
                  [3 - 2j, 4, 3, 5],
                  [1, 0, 2, 3],
                  [0, 0, 1, 0]], complex)
    z = csc_matrix(x)

    for solver in [None, 'arpack', 'lobpcg']:
        for m in [x, x.T.conjugate(), x.T, y, y.conjugate(), z, z.T]:
            for k in range(1, min(m.shape) - 1):
                u, s, vh = sorted_svd(m, k)
                su, ss, svh = svds(m, k, solver=solver)

                m_hat = svd_estimate(u, s, vh)
                sm_hat = svd_estimate(su, ss, svh)

                assert_array_almost_equal_nulp(m_hat, sm_hat, nulp=1000)


def test_svd_maxiter():
    # check that maxiter works as expected
    x = hilbert(6)
    # ARPACK shouldn't converge on such an ill-conditioned matrix with just
    # one iteration
    assert_raises(ArpackNoConvergence, svds, x, 1, maxiter=1, ncv=3)
    # but 100 iterations should be more than enough
    u, s, vt = svds(x, 1, maxiter=100, ncv=3)
    assert_allclose(s, [1.7], atol=0.5)


def test_svd_return():
    # check that the return_singular_vectors parameter works as expected
    x = hilbert(6)
    _, s, _ = sorted_svd(x, 2)
    ss = svds(x, 2, return_singular_vectors=False)
    assert_allclose(s, ss)


def test_svd_which():
    # check that the which parameter works as expected
    x = hilbert(6)
    for which in ['LM', 'SM']:
        _, s, _ = sorted_svd(x, 2, which=which)
        for solver in [None, 'arpack', 'lobpcg']:
            ss = svds(x, 2, which=which, return_singular_vectors=False,
                      solver=solver)
            ss.sort()
            assert_allclose(s, ss, atol=np.sqrt(1e-15))


def test_svd_v0():
    # check that the v0 parameter works as expected
    x = np.array([[1, 2, 3, 4], [5, 6, 7, 8]], float)

    for solver in [None, 'arpack', 'lobpcg']:
        u, s, vh = svds(x, 1, solver=solver)
        u2, s2, vh2 = svds(x, 1, v0=u[:, 0], solver=solver)

        assert_allclose(s, s2, atol=np.sqrt(1e-15))


def _check_svds(A, k, U, s, VH):
    n, m = A.shape

    # Check shapes.
    assert_equal(U.shape, (n, k))
    assert_equal(s.shape, (k,))
    assert_equal(VH.shape, (k, m))

    # Check that the original matrix can be reconstituted.
    A_rebuilt = (U*s).dot(VH)
    assert_equal(A_rebuilt.shape, A.shape)
    assert_allclose(A_rebuilt, A)

    # Check that U is a semi-orthogonal matrix.
    UH_U = np.dot(U.T.conj(), U)
    assert_equal(UH_U.shape, (k, k))
    assert_allclose(UH_U, np.identity(k), atol=1e-12)

    # Check that V is a semi-orthogonal matrix.
    VH_V = np.dot(VH, VH.T.conj())
    assert_equal(VH_V.shape, (k, k))
    assert_allclose(VH_V, np.identity(k), atol=1e-12)


def test_svd_LM_ones_matrix():
    # Check that svds can deal with matrix_rank less than k in LM mode.
    k = 3
    for n, m in (6, 5), (5, 5), (5, 6):
        for t in float, complex:
            A = np.ones((n, m), dtype=t)
            for solver in [None, 'arpack', 'lobpcg']:
                U, s, VH = svds(A, k, solver=solver)

                # Check some generic properties of svd.
                _check_svds(A, k, U, s, VH)

                # Check that the largest singular value is near sqrt(n*m)
                # and the other singular values have been forced to zero.
                assert_allclose(np.max(s), np.sqrt(n*m))
                assert_array_equal(sorted(s)[:-1], 0)


def test_svd_LM_zeros_matrix():
    # Check that svds can deal with matrices containing only zeros.
    k = 1
    for n, m in (3, 4), (4, 4), (4, 3):
        for t in float, complex:
            A = np.zeros((n, m), dtype=t)
            for solver in [None, 'arpack', 'lobpcg']:
                U, s, VH = svds(A, k, solver=solver)

                # Check some generic properties of svd.
                _check_svds(A, k, U, s, VH)

                # Check that the singular values are zero.
                assert_array_equal(s, 0)


def test_svd_LM_zeros_matrix_gh_3452():
    # Regression test for a github issue.
    # https://github.com/scipy/scipy/issues/3452
    # Note that for complex dype the size of this matrix is too small for k=1.
    n, m, k = 4, 2, 1
    A = np.zeros((n, m))
    for solver in [None, 'arpack', 'lobpcg']:
        U, s, VH = svds(A, k, solver=solver)

        # Check some generic properties of svd.
        _check_svds(A, k, U, s, VH)

        # Check that the singular values are zero.
        assert_array_equal(s, 0)


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


def test_svd_linop():
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

        v0 = np.ones(min(A.shape))

        for solver in [None, 'arpack', 'lobpcg']:
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

        for solver in [None, 'arpack', 'lobpcg']:
            U1, s1, VH1 = reorder(svds(A, k, which="SM", solver=solver))
            U2, s2, VH2 = reorder(svds(L, k, which="SM", solver=solver))

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

                for solver in [None, 'arpack', 'lobpcg']:
                    U1, s1, VH1 = reorder(svds(A, k, which="LM",
                                               solver=solver))
                    U2, s2, VH2 = reorder(svds(L, k, which="LM",
                                               solver=solver))

                    assert_allclose(np.abs(U1), np.abs(U2), rtol=eps)
                    assert_allclose(s1, s2, rtol=eps)
                    assert_allclose(np.abs(VH1), np.abs(VH2), rtol=eps)
                    assert_allclose(np.dot(U1, np.dot(np.diag(s1), VH1)),
                                    np.dot(U2, np.dot(np.diag(s2), VH2)),
                                    rtol=eps)


def test_svds_partial_return():
    x = np.array([[1, 2, 3],
                  [3, 4, 3],
                  [1, 0, 2],
                  [0, 0, 1]], float)
    # test vertical matrix
    z = csr_matrix(x)
    vh_full = svds(z, 2)[-1]
    vh_partial = svds(z, 2, return_singular_vectors='vh')[-1]
    dvh = np.linalg.norm(np.abs(vh_full) - np.abs(vh_partial))
    if dvh > 1e-10:
        raise AssertionError('right eigenvector matrices differ when using '
                             'return_singular_vectors parameter')
    if svds(z, 2, return_singular_vectors='vh')[0] is not None:
        raise AssertionError('left eigenvector matrix was computed when it '
                             'should not have been')
    # test horizontal matrix
    z = csr_matrix(x.T)
    u_full = svds(z, 2)[0]
    u_partial = svds(z, 2, return_singular_vectors='vh')[0]
    du = np.linalg.norm(np.abs(u_full) - np.abs(u_partial))
    if du > 1e-10:
        raise AssertionError('left eigenvector matrices differ when using '
                             'return_singular_vectors parameter')
    if svds(z, 2, return_singular_vectors='u')[-1] is not None:
        raise AssertionError('right eigenvector matrix was computed when it '
                             'should not have been')


def test_svds_wrong_eigen_type():
    # Regression test for a github issue.
    # https://github.com/scipy/scipy/issues/4590
    # Function was not checking for eigenvalue type and unintended
    # values could be returned.
    x = np.array([[1, 2, 3],
                  [3, 4, 3],
                  [1, 0, 2],
                  [0, 0, 1]], float)
    assert_raises(ValueError, svds, x, 1, which='LA')
