from numpy.testing import assert_array_almost_equal, assert_array_equal
from pytest import raises as assert_raises

from numpy import array, dot, zeros_like, empty
from numpy.random import random
from scipy.linalg import (
    cholesky, cholesky_banded, cho_solve_banded, cho_factor, cho_solve
)

from scipy.linalg._testutils import assert_no_overwrite

from scipy.conftest import array_api_compatible, skip_if_array_api_backend
from scipy._lib._array_api import xp_assert_close


class TestCholesky:

    # integer dtypes not accepted
    @skip_if_array_api_backend('numpy.array_api')
    @array_api_compatible
    def test_simple(self, xp):
        a = xp.asarray([[8, 2, 3], [2, 9, 3], [3, 3, 6]])
        c = cholesky(a)
        xp_assert_close(c.T @ c, a, rtol=1e-6, check_dtype=False)
        c = c.T
        a = c @ c.T
        xp_assert_close(cholesky(a, lower=True), c)

    # integer dtypes not accepted
    @skip_if_array_api_backend('numpy.array_api')
    @array_api_compatible
    def test_check_finite(self, xp):
        a = xp.asarray([[8, 2, 3], [2, 9, 3], [3, 3, 6]])
        c = cholesky(a, check_finite=False)
        xp_assert_close(c.T @ c, a, rtol=1e-6, check_dtype=False)
        c = c.T
        a = c @ c.T
        xp_assert_close(cholesky(a, lower=True, check_finite=False), c)

    # https://github.com/numpy/numpy/issues/24451
    @skip_if_array_api_backend('numpy.array_api')
    # https://github.com/data-apis/array-api-compat/issues/54
    @skip_if_array_api_backend('cupy')
    @array_api_compatible
    def test_simple_complex(self, xp):
        m = xp.asarray([[3+1j, 3+4j, 5], [0, 2+2j, 2+7j], [0, 0, 7+4j]])
        a = xp.conj(m).T @ m
        c = cholesky(a)
        a1 = xp.conj(c).T @ c
        xp_assert_close(a, a1)
        c = c.T
        a = c @ xp.conj(c).T
        xp_assert_close(cholesky(a, lower=True), c)

    @array_api_compatible
    def test_random(self, xp):
        n = 20
        for k in range(2):
            m = random([n, n])
            for i in range(n):
                m[i, i] = 20*(.1+m[i, i])
            m = xp.asarray(m)
            a = m.T @ m
            c = cholesky(a)
            a1 = c.T @ c
            xp_assert_close(a, a1)
            c = c.T
            a = c @ c.T
            xp_assert_close(cholesky(a, lower=True), c)

    # https://github.com/numpy/numpy/issues/24451
    @skip_if_array_api_backend('numpy.array_api')
    # https://github.com/data-apis/array-api-compat/issues/54
    @skip_if_array_api_backend('cupy')
    @array_api_compatible
    def test_random_complex(self, xp):
        n = 20
        for k in range(2):
            m = random([n, n])+1j*random([n, n])
            for i in range(n):
                m[i, i] = 20*(.1+abs(m[i, i]))
            m = xp.asarray(m)
            a = xp.conj(m).T @ m
            c = cholesky(a)
            a1 = xp.conj(c).T @ c
            xp_assert_close(a, a1)
            c = c.T
            a = c @ xp.conj(c).T
            xp_assert_close(cholesky(a, lower=True), c)


class TestCholeskyBanded:
    """Tests for cholesky_banded() and cho_solve_banded."""

    def test_check_finite(self):
        # Symmetric positive definite banded matrix `a`
        a = array([[4.0, 1.0, 0.0, 0.0],
                   [1.0, 4.0, 0.5, 0.0],
                   [0.0, 0.5, 4.0, 0.2],
                   [0.0, 0.0, 0.2, 4.0]])
        # Banded storage form of `a`.
        ab = array([[-1.0, 1.0, 0.5, 0.2],
                    [4.0, 4.0, 4.0, 4.0]])
        c = cholesky_banded(ab, lower=False, check_finite=False)
        ufac = zeros_like(a)
        ufac[list(range(4)), list(range(4))] = c[-1]
        ufac[(0, 1, 2), (1, 2, 3)] = c[0, 1:]
        assert_array_almost_equal(a, dot(ufac.T, ufac))

        b = array([0.0, 0.5, 4.2, 4.2])
        x = cho_solve_banded((c, False), b, check_finite=False)
        assert_array_almost_equal(x, [0.0, 0.0, 1.0, 1.0])

    def test_upper_real(self):
        # Symmetric positive definite banded matrix `a`
        a = array([[4.0, 1.0, 0.0, 0.0],
                   [1.0, 4.0, 0.5, 0.0],
                   [0.0, 0.5, 4.0, 0.2],
                   [0.0, 0.0, 0.2, 4.0]])
        # Banded storage form of `a`.
        ab = array([[-1.0, 1.0, 0.5, 0.2],
                    [4.0, 4.0, 4.0, 4.0]])
        c = cholesky_banded(ab, lower=False)
        ufac = zeros_like(a)
        ufac[list(range(4)), list(range(4))] = c[-1]
        ufac[(0, 1, 2), (1, 2, 3)] = c[0, 1:]
        assert_array_almost_equal(a, dot(ufac.T, ufac))

        b = array([0.0, 0.5, 4.2, 4.2])
        x = cho_solve_banded((c, False), b)
        assert_array_almost_equal(x, [0.0, 0.0, 1.0, 1.0])

    def test_upper_complex(self):
        # Hermitian positive definite banded matrix `a`
        a = array([[4.0, 1.0, 0.0, 0.0],
                   [1.0, 4.0, 0.5, 0.0],
                   [0.0, 0.5, 4.0, -0.2j],
                   [0.0, 0.0, 0.2j, 4.0]])
        # Banded storage form of `a`.
        ab = array([[-1.0, 1.0, 0.5, -0.2j],
                    [4.0, 4.0, 4.0, 4.0]])
        c = cholesky_banded(ab, lower=False)
        ufac = zeros_like(a)
        ufac[list(range(4)), list(range(4))] = c[-1]
        ufac[(0, 1, 2), (1, 2, 3)] = c[0, 1:]
        assert_array_almost_equal(a, dot(ufac.conj().T, ufac))

        b = array([0.0, 0.5, 4.0-0.2j, 0.2j + 4.0])
        x = cho_solve_banded((c, False), b)
        assert_array_almost_equal(x, [0.0, 0.0, 1.0, 1.0])

    def test_lower_real(self):
        # Symmetric positive definite banded matrix `a`
        a = array([[4.0, 1.0, 0.0, 0.0],
                   [1.0, 4.0, 0.5, 0.0],
                   [0.0, 0.5, 4.0, 0.2],
                   [0.0, 0.0, 0.2, 4.0]])
        # Banded storage form of `a`.
        ab = array([[4.0, 4.0, 4.0, 4.0],
                    [1.0, 0.5, 0.2, -1.0]])
        c = cholesky_banded(ab, lower=True)
        lfac = zeros_like(a)
        lfac[list(range(4)), list(range(4))] = c[0]
        lfac[(1, 2, 3), (0, 1, 2)] = c[1, :3]
        assert_array_almost_equal(a, dot(lfac, lfac.T))

        b = array([0.0, 0.5, 4.2, 4.2])
        x = cho_solve_banded((c, True), b)
        assert_array_almost_equal(x, [0.0, 0.0, 1.0, 1.0])

    def test_lower_complex(self):
        # Hermitian positive definite banded matrix `a`
        a = array([[4.0, 1.0, 0.0, 0.0],
                   [1.0, 4.0, 0.5, 0.0],
                   [0.0, 0.5, 4.0, -0.2j],
                   [0.0, 0.0, 0.2j, 4.0]])
        # Banded storage form of `a`.
        ab = array([[4.0, 4.0, 4.0, 4.0],
                    [1.0, 0.5, 0.2j, -1.0]])
        c = cholesky_banded(ab, lower=True)
        lfac = zeros_like(a)
        lfac[list(range(4)), list(range(4))] = c[0]
        lfac[(1, 2, 3), (0, 1, 2)] = c[1, :3]
        assert_array_almost_equal(a, dot(lfac, lfac.conj().T))

        b = array([0.0, 0.5j, 3.8j, 3.8])
        x = cho_solve_banded((c, True), b)
        assert_array_almost_equal(x, [0.0, 0.0, 1.0j, 1.0])


class TestOverwrite:
    # do we need a version of this which works with xp arrays?
    def test_cholesky(self):
        assert_no_overwrite(cholesky, [(3, 3)])

    def test_cho_factor(self):
        assert_no_overwrite(cho_factor, [(3, 3)])

    def test_cho_solve(self):
        x = array([[2, -1, 0], [-1, 2, -1], [0, -1, 2]])
        xcho = cho_factor(x)
        assert_no_overwrite(lambda b: cho_solve(xcho, b), [(3,)])

    def test_cholesky_banded(self):
        assert_no_overwrite(cholesky_banded, [(2, 3)])

    def test_cho_solve_banded(self):
        x = array([[0, -1, -1], [2, 2, 2]])
        xcho = cholesky_banded(x)
        assert_no_overwrite(lambda b: cho_solve_banded((xcho, False), b),
                            [(3,)])


class TestEmptyArray:
    def test_cho_factor_empty_square(self):
        a = empty((0, 0))
        b = array([])
        c = array([[]])
        d = []
        e = [[]]

        x, _ = cho_factor(a)
        assert_array_equal(x, a)

        for x in ([b, c, d, e]):
            assert_raises(ValueError, cho_factor, x)
