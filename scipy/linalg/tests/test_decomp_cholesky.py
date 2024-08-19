import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_allclose, assert_equal
from pytest import raises as assert_raises

from numpy import array, transpose, dot, conjugate, zeros_like, empty
from numpy.random import random
from scipy.linalg import (cholesky, cholesky_banded, cho_solve_banded,
     cho_factor, cho_solve, cholesky_update)

from scipy.linalg._testutils import assert_no_overwrite


class TestCholesky:

    def test_simple(self):
        a = [[8, 2, 3], [2, 9, 3], [3, 3, 6]]
        c = cholesky(a)
        assert_array_almost_equal(dot(transpose(c), c), a)
        c = transpose(c)
        a = dot(c, transpose(c))
        assert_array_almost_equal(cholesky(a, lower=1), c)

    def test_check_finite(self):
        a = [[8, 2, 3], [2, 9, 3], [3, 3, 6]]
        c = cholesky(a, check_finite=False)
        assert_array_almost_equal(dot(transpose(c), c), a)
        c = transpose(c)
        a = dot(c, transpose(c))
        assert_array_almost_equal(cholesky(a, lower=1, check_finite=False), c)

    def test_simple_complex(self):
        m = array([[3+1j, 3+4j, 5], [0, 2+2j, 2+7j], [0, 0, 7+4j]])
        a = dot(transpose(conjugate(m)), m)
        c = cholesky(a)
        a1 = dot(transpose(conjugate(c)), c)
        assert_array_almost_equal(a, a1)
        c = transpose(c)
        a = dot(c, transpose(conjugate(c)))
        assert_array_almost_equal(cholesky(a, lower=1), c)

    def test_random(self):
        n = 20
        for k in range(2):
            m = random([n, n])
            for i in range(n):
                m[i, i] = 20*(.1+m[i, i])
            a = dot(transpose(m), m)
            c = cholesky(a)
            a1 = dot(transpose(c), c)
            assert_array_almost_equal(a, a1)
            c = transpose(c)
            a = dot(c, transpose(c))
            assert_array_almost_equal(cholesky(a, lower=1), c)

    def test_random_complex(self):
        n = 20
        for k in range(2):
            m = random([n, n])+1j*random([n, n])
            for i in range(n):
                m[i, i] = 20*(.1+abs(m[i, i]))
            a = dot(transpose(conjugate(m)), m)
            c = cholesky(a)
            a1 = dot(transpose(conjugate(c)), c)
            assert_array_almost_equal(a, a1)
            c = transpose(c)
            a = dot(c, transpose(conjugate(c)))
            assert_array_almost_equal(cholesky(a, lower=1), c)

    @pytest.mark.xslow
    def test_int_overflow(self):
       # regression test for
       # https://github.com/scipy/scipy/issues/17436
       # the problem was an int overflow in zeroing out
       # the unused triangular part
       n = 47_000
       x = np.eye(n, dtype=np.float64, order='F')
       x[:4, :4] = np.array([[4, -2, 3, -1],
                             [-2, 4, -3, 1],
                             [3, -3, 5, 0],
                             [-1, 1, 0, 5]])

       cholesky(x, check_finite=False, overwrite_a=True)  # should not segfault

    @pytest.mark.parametrize('dt', [int, float, np.float32, complex, np.complex64])
    @pytest.mark.parametrize('dt_b', [int, float, np.float32, complex, np.complex64])
    def test_empty(self, dt, dt_b):
        a = empty((0, 0), dtype=dt)

        c = cholesky(a)
        assert c.shape == (0, 0)
        assert c.dtype == cholesky(np.eye(2, dtype=dt)).dtype

        c_and_lower = (c, True)
        b = np.asarray([], dtype=dt_b)
        x = cho_solve(c_and_lower, b)
        assert x.shape == (0,)
        assert x.dtype == cho_solve((np.eye(2, dtype=dt), True),
                                     np.ones(2, dtype=dt_b)).dtype

        b = empty((0, 0), dtype=dt_b)
        x = cho_solve(c_and_lower, b)
        assert x.shape == (0, 0)
        assert x.dtype == cho_solve((np.eye(2, dtype=dt), True),
                                     np.ones(2, dtype=dt_b)).dtype

        a1 = array([])
        a2 = array([[]])
        a3 = []
        a4 = [[]]
        for x in ([a1, a2, a3, a4]):
            assert_raises(ValueError, cholesky, x)


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

    @pytest.mark.parametrize('dt', [int, float, np.float32, complex, np.complex64])
    @pytest.mark.parametrize('dt_b', [int, float, np.float32, complex, np.complex64])
    def test_empty(self, dt, dt_b):
        ab = empty((0, 0), dtype=dt)

        cb = cholesky_banded(ab)
        assert cb.shape == (0, 0)

        m = cholesky_banded(np.array([[0, 0], [1, 1]], dtype=dt))
        assert cb.dtype == m.dtype

        cb_and_lower = (cb, True)
        b = np.asarray([], dtype=dt_b)
        x = cho_solve_banded(cb_and_lower, b)
        assert x.shape == (0,)

        dtype_nonempty = cho_solve_banded((m, True), np.ones(2, dtype=dt_b)).dtype
        assert x.dtype == dtype_nonempty

        b = empty((0, 0), dtype=dt_b)
        x = cho_solve_banded(cb_and_lower, b)
        assert x.shape == (0, 0)
        assert x.dtype == dtype_nonempty


class TestOverwrite:
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

class TestChoFactor:
    @pytest.mark.parametrize('dt', [int, float, np.float32, complex, np.complex64])
    def test_empty(self, dt):
        a = np.empty((0, 0), dtype=dt)
        x, lower = cho_factor(a)

        assert x.shape == (0, 0)

        xx, lower = cho_factor(np.eye(2, dtype=dt))
        assert x.dtype == xx.dtype


def gen_input(n=1000, seed=72479):
    rng = np.random.default_rng(seed)
    A = rng.random((n, n))
    A = A + A.T + n * np.eye(n)
    R = cholesky(A)
    z = rng.random(n)
    # Uncomment the following to make it occasionally trip up
    #z *= 0.95 * np.sqrt(n)
    return R, z


class TestCholeskyUpdate:
    def test_update(self):
        # Test rank-1 update
        n = 1000
        R, z = gen_input(n)
        U = cholesky_update(R, z)
        tol = n * np.spacing(100.)
        assert_allclose(U.T.dot(U) - (R.T.dot(R) + np.outer(z, z)),
                        np.zeros([n, n]), atol=tol)

    def test_downdate(self):
        # Test rank-1 downdate
        n = 1000
        R, z = gen_input(n)
        D = cholesky_update(R, z, downdate=True)
        tol = n * np.spacing(100.)
        assert_allclose(D.T.dot(D) - (R.T.dot(R) - np.outer(z, z)),
                        np.zeros([n, n]), atol=tol)

    def test_list_input(self):
        # Test list as input
        n = 1000
        R, z = gen_input(n)
        D = cholesky_update(list(R), list(z))
        assert_allclose(D.T.dot(D), (R.T.dot(R) + np.outer(z, z)))

    @pytest.mark.parametrize('bad_value', [np.nan, np.inf])
    def test_finite(self, bad_value):
        # Test with nan and inf in the input
        n = 10
        R0, z = gen_input(n)
        message = 'array must not contain infs or NaNs'

        R = R0.copy()
        R[0, 0] = bad_value
        with pytest.raises(ValueError, match=message):
            cholesky_update(R, z)

        z[0] = bad_value
        with pytest.raises(ValueError, match=message):
            cholesky_update(R0, z)

    def test_dimensions(self):
        # Tests related to dimension of both array and vector
        n = 10
        R, z = gen_input(n)

        message = "Expected 2D array to be updated."
        with pytest.raises(ValueError, match=message):  # too many dims in R
            cholesky_update(np.expand_dims(R, 0), z)
        with pytest.raises(ValueError, match=message):  # too few dims in R
            cholesky_update(z, z)

        message = "Expected 1D update vector."
        with pytest.raises(ValueError, match=message):
            cholesky_update(R, R)

        message = "Input needs to be a square matrix."
        with pytest.raises(ValueError, match=message):
            cholesky_update(R[:, :-1], z)

        message = "Input z has to have same length as the number of rows as input R."
        with pytest.raises(ValueError, match=message):
            cholesky_update(R[:-1, :-1], z)

    def test_param_overwrite(self):
        # Test param overwrite
        n = 100

        # Test with overwrite_R=False and overwrite_z=False
        R, z = gen_input(n)
        R_copy = R.copy()
        z_copy = z.copy()
        cholesky_update(R, z)
        assert_equal(R, R_copy)
        assert_equal(z, z_copy)

        # Test with overwrite_R=True and overwrite_z=True
        R, z = gen_input(n)
        R_copy = R.copy()
        z_copy = z.copy()
        cholesky_update(R, z, overwrite_R=True, overwrite_z=True)
        assert_raises(AssertionError, assert_equal, R, R_copy)
        assert_raises(AssertionError, assert_equal, z, z_copy)

    def test_param_lower(self):
        # Test the 'lower' parameter
        n = 100
        R, z = gen_input(n)
        L = R.T
        assert_allclose(cholesky_update(R, z),
                        cholesky_update(L, z, lower=True).T)

    def test_param_eps(self):
        # Test the 'eps' parameter
        n = 100
        R, z = gen_input(n)
        eps_val = n * np.spacing(1.)
        assert_allclose(cholesky_update(R, z),
                        cholesky_update(R, z, eps=eps_val))
