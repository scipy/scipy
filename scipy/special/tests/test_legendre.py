import numpy as np

import pytest
from numpy.testing import (assert_equal, assert_almost_equal, assert_array_almost_equal,
    assert_, assert_allclose)

from scipy import special

# Base polynomials come from Abrahmowitz and Stegan
class TestLegendre:
    def test_legendre(self):
        leg0 = special.legendre(0)
        leg1 = special.legendre(1)
        leg2 = special.legendre(2)
        leg3 = special.legendre(3)
        leg4 = special.legendre(4)
        leg5 = special.legendre(5)
        assert_equal(leg0.c, [1])
        assert_equal(leg1.c, [1,0])
        assert_almost_equal(leg2.c, np.array([3,0,-1])/2.0, decimal=13)
        assert_almost_equal(leg3.c, np.array([5,0,-3,0])/2.0)
        assert_almost_equal(leg4.c, np.array([35,0,-30,0,3])/8.0)
        assert_almost_equal(leg5.c, np.array([63,0,-70,0,15,0])/8.0)

    @pytest.mark.parametrize('n', [1, 2, 3, 4, 5])
    @pytest.mark.parametrize('zr', [0.5241717, 12.80232, -9.699001,
                                    0.5122437, 0.1714377])
    @pytest.mark.parametrize('zi', [9.766818, 0.2999083, 8.24726, -22.84843,
                                    -0.8792666])
    def test_lpn_against_clpmn(self, n, zr, zi):
        reslpn = special.lpn(n, zr + zi*1j)
        resclpmn = special.clpmn(0, n, zr+zi*1j)
        assert_allclose(reslpn[0], resclpmn[0][0])
        assert_allclose(reslpn[1], resclpmn[1][0])


class TestLegendreFunctions:
    def test_clpmn(self):
        z = 0.5+0.3j
        clp = special.clpmn(2, 2, z, 3)
        assert_array_almost_equal(clp,
                   (np.array([[1.0000, z, 0.5*(3*z*z-1)],
                           [0.0000, np.sqrt(z*z-1), 3*z*np.sqrt(z*z-1)],
                           [0.0000, 0.0000, 3*(z*z-1)]]),
                    np.array([[0.0000, 1.0000, 3*z],
                           [0.0000, z/np.sqrt(z*z-1), 3*(2*z*z-1)/np.sqrt(z*z-1)],
                           [0.0000, 0.0000, 6*z]])),
                    7)

    def test_clpmn_close_to_real_2(self):
        eps = 1e-10
        m = 1
        n = 3
        x = 0.5
        clp_plus = special.clpmn(m, n, x+1j*eps, 2)[0][m, n]
        clp_minus = special.clpmn(m, n, x-1j*eps, 2)[0][m, n]
        assert_array_almost_equal(np.array([clp_plus, clp_minus]),
                                  np.array([special.lpmv(m, n, x),
                                         special.lpmv(m, n, x)]),
                                  7)

    def test_clpmn_close_to_real_3(self):
        eps = 1e-10
        m = 1
        n = 3
        x = 0.5
        clp_plus = special.clpmn(m, n, x+1j*eps, 3)[0][m, n]
        clp_minus = special.clpmn(m, n, x-1j*eps, 3)[0][m, n]
        assert_array_almost_equal(np.array([clp_plus, clp_minus]),
                                  np.array([special.lpmv(m, n, x)*np.exp(-0.5j*m*np.pi),
                                         special.lpmv(m, n, x)*np.exp(0.5j*m*np.pi)]),
                                  7)

    def test_clpmn_across_unit_circle(self):
        eps = 1e-7
        m = 1
        n = 1
        x = 1j
        for type in [2, 3]:
            assert_almost_equal(special.clpmn(m, n, x+1j*eps, type)[0][m, n],
                            special.clpmn(m, n, x-1j*eps, type)[0][m, n], 6)

    def test_inf(self):
        for z in (1, -1):
            for n in range(4):
                for m in range(1, n):
                    lp = special.clpmn(m, n, z)
                    assert_(np.isinf(lp[1][1,1:]).all())
                    lp = special.lpmn(m, n, z)
                    assert_(np.isinf(lp[1][1,1:]).all())

    def test_deriv_clpmn(self):
        # data inside and outside of the unit circle
        zvals = [0.5+0.5j, -0.5+0.5j, -0.5-0.5j, 0.5-0.5j,
                 1+1j, -1+1j, -1-1j, 1-1j]
        m = 2
        n = 3
        for type in [2, 3]:
            for z in zvals:
                for h in [1e-3, 1e-3j]:
                    approx_derivative = (special.clpmn(m, n, z+0.5*h, type)[0]
                                         - special.clpmn(m, n, z-0.5*h, type)[0])/h
                    assert_allclose(special.clpmn(m, n, z, type)[1],
                                    approx_derivative,
                                    rtol=1e-4)

    @pytest.mark.parametrize("shape", [(10,), (4, 9), (3, 5, 7, 10)])
    @pytest.mark.parametrize("m_max", [5, 4])
    @pytest.mark.parametrize("n_max", [7, 10])
    def test_lpmn(self, shape, n_max, m_max):
        rng = np.random.default_rng(1234)

        x = rng.uniform(-0.99, 0.99, shape)
        p_all, p_all_jac, p_all_hess = special.lpmn_all(m_max, n_max, x, diff_n = 2)

        n = np.arange(n_max + 1)
        n = np.expand_dims(n, axis = (0,) + tuple(range(2, x.ndim + 2)))

        m = np.concatenate([np.arange(m_max + 1), np.arange(-m_max, 0)])
        m = np.expand_dims(m, axis = tuple(range(1, x.ndim + 2)))

        x = np.expand_dims(x, axis = (0, 1))
        p, p_jac, p_hess = special.lpmn(m, n, x, diff_n = 2, legacy = False)

        np.testing.assert_allclose(p, p_all)
        np.testing.assert_allclose(p_jac, p_all_jac)
        np.testing.assert_allclose(p_hess, p_all_hess)

    @pytest.mark.parametrize("shape", [(10,), (4, 9), (3, 5, 7, 10)])
    def test_lpmn_ode(self, shape):
        rng = np.random.default_rng(1234)

        n = rng.integers(0, 10, shape)
        m = rng.integers(-10, 10, shape)
        x = rng.uniform(-1, 1, shape)

        p, p_jac, p_hess = special.lpmn(m, n, x, diff_n = 2, legacy = False)

        assert p.shape == shape
        assert p_jac.shape == p.shape
        assert p_hess.shape == p_jac.shape

        np.testing.assert_allclose((1 - x * x) * p_hess,
            2 * x * p_jac - (n * (n + 1) - m * m / (1 - x * x)) * p,
            rtol = 1e-05, atol = 1e-08)

    @pytest.mark.parametrize("shape", [(10,), (4, 9), (3, 5, 7)])
    def test_lpmn_all(self, shape):
        rng = np.random.default_rng(1234)

        n_max = 20
        m_max = 20

        x = rng.uniform(-0.99, 0.99, shape)

        p, p_jac, p_hess = special.lpmn_all(m_max, n_max, x, diff_n = 2)

        m = np.concatenate([np.arange(m_max + 1), np.arange(-m_max, 0)])
        n = np.arange(n_max + 1)

        m = np.expand_dims(m, axis = tuple(range(1, x.ndim + 2)))
        n = np.expand_dims(n, axis = (0,) + tuple(range(2, x.ndim + 2)))
        np.testing.assert_allclose((1 - x * x) * p_hess,
            2 * x * p_jac - (n * (n + 1) - m * m / (1 - x * x)) * p,
            rtol = 1e-05, atol = 1e-08)

    @pytest.mark.parametrize("shape", [(10,), (4, 9), (3, 5, 7)])
    def test_lpmn_all_specific(self, shape):
        rng = np.random.default_rng(1234)

        x = rng.uniform(-5, 5, shape)

        p, p_jac = special.lpmn_all(4, 4, x, diff_n = 1)

        np.testing.assert_allclose(p[0, 0], lpmn_desired(0, 0, x))

        np.testing.assert_allclose(p[0, 1], lpmn_desired(0, 1, x))
        np.testing.assert_allclose(p[1, 1], lpmn_desired(1, 1, x))
        np.testing.assert_allclose(p[-1, 1], lpmn_desired(-1, 1, x))

        np.testing.assert_allclose(p[0, 2], lpmn_desired(0, 2, x))
        np.testing.assert_allclose(p[1, 2], lpmn_desired(1, 2, x))
        np.testing.assert_allclose(p[2, 2], lpmn_desired(2, 2, x))
        np.testing.assert_allclose(p[-2, 2], lpmn_desired(-2, 2, x))
        np.testing.assert_allclose(p[-1, 2], lpmn_desired(-1, 2, x))

        np.testing.assert_allclose(p[0, 3], lpmn_desired(0, 3, x))
        np.testing.assert_allclose(p[1, 3], lpmn_desired(1, 3, x))
        np.testing.assert_allclose(p[2, 3], lpmn_desired(2, 3, x))
        np.testing.assert_allclose(p[3, 3], lpmn_desired(3, 3, x))
        np.testing.assert_allclose(p[-3, 3], lpmn_desired(-3, 3, x))
        np.testing.assert_allclose(p[-2, 3], lpmn_desired(-2, 3, x))
        np.testing.assert_allclose(p[-1, 3], lpmn_desired(-1, 3, x))

        np.testing.assert_allclose(p[0, 4], lpmn_desired(0, 4, x))
        np.testing.assert_allclose(p[1, 4], lpmn_desired(1, 4, x))
        np.testing.assert_allclose(p[2, 4], lpmn_desired(2, 4, x))
        np.testing.assert_allclose(p[3, 4], lpmn_desired(3, 4, x))
        np.testing.assert_allclose(p[4, 4], lpmn_desired(4, 4, x))
        np.testing.assert_allclose(p[-4, 4], lpmn_desired(-4, 4, x))
        np.testing.assert_allclose(p[-3, 4], lpmn_desired(-3, 4, x))
        np.testing.assert_allclose(p[-2, 4], lpmn_desired(-2, 4, x))
        np.testing.assert_allclose(p[-1, 4], lpmn_desired(-1, 4, x))

        np.testing.assert_allclose(p_jac[0, 1], lpmn_jac_desired(0, 1, x))
        np.testing.assert_allclose(p_jac[1, 1], lpmn_jac_desired(1, 1, x))
        np.testing.assert_allclose(p_jac[-1, 1], lpmn_jac_desired(-1, 1, x))

        np.testing.assert_allclose(p_jac[0, 2], lpmn_jac_desired(0, 2, x))
        np.testing.assert_allclose(p_jac[1, 2], lpmn_jac_desired(1, 2, x))
        np.testing.assert_allclose(p_jac[2, 2], lpmn_jac_desired(2, 2, x))
        np.testing.assert_allclose(p_jac[-2, 2], lpmn_jac_desired(-2, 2, x))
        np.testing.assert_allclose(p_jac[-1, 2], lpmn_jac_desired(-1, 2, x))

        np.testing.assert_allclose(p_jac[0, 3], lpmn_jac_desired(0, 3, x))
        np.testing.assert_allclose(p_jac[1, 3], lpmn_jac_desired(1, 3, x))
        np.testing.assert_allclose(p_jac[2, 3], lpmn_jac_desired(2, 3, x))
        np.testing.assert_allclose(p_jac[3, 3], lpmn_jac_desired(3, 3, x))
        np.testing.assert_allclose(p_jac[-3, 3], lpmn_jac_desired(-3, 3, x))
        np.testing.assert_allclose(p_jac[-2, 3], lpmn_jac_desired(-2, 3, x))
        np.testing.assert_allclose(p_jac[-1, 3], lpmn_jac_desired(-1, 3, x))

        np.testing.assert_allclose(p_jac[0, 4], lpmn_jac_desired(0, 4, x))
        np.testing.assert_allclose(p_jac[1, 4], lpmn_jac_desired(1, 4, x))
        np.testing.assert_allclose(p_jac[2, 4], lpmn_jac_desired(2, 4, x))
        np.testing.assert_allclose(p_jac[3, 4], lpmn_jac_desired(3, 4, x))
        np.testing.assert_allclose(p_jac[4, 4], lpmn_jac_desired(4, 4, x))
        np.testing.assert_allclose(p_jac[-4, 4], lpmn_jac_desired(-4, 4, x))
        np.testing.assert_allclose(p_jac[-3, 4], lpmn_jac_desired(-3, 4, x))
        np.testing.assert_allclose(p_jac[-2, 4], lpmn_jac_desired(-2, 4, x))
        np.testing.assert_allclose(p_jac[-1, 4], lpmn_jac_desired(-1, 4, x))

    @pytest.mark.parametrize("m_max", [7])
    @pytest.mark.parametrize("n_max", [10])
    @pytest.mark.parametrize("x", [1, -1])
    def test_lpmn_all_limits(self, m_max, n_max, x):
        rng = np.random.default_rng(1234)

        p, p_jac = special.lpmn_all(m_max, n_max, x, diff_n = 1)

        n = np.arange(n_max + 1)

        np.testing.assert_allclose(p_jac[0], pow(x, n + 1) * n * (n + 1) / 2)
        np.testing.assert_allclose(p_jac[1], np.where(n >= 1, pow(x, n) * np.inf, 0))
        np.testing.assert_allclose(p_jac[2], np.where(n >= 2, -pow(x, n + 1) * (n + 2) * (n + 1) * n * (n - 1) / 4, 0))
        np.testing.assert_allclose(p_jac[-2], np.where(n >= 2, -pow(x, n + 1) / 4, 0))
        np.testing.assert_allclose(p_jac[-1], np.where(n >= 1, -pow(x, n) * np.inf, 0))

        for m in range(3, m_max + 1):
            np.testing.assert_allclose(p_jac[m], 0)
            np.testing.assert_allclose(p_jac[-m], 0)

    @pytest.mark.parametrize("m_max", [3, 5, 10])
    @pytest.mark.parametrize("n_max", [10])
    def test_lpmn_legacy(self, m_max, n_max):
        x = 0.5
        p, p_jac = special.lpmn_all(m_max, n_max, x, diff_n = 1)

        p_legacy, p_jac_legacy = special.lpmn(m_max, n_max, x)
        for m in range(m_max + 1):
            np.testing.assert_allclose(p_legacy[m], p[m])

        p_legacy, p_jac_legacy = special.lpmn(-m_max, n_max, x)
        for m in range(m_max + 1):
            np.testing.assert_allclose(p_legacy[m], p[-m])

    @pytest.mark.parametrize("shape", [(1000,), (4, 9), (3, 5, 7)])
    @pytest.mark.parametrize("type", [2, 3])
    @pytest.mark.parametrize("z_min, z_max", [(-10 - 10j, 10 + 10j), (-1, 1), (-10j, 10j)])
    def test_clpmn_all_specific(self, shape, type, z_min, z_max):
        rng = np.random.default_rng(1234)

        z = rng.uniform(z_min.real, z_max.real, shape) + 1j * rng.uniform(z_min.imag, z_max.imag, shape)

        p, p_jac = special.clpmn_all(4, 4, type, z, diff_n = 1)

        np.testing.assert_allclose(p[0, 0], clpmn_desired(0, 0, type, z))

        np.testing.assert_allclose(p[0, 1], clpmn_desired(0, 1, type, z))
        np.testing.assert_allclose(p[1, 1], clpmn_desired(1, 1, type, z))
        np.testing.assert_allclose(p[-1, 1], clpmn_desired(-1, 1, type, z))

        np.testing.assert_allclose(p[0, 2], clpmn_desired(0, 2, type, z))
        np.testing.assert_allclose(p[1, 2], clpmn_desired(1, 2, type, z))
        np.testing.assert_allclose(p[2, 2], clpmn_desired(2, 2, type, z))
        np.testing.assert_allclose(p[-2, 2], clpmn_desired(-2, 2, type, z))
        np.testing.assert_allclose(p[-1, 2], clpmn_desired(-1, 2, type, z))
 
        np.testing.assert_allclose(p[0, 3], clpmn_desired(0, 3, type, z))
        np.testing.assert_allclose(p[1, 3], clpmn_desired(1, 3, type, z))
        np.testing.assert_allclose(p[2, 3], clpmn_desired(2, 3, type, z))
        np.testing.assert_allclose(p[3, 3], clpmn_desired(3, 3, type, z))
        np.testing.assert_allclose(p[-3, 3], clpmn_desired(-3, 3, type, z))
        np.testing.assert_allclose(p[-2, 3], clpmn_desired(-2, 3, type, z))
        np.testing.assert_allclose(p[-1, 3], clpmn_desired(-1, 3, type, z))

        np.testing.assert_allclose(p[0, 4], clpmn_desired(0, 4, type, z))
        np.testing.assert_allclose(p[1, 4], clpmn_desired(1, 4, type, z))
        np.testing.assert_allclose(p[2, 4], clpmn_desired(2, 4, type, z))
        np.testing.assert_allclose(p[3, 4], clpmn_desired(3, 4, type, z))
        np.testing.assert_allclose(p[4, 4], clpmn_desired(4, 4, type, z))
        np.testing.assert_allclose(p[-4, 4], clpmn_desired(-4, 4, type, z))
        np.testing.assert_allclose(p[-3, 4], clpmn_desired(-3, 4, type, z))
        np.testing.assert_allclose(p[-2, 4], clpmn_desired(-2, 4, type, z))
        np.testing.assert_allclose(p[-1, 4], clpmn_desired(-1, 4, type, z))

        np.testing.assert_allclose(p_jac[0, 0], clpmn_jac_desired(0, 0, type, z))

        np.testing.assert_allclose(p_jac[0, 1], clpmn_jac_desired(0, 1, type, z))
        np.testing.assert_allclose(p_jac[1, 1], clpmn_jac_desired(1, 1, type, z))
        np.testing.assert_allclose(p_jac[-1, 1], clpmn_jac_desired(-1, 1, type, z))

        np.testing.assert_allclose(p_jac[0, 2], clpmn_jac_desired(0, 2, type, z))
        np.testing.assert_allclose(p_jac[1, 2], clpmn_jac_desired(1, 2, type, z))
        np.testing.assert_allclose(p_jac[2, 2], clpmn_jac_desired(2, 2, type, z))
        np.testing.assert_allclose(p_jac[-2, 2], clpmn_jac_desired(-2, 2, type, z))
        np.testing.assert_allclose(p_jac[-1, 2], clpmn_jac_desired(-1, 2, type, z))

        np.testing.assert_allclose(p_jac[0, 3], clpmn_jac_desired(0, 3, type, z))
        np.testing.assert_allclose(p_jac[1, 3], clpmn_jac_desired(1, 3, type, z))
        np.testing.assert_allclose(p_jac[2, 3], clpmn_jac_desired(2, 3, type, z))
        np.testing.assert_allclose(p_jac[3, 3], clpmn_jac_desired(3, 3, type, z))
        np.testing.assert_allclose(p_jac[-3, 3], clpmn_jac_desired(-3, 3, type, z))
        np.testing.assert_allclose(p_jac[-2, 3], clpmn_jac_desired(-2, 3, type, z))
        np.testing.assert_allclose(p_jac[-1, 3], clpmn_jac_desired(-1, 3, type, z))

        np.testing.assert_allclose(p_jac[0, 4], clpmn_jac_desired(0, 4, type, z))
        np.testing.assert_allclose(p_jac[1, 4], clpmn_jac_desired(1, 4, type, z))
        np.testing.assert_allclose(p_jac[2, 4], clpmn_jac_desired(2, 4, type, z))
        np.testing.assert_allclose(p_jac[3, 4], clpmn_jac_desired(3, 4, type, z))
        np.testing.assert_allclose(p_jac[4, 4], clpmn_jac_desired(4, 4, type, z))
        np.testing.assert_allclose(p_jac[-4, 4], clpmn_jac_desired(-4, 4, type, z))
        np.testing.assert_allclose(p_jac[-3, 4], clpmn_jac_desired(-3, 4, type, z))
        np.testing.assert_allclose(p_jac[-2, 4], clpmn_jac_desired(-2, 4, type, z))
        np.testing.assert_allclose(p_jac[-1, 4], clpmn_jac_desired(-1, 4, type, z))

    """
    @pytest.mark.parametrize("m_max", [3])
    @pytest.mark.parametrize("n_max", [5])
    @pytest.mark.parametrize("z", [-1])
    def test_clpmn_all_limits(self, m_max, n_max, z):
        rng = np.random.default_rng(1234)

        type = 2

#        p, p_jac = special.clpmn(m_max, n_max, z, type = type, legacy = False)
        p, p_jac = special.clpmn_all(m_max, n_max, type, z, diff_n = 1)

        n = np.arange(n_max + 1)

        np.testing.assert_allclose(p_jac[0], pow(z, n + 1) * n * (n + 1) / 2)
        np.testing.assert_allclose(p_jac[1], np.where(n >= 1, pow(z, n) * np.inf, 0))
        np.testing.assert_allclose(p_jac[2], np.where(n >= 2,
            -pow(z, n + 1) * (n + 2) * (n + 1) * n * (n - 1) / 4, 0))
        np.testing.assert_allclose(p_jac[-2], np.where(n >= 2, -pow(z, n + 1) / 4, 0))
        np.testing.assert_allclose(p_jac[-1], np.where(n >= 1, -pow(z, n) * np.inf, 0))

        for m in range(3, m_max + 1):
            np.testing.assert_allclose(p_jac[m], 0)
            np.testing.assert_allclose(p_jac[-m], 0)
    """

    @pytest.mark.parametrize("shape", [(10,), (4, 9), (3, 5, 7)])
    def test_lpn_ode(self, shape):
        rng = np.random.default_rng(1234)

        n = rng.integers(0, 100, shape)
        x = rng.uniform(-1, 1, shape)

        p, p_jac, p_hess = special.lpn(n, x, diff_n = 2, legacy = False)

        assert p.shape == shape
        assert p_jac.shape == p.shape
        assert p_hess.shape == p_jac.shape

        err = (1 - x * x) * p_hess - 2 * x * p_jac + n * (n + 1) * p
        np.testing.assert_allclose(err, 0, atol = 1e-10)

    @pytest.mark.parametrize("n_max", [1, 2, 4, 8, 16, 32])
    @pytest.mark.parametrize("x_shape", [(10,), (4, 9), (3, 5, 7)])
    def test_lpn_all_ode(self, n_max, x_shape):
        rng = np.random.default_rng(1234)

        x = rng.uniform(-1, 1, x_shape)
        p, p_jac, p_hess = special.lpn_all(n_max, x, diff_n = 2)

        n = np.arange(n_max + 1)
        n = np.expand_dims(n, axis = tuple(range(1, x.ndim + 1)))

        assert p.shape == (len(n),) + x.shape
        assert p_jac.shape == p.shape
        assert p_hess.shape == p_jac.shape

        err = (1 - x * x) * p_hess - 2 * x * p_jac + n * (n + 1) * p
        np.testing.assert_allclose(err, 0, atol = 1e-10)

    def test_lpn_legacy(self):
        p, pd = special.lpn(2, 0.5)
        assert_array_almost_equal(p, [1.00000, 0.50000, -0.12500], 4)
        assert_array_almost_equal(pd, [0.00000, 1.00000, 1.50000], 4)

    def test_lpmv(self):
        lp = special.lpmv(0,2,.5)
        assert_almost_equal(lp,-0.125,7)
        lp = special.lpmv(0,40,.001)
        assert_almost_equal(lp,0.1252678976534484,7)

        # XXX: this is outside the domain of the current implementation,
        #      so ensure it returns a NaN rather than a wrong answer.
        with np.errstate(all='ignore'):
            lp = special.lpmv(-1,-1,.001)
        assert_(lp != 0 or np.isnan(lp))

    def test_lqmn(self):
        lqmnf = special.lqmn(0,2,.5)
        lqf = special.lqn(2,.5)
        assert_array_almost_equal(lqmnf[0][0],lqf[0],4)
        assert_array_almost_equal(lqmnf[1][0],lqf[1],4)

    def test_lqmn_gt1(self):
        """algorithm for real arguments changes at 1.0001
           test against analytical result for m=2, n=1
        """
        x0 = 1.0001
        delta = 0.00002
        for x in (x0-delta, x0+delta):
            lq = special.lqmn(2, 1, x)[0][-1, -1]
            expected = 2/(x*x-1)
            assert_almost_equal(lq, expected)

    def test_lqmn_shape(self):
        a, b = special.lqmn(4, 4, 1.1)
        assert_equal(a.shape, (5, 5))
        assert_equal(b.shape, (5, 5))

        a, b = special.lqmn(4, 0, 1.1)
        assert_equal(a.shape, (5, 1))
        assert_equal(b.shape, (5, 1))

    def test_lqn(self):
        lqf = special.lqn(2,.5)
        assert_array_almost_equal(lqf,(np.array([0.5493, -0.7253, -0.8187]),
                                       np.array([1.3333, 1.216, -0.8427])),4)

    @pytest.mark.parametrize("function", [special.lpn, special.lqn])
    @pytest.mark.parametrize("n", [1, 2, 4, 8, 16, 32])
    @pytest.mark.parametrize("z_complex", [False, True])
    @pytest.mark.parametrize("z_inexact", [False, True])
    @pytest.mark.parametrize(
        "input_shape",
        [
            (), (1, ), (2, ), (2, 1), (1, 2), (2, 2), (2, 2, 1), (2, 2, 2)
        ]
    )
    def test_array_inputs_lxn(self, function, n, z_complex, z_inexact, input_shape):
        """Tests for correct output shapes."""
        rng = np.random.default_rng(1234)
        if z_inexact:
            z = rng.integers(-3, 3, size=input_shape)
        else:
            z = rng.uniform(-1, 1, size=input_shape)

        if z_complex:
            z = 1j * z + 0.5j * z

        P_z, P_d_z = function(n, z)
        assert P_z.shape == (n + 1, ) + input_shape
        assert P_d_z.shape == (n + 1, ) + input_shape

    @pytest.mark.parametrize("function", [special.lqmn])
    @pytest.mark.parametrize(
        "m,n",
        [(0, 1), (1, 2), (1, 4), (3, 8), (11, 16), (19, 32)]
    )
    @pytest.mark.parametrize("z_inexact", [False, True])
    @pytest.mark.parametrize(
        "input_shape", [
            (), (1, ), (2, ), (2, 1), (1, 2), (2, 2), (2, 2, 1)
        ]
    )
    def test_array_inputs_lxmn(self, function, m, n, z_inexact, input_shape):
        """Tests for correct output shapes and dtypes."""
        rng = np.random.default_rng(1234)
        if z_inexact:
            z = rng.integers(-3, 3, size=input_shape)
        else:
            z = rng.uniform(-1, 1, size=input_shape)

        P_z, P_d_z = function(m, n, z)
        assert P_z.shape == (m + 1, n + 1) + input_shape
        assert P_d_z.shape == (m + 1, n + 1) + input_shape

    @pytest.mark.parametrize("function", [special.clpmn, special.lqmn])
    @pytest.mark.parametrize(
        "m,n",
        [(0, 1), (1, 2), (1, 4), (3, 8), (11, 16), (19, 32)]
    )
    @pytest.mark.parametrize(
        "input_shape", [
            (), (1, ), (2, ), (2, 1), (1, 2), (2, 2), (2, 2, 1)
        ]
    )
    def test_array_inputs_clxmn(self, function, m, n, input_shape):
        """Tests for correct output shapes and dtypes."""
        rng = np.random.default_rng(1234)
        z = rng.uniform(-1, 1, size=input_shape)
        z = 1j * z + 0.5j * z

        P_z, P_d_z = function(m, n, z)
        assert P_z.shape == (m + 1, n + 1) + input_shape
        assert P_d_z.shape == (m + 1, n + 1) + input_shape

def clpmn_desired(m, n, type, z):
    type_sign = np.where(type == 3, -1, 1)
    branch_sign = np.where(type == 3, np.where(np.signbit(np.real(z)), 1, -1), -1)

    out11 = type_sign * branch_sign * np.sqrt(np.where(type == 3, z * z - 1, 1 - z * z))

    if (n == 0):
        if (m == 0):
            return np.ones_like(z)

    if (n == 1):
        if (m == 0):
            return z

        if (m == 1):
            return out11

        if (m == -1):
            return -type_sign * out11 / 2

    if (n == 2):
        if (m == 0):
            return (3 * z * z - 1) / 2

        if (m == 1):
            return 3 * z * out11

        if (m == 2):
            return 3 * type_sign * (1 - z * z)

        if (m == -2):
            return type_sign * (1 - z * z) / 8

        if (m == -1):
            return -type_sign * z * out11 / 2

    if (n == 3):
        if (m == 0):
            return (5 * z * z - 3) * z / 2

        if (m == 1):
            return 3 * (5 * z * z - 1) * out11 / 2

        if (m == 2):
            return 15 * type_sign * (1 - z * z) * z

        if (m == 3):
            return 15 * type_sign * (1 - z * z) * out11

        if (m == -3):
            return (z * z - 1) * out11 / 48

        if (m == -2):
            return type_sign * (1 - z * z) * z / 8

        if (m == -1):
            return type_sign * (1 - 5 * z * z) * out11 / 8

    if (n == 4):
        if (m == 0):
            return ((35 * z * z - 30) * z * z + 3) / 8

        if (m == 1):
            return 5 * (7 * z * z - 3) * z * out11 / 2

        if (m == 2):
            return 15 * type_sign * ((8 - 7 * z * z) * z * z - 1) / 2

        if (m == 3):
            return 105 * type_sign * (1 - z * z) * z * out11

        if (m == 4):
            return 105 * np.square(z * z - 1)

        if (m == -1):
            return type_sign * (3 - 7 * z * z) * z * out11 / 8

        if (m == -2):
            return type_sign * ((8 - 7 * z * z) * z * z - 1) / 48

        if (m == -3):
            return (z * z - 1) * z * out11 / 48

        if (m == -4):
            return np.square(z * z - 1) / 384

    raise NotImplementedError

def lpmn_desired(m, n, z):
    type = np.where(np.abs(z) <= 1, 2, 3)
    return clpmn_desired(m, n, type, z)

def clpmn_jac_desired(m, n, type, z):
    type_sign = np.where(type == 3, -1, 1)
    branch_sign = np.where(type == 3, np.where(np.signbit(np.real(z)), 1, -1), -1)

    out11_div_z = -branch_sign / np.sqrt(np.where(type == 3, z * z - 1, 1 - z * z))

    if (n == 0):
        if (m == 0):
            return np.zeros_like(z)

    if (n == 1):
        if (m == 0):
            return np.ones_like(z)

        if (m == 1):
            return z * out11_div_z

        if (m == -1):
            return -type_sign * z * out11_div_z / 2

    if (n == 2):
        if (m == 0):
            return 3 * z

        if (m == 1):
            return 3 * (2 * z * z - 1) * out11_div_z

        if (m == 2):
            return -6 * type_sign * z

        if (m == -2):
            return -type_sign * z / 4

        if (m == -1):
            return type_sign * (1 - 2 * z * z) * out11_div_z / 2

    if (n == 3):
        if (m == 0):
            return 3 * (5 * z * z - 1) / 2

        if (m == 1):
            return 3 * (15 * z * z - 11) * z * out11_div_z / 2

        if (m == 2):
            return 15 * type_sign * (1 - 3 * z * z)

        if (m == 3):
            return 45 * type_sign * (1 - z * z) * z * out11_div_z

        if (m == -3):
            return (z * z - 1) * z * out11_div_z / 16

        if (m == -2):
            return type_sign * (1 - 3 * z * z) / 8

        if (m == -1):
            return type_sign * (11 - 15 * z * z) * z * out11_div_z / 8

    if (n == 4):
        if (m == 0):
            return 5 * (7 * z * z - 3) * z / 2

        if (m == 1):
            return 5 * ((28 * z * z - 27) * z * z + 3) * out11_div_z / 2

        if (m == 2):
            return 30 * type_sign * (4 - 7 * z * z) * z

        if (m == 3):
            return 105 * type_sign * ((5 - 4 * z * z) * z * z - 1) * out11_div_z

        if (m == 4):
            return 420 * (z * z - 1) * z

        if (m == -4):
            return (z * z - 1) * z / 96

        if (m == -3):
            return ((4 * z * z - 5) * z * z + 1) * out11_div_z / 48

        if (m == -2):
            return type_sign * (4 - 7 * z * z) * z / 12

        if (m == -1):
            return type_sign * ((27 - 28 * z * z) * z * z - 3) * out11_div_z / 8

    raise NotImplementedError

def lpmn_jac_desired(m, n, z):
    type = np.where(np.abs(z) <= 1, 2, 3)
    return clpmn_jac_desired(m, n, type, z)