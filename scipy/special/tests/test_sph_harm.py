import numpy as np
import pytest

from scipy._lib._array_api import make_xp_test_case, xp_assert_close
import scipy.special as sc


lazy_xp_modules = [sc]


class TestSphHarm:
    @pytest.mark.parametrize("n, m", [(0, 0), (1, 1), (3, -2), (7, 4)])
    @make_xp_test_case(sc.sph_harm_y)
    def test_p(self, n, m, xp):
        theta_np = np.linspace(0, np.pi, 20)[:, None]
        phi_np = np.linspace(0, 2 * np.pi, 25)[None, :]
        theta, phi = xp.asarray(theta_np), xp.asarray(phi_np)

        y_actual = sc.sph_harm_y(n, m, theta, phi)

        p = sc.sph_legendre_p(n, m, theta_np)
        y_desired = p * np.exp(1j * m * phi_np)

        xp_assert_close(y_actual, xp.asarray(y_desired), rtol=1e-7)

    @pytest.mark.slow
    @make_xp_test_case(sc.sph_harm_y_all)
    def test_p_all(self, xp):
        m_max = 20
        n_max = 10

        theta_np = np.linspace(0, np.pi)
        phi_np = np.linspace(0, 2 * np.pi)
        theta_np, phi_np = np.meshgrid(theta_np, phi_np)
        theta, phi = xp.asarray(theta_np), xp.asarray(phi_np)

        y, y_jac, y_hess = sc.sph_harm_y_all(
            n_max, m_max, theta, phi, diff_n=2
        )
        p, p_jac, p_hess = sc.sph_legendre_p_all(
            n_max, m_max, theta_np, diff_n=2
        )

        m = np.concatenate([np.arange(m_max + 1), np.arange(-m_max, 0)])
        m = np.expand_dims(m, axis=(0,) + tuple(range(2, theta_np.ndim + 2)))
        phase = np.exp(1j * m * phi_np)

        xp_assert_close(y, xp.asarray(p * phase), rtol=1e-7)
        xp_assert_close(y_jac[..., 0], xp.asarray(p_jac * phase), rtol=1e-7)
        xp_assert_close(
            y_jac[..., 1], xp.asarray(1j * m * p * phase), rtol=1e-7
        )
        xp_assert_close(
            y_hess[..., 0, 0], xp.asarray(p_hess * phase), rtol=1e-7
        )
        xp_assert_close(
            y_hess[..., 0, 1], xp.asarray(1j * m * p_jac * phase), rtol=1e-7
        )
        xp_assert_close(y_hess[..., 1, 0], y_hess[..., 0, 1], rtol=1e-7)
        xp_assert_close(
            y_hess[..., 1, 1], xp.asarray(-m * m * p * phase), rtol=1e-7
        )

    @pytest.mark.parametrize("n_max", [7, 10, 50])
    @pytest.mark.parametrize("m_max", [1, 4, 5, 9, 14])
    @make_xp_test_case(sc.sph_harm_y_all)
    def test_all(self, n_max, m_max, xp):
        theta_np = np.linspace(0, np.pi)
        phi_np = np.linspace(0, 2 * np.pi)
        theta, phi = xp.asarray(theta_np), xp.asarray(phi_np)

        n = np.arange(n_max + 1)
        n = np.expand_dims(n, axis=tuple(range(1, theta_np.ndim + 2)))

        m = np.concatenate([np.arange(m_max + 1), np.arange(-m_max, 0)])
        m = np.expand_dims(m, axis=(0,) + tuple(range(2, theta_np.ndim + 2)))

        y_actual = sc.sph_harm_y_all(n_max, m_max, theta, phi)
        y_desired = sc.sph_harm_y(n, m, theta_np, phi_np)

        xp_assert_close(y_actual, xp.asarray(y_desired), rtol=1e-5)
