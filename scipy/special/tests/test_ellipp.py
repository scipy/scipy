import numpy as np
from numpy.testing import assert_allclose
from scipy.special import ellipp, ellippinc


def test_ellipp_sympy():
    n = np.array([-0.5, 0.0, 0.3, 1.3, -0.7])
    m = np.array([0.0, 0.2, 0.4, 0.8, 0.5])
    p_scipy = ellipp(n, m)
    # p_sympy = np.fromiter((sy.re(sy.elliptic_pi(ni, mi)).evalf()
    #                        for ni, mi in zip(n, m)), dtype=float)
    p_sympy = np.array([1.28254983, 1.6596236, 2.14879542, -1.7390616, 1.3902519])
    assert_allclose(p_scipy, p_sympy)


def test_ellippinc_sympy():
    n = np.array([-0.5, 0.0, 0.3, 1.3, -0.7])
    m = np.array([0.0, 0.2, 0.4, 0.8, 0.5])
    phi = np.array([-2.5, 0.0, 0.5, 5.5, 7.2])
    p_scipy = ellippinc(phi, n, m)
    # p_sympy = np.fromiter((sy.re(sy.elliptic_pi(ni, phii, mi)).evalf()
    #                        for ni, phii, mi in zip(n, phi, m)), dtype=float)
    p_sympy = np.array([0.60501809, 0.0, 0.52106308, -1.24837946, 0.84754152])
    assert_allclose(p_scipy, p_sympy)
