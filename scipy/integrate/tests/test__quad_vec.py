from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_allclose, assert_equal

from scipy.integrate import quad_vec


def test_quad_vec_simple():
    n = np.arange(10)
    f = lambda x: x**n
    for epsabs in [0.1, 1e-3, 1e-6]:
        exact = 2**(n+1)/(n + 1)

        res, err = quad_vec(f, 0, 2, norm='max', epsabs=epsabs)
        assert_allclose(res, exact, rtol=0, atol=epsabs)

        res, err = quad_vec(f, 0, 2, norm='2', epsabs=epsabs)
        assert np.linalg.norm(res - exact) < epsabs

        res, err = quad_vec(f, 0, 2, norm='max', epsabs=epsabs, points=(0.5, 1.0))
        assert_allclose(res, exact, rtol=0, atol=epsabs)

        res, err, *rest = quad_vec(f, 0, 2, norm='max',
                                   epsrel=1e-8, epsabs=epsabs,
                                   full_output=True,
                                   limit=10000, quadrature='trapz')
        assert_allclose(res, exact, rtol=0, atol=epsabs)


def test_quad_vec_simple_inf():
    f = lambda x: 1 / (1 + np.float64(x)**2)

    for epsabs in [0.1, 1e-3, 1e-6]:
        res, err = quad_vec(f, 0, np.inf, norm='max', epsabs=epsabs)
        assert_allclose(res, np.pi/2, rtol=0, atol=max(epsabs, err))

        res, err = quad_vec(f, 0, -np.inf, norm='max', epsabs=epsabs)
        assert_allclose(res, -np.pi/2, rtol=0, atol=max(epsabs, err))

        res, err = quad_vec(f, -np.inf, 0, norm='max', epsabs=epsabs)
        assert_allclose(res, np.pi/2, rtol=0, atol=max(epsabs, err))

        res, err = quad_vec(f, np.inf, 0, norm='max', epsabs=epsabs)
        assert_allclose(res, -np.pi/2, rtol=0, atol=max(epsabs, err))

        res, err = quad_vec(f, -np.inf, np.inf, norm='max', epsabs=epsabs)
        assert_allclose(res, np.pi, rtol=0, atol=max(epsabs, err))

        res, err = quad_vec(f, np.inf, -np.inf, norm='max', epsabs=epsabs)
        assert_allclose(res, -np.pi, rtol=0, atol=max(epsabs, err))

        res, err = quad_vec(f, np.inf, np.inf, norm='max', epsabs=epsabs)
        assert_allclose(res, 0, rtol=0, atol=max(epsabs, err))

        res, err = quad_vec(f, -np.inf, -np.inf, norm='max', epsabs=epsabs)
        assert_allclose(res, 0, rtol=0, atol=max(epsabs, err))

        res, err = quad_vec(f, 0, np.inf, norm='max', epsabs=epsabs, points=(1.0, 2.0))
        assert_allclose(res, np.pi/2, rtol=0, atol=max(epsabs, err))

    f = lambda x: np.sin(x + 2) / (1 + x**2)
    exact = np.pi / np.e * np.sin(2)
    epsabs = 1e-5

    res, err, info = quad_vec(f, -np.inf, np.inf, limit=1000, norm='max', epsabs=epsabs,
                              full_output=True)
    assert_allclose(res, exact, rtol=0, atol=max(epsabs, err))


def _lorenzian(x):
    return 1 / (1 + x**2)


def test_quad_vec_pool():
    from multiprocessing.dummy import Pool

    f = _lorenzian
    res, err = quad_vec(f, -np.inf, np.inf, norm='max', epsabs=1e-4, workers=4)
    assert_allclose(res, np.pi, rtol=0, atol=1e-4)

    with Pool(10) as pool:
        f = lambda x: 1 / (1 + x**2)
        res, err = quad_vec(f, -np.inf, np.inf, norm='max', epsabs=1e-4, workers=pool.map)
        assert_allclose(res, np.pi, rtol=0, atol=1e-4)


def test_num_eval():
    def f(x):
        count[0] += 1
        return x**5

    for q in ['gk21', 'trapz']:
        count = [0]
        res = quad_vec(f, 0, 1, norm='max', full_output=True, quadrature=q)
        assert res[2].neval == count[0]


def test_info():
    def f(x):
        return np.ones((3, 2, 1))

    res, err, info = quad_vec(f, 0, 1, norm='max', full_output=True)

    assert info.success == True
    assert info.status == 0
    assert info.message == 'Target precision reached.'
    assert info.neval > 0
    assert info.intervals.shape[1] == 2
    assert info.integrals.shape == (info.intervals.shape[0], 3, 2, 1)
    assert info.errors.shape == (info.intervals.shape[0],)


def test_nan_inf():
    def f_nan(x):
        return np.nan

    def f_inf(x):
        return np.inf if x < 0.1 else 1/x

    res, err, info = quad_vec(f_nan, 0, 1, full_output=True)
    assert info.status == 3

    res, err, info = quad_vec(f_inf, 0, 1, full_output=True)
    assert info.status == 3
