import pytest

import numpy as np
from numpy.testing import assert_allclose, assert_equal

from scipy.integrate import quad_vec

from multiprocessing.dummy import Pool


quadrature_params = pytest.mark.parametrize(
    'quadrature', [None, "gk15", "gk21", "trapezoid"])


@quadrature_params
def test_quad_vec_simple(quadrature):
    n = np.arange(10)
    def f(x):
        return x ** n
    for epsabs in [0.1, 1e-3, 1e-6]:
        if quadrature == 'trapezoid' and epsabs < 1e-4:
            # slow: skip
            continue

        kwargs = dict(epsabs=epsabs, quadrature=quadrature)

        exact = 2**(n+1)/(n + 1)

        res, err = quad_vec(f, 0, 2, norm='max', **kwargs)
        assert_allclose(res, exact, rtol=0, atol=epsabs)

        res, err = quad_vec(f, 0, 2, norm='2', **kwargs)
        assert np.linalg.norm(res - exact) < epsabs

        res, err = quad_vec(f, 0, 2, norm='max', points=(0.5, 1.0), **kwargs)
        assert_allclose(res, exact, rtol=0, atol=epsabs)

        res = quad_vec(f, 0, 2, norm='max',
                                   epsrel=1e-8,
                                   limit=10000,
                                   **kwargs)
        assert_allclose(res.integral, exact, rtol=0, atol=epsabs)


@quadrature_params
def test_quad_vec_simple_inf(quadrature):
    def f(x):
        return 1 / (1 + np.float64(x) ** 2)

    for epsabs in [0.1, 1e-3, 1e-6]:
        if quadrature == 'trapezoid' and epsabs < 1e-4:
            # slow: skip
            continue

        kwargs = dict(norm='max', epsabs=epsabs, quadrature=quadrature)

        res, err = quad_vec(f, 0, np.inf, **kwargs)
        assert_allclose(res, np.pi/2, rtol=0, atol=max(epsabs, err))

        res, err = quad_vec(f, 0, -np.inf, **kwargs)
        assert_allclose(res, -np.pi/2, rtol=0, atol=max(epsabs, err))

        res, err = quad_vec(f, -np.inf, 0, **kwargs)
        assert_allclose(res, np.pi/2, rtol=0, atol=max(epsabs, err))

        res, err = quad_vec(f, np.inf, 0, **kwargs)
        assert_allclose(res, -np.pi/2, rtol=0, atol=max(epsabs, err))

        res, err = quad_vec(f, -np.inf, np.inf, **kwargs)
        assert_allclose(res, np.pi, rtol=0, atol=max(epsabs, err))

        res, err = quad_vec(f, np.inf, -np.inf, **kwargs)
        assert_allclose(res, -np.pi, rtol=0, atol=max(epsabs, err))

        res, err = quad_vec(f, np.inf, np.inf, **kwargs)
        assert_allclose(res, 0, rtol=0, atol=max(epsabs, err))

        res, err = quad_vec(f, -np.inf, -np.inf, **kwargs)
        assert_allclose(res, 0, rtol=0, atol=max(epsabs, err))

        res, err = quad_vec(f, 0, np.inf, points=(1.0, 2.0), **kwargs)
        assert_allclose(res, np.pi/2, rtol=0, atol=max(epsabs, err))

    def f(x):
        return np.sin(x + 2) / (1 + x ** 2)
    exact = np.pi / np.e * np.sin(2)
    epsabs = 1e-5

    res = quad_vec(f, -np.inf, np.inf, limit=1000, norm='max', epsabs=epsabs,
                   quadrature=quadrature)
    assert res.status == 1
    assert_allclose(res.integral, exact, rtol=0, atol=max(epsabs, 1.5 * res.abserr))


def test_quad_vec_args():
    def f(x, a):
        return x * (x + a) * np.arange(3)
    a = 2
    exact = np.array([0, 4/3, 8/3])

    res, err = quad_vec(f, 0, 1, args=(a,))
    assert_allclose(res, exact, rtol=0, atol=1e-4)


def _lorenzian(x):
    return 1 / (1 + x**2)


def test_quad_vec_pool():
    f = _lorenzian
    res, err = quad_vec(f, -np.inf, np.inf, norm='max', epsabs=1e-4, workers=4)
    assert_allclose(res, np.pi, rtol=0, atol=1e-4)

    with Pool(10) as pool:
        def f(x):
            return 1 / (1 + x ** 2)
        res, err = quad_vec(f, -np.inf, np.inf, norm='max', epsabs=1e-4, workers=pool.map)
        assert_allclose(res, np.pi, rtol=0, atol=1e-4)


def _func_with_args(x, a):
    return x * (x + a) * np.arange(3)


@pytest.mark.parametrize('extra_args', [2, (2,)])
@pytest.mark.parametrize('workers', [1, 10])
def test_quad_vec_pool_args(extra_args, workers):
    f = _func_with_args
    exact = np.array([0, 4/3, 8/3])

    res, err = quad_vec(f, 0, 1, args=extra_args, workers=workers)
    assert_allclose(res, exact, rtol=0, atol=1e-4)

    with Pool(workers) as pool:
        res, err = quad_vec(f, 0, 1, args=extra_args, workers=pool.map)
        assert_allclose(res, exact, rtol=0, atol=1e-4)


@quadrature_params
def test_num_eval(quadrature):
    def f(x):
        count[0] += 1
        return x**5

    count = [0]
    res = quad_vec(f, 0, 1, norm='max', quadrature=quadrature)
    assert res.neval == count[0]


def test_info():
    def f(x):
        return np.ones((3, 2, 1))

    res = quad_vec(f, 0, 1, norm='max')

    assert res.success is True
    assert res.status == 0
    assert res.message == 'Target precision reached.'
    assert res.neval > 0
    assert res.intervals.shape[1] == 2
    assert res.integrals.shape == (res.intervals.shape[0], 3, 2, 1)
    assert res.errors.shape == (res.intervals.shape[0],)


def test_nan_inf():
    def f_nan(x):
        return np.nan

    def f_inf(x):
        return np.inf if x < 0.1 else 1/x

    res = quad_vec(f_nan, 0, 1)
    assert res.status == 3

    res = quad_vec(f_inf, 0, 1)
    assert res.status == 3


@pytest.mark.parametrize('a,b', [(0, 1), (0, np.inf), (np.inf, 0),
                                 (-np.inf, np.inf), (np.inf, -np.inf)])
def test_points(a, b):
    # Check that initial interval splitting is done according to
    # `points`, by checking that consecutive sets of 15 point (for
    # gk15) function evaluations lie between `points`

    points = (0, 0.25, 0.5, 0.75, 1.0)
    points += tuple(-x for x in points)

    quadrature_points = 15
    interval_sets = []
    count = 0

    def f(x):
        nonlocal count

        if count % quadrature_points == 0:
            interval_sets.append(set())

        count += 1
        interval_sets[-1].add(float(x))
        return 0.0

    quad_vec(f, a, b, points=points, quadrature='gk15', limit=0)

    # Check that all point sets lie in a single `points` interval
    for p in interval_sets:
        j = np.searchsorted(sorted(points), tuple(p))
        assert np.all(j == j[0])

def test_result_object():
        # Check that result object contains attributes 'success', 'status',
        # 'neval', 'intervals', 'integrals', 'errors', 'message'.
        # During the `full_output` deprecation period, also check
        # that specifying `full_output` produces a warning and that values
        # are the same whether `full_output` is True, False, or unspecified.
        def func(x):
            return x**2 + 1

        res = quad_vec(func, 0, 4)
        with np.testing.assert_warns(DeprecationWarning,
                                     match="'full output"):
            res2 = quad_vec(func, 0, 4, full_output=False)
        with np.testing.assert_warns(DeprecationWarning,
                                     match="'full output"):
            res3 = quad_vec(func, 0, 4, full_output=True)

        assert_equal(res, res2)
        assert res.success == res2.success
        assert res.status == res2.status
        assert res.neval == res2.neval
        assert_equal(res.intervals, res2.intervals)
        assert_equal(res.integrals, res2.integrals)
        assert_equal(res.errors, res2.errors)
        assert_equal(res.message, res2.message)

        assert_equal(res, res3[:2])
        assert_equal(res.success, res3[2].success)
        assert_equal(res.status, res3[2].status)
        assert_equal(res.neval, res3[2].neval)
        assert_equal(res.intervals, res3[2].intervals)
        assert_equal(res.integrals, res3[2].integrals)
        assert_equal(res.errors, res3[2].errors)
        assert_equal(res.message, res3[2].message)
