from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (run_module_suite, assert_equal, assert_almost_equal,
                           assert_allclose)
from scipy.special import boxcox, boxcox1p


# mpmath functions used to generate the expected values:
"""
def mp_boxcox(x, lmbda, dps=200):
    with mp.workdps(dps):
        x = mp.mpf(x)
        lmbda = mp.mpf(lmbda)
        if lmbda == 0:
            return mp.log(x)
        else:
            return mp.powm1(x, lmbda) / lmbda


def mp_boxcox1p(x, lmbda, dps=200):
    with mp.workdps(dps):
        x = mp.mpf(x)
        lmbda = mp.mpf(lmbda)
        one = mp.mpf(1)
        if lmbda == 0:
            return mp.log(one + x)
        else:
            return mp.powm1(one + x, lmbda) / lmbda
"""


def test_boxcox_basic():
    x = np.array([1,2,3])
    y = boxcox(x, 0)
    yield assert_almost_equal, y, np.log(x)
    y = boxcox(x, 1)
    yield assert_almost_equal, y, x - 1
    y = boxcox(x, 2)
    yield assert_almost_equal, y, 0.5*(x**2 - 1)
    lam = np.array([0.5, 1, 2])

    y = boxcox(0, lam)
    yield assert_almost_equal, y, -1.0 / lam


def test_boxcox_nonfinite():
    x = np.array([-1, -1, -0.5])
    y = boxcox(x, [0.5, 2.0, -1.5])
    yield assert_equal, y, np.array([np.nan, np.nan, np.nan])
    x = 0
    y = boxcox(x, [-2.5, 0])
    yield assert_equal, y, np.array([-np.inf, -np.inf])


def test_boxcox_vs_mp():
    lams = np.array([
        -10.0,
        -1.0,
        -1.0000000000000001e-05,
        -1e-10,
        0.0,
        1e-10,
        1.0000000000000001e-05,
        1.0,
        10.0,
    ])
    xs = np.array([
        1.0000000000000001e-30,
        1e-08,
        1.0,
        100.0,
    ])
    # `expected` was computed with high precision using mpmath.
    expected = np.array([
        [-9.9999999999999931e+298,  -9.999999999999998e+78,
                              0.0,    0.099999999999999992,],
        [                  -1e+30,             -99999999.0,
                              0.0,     0.98999999999999999,],
        [     -69.101416825899577,     -18.422377455528061,
                              0.0,      4.6050641496536056,],
        [     -69.077553028406797,      -18.42068076091844,
                              0.0,      4.6051701849277116,],
        [     -69.077552789821382,     -18.420680743952367,
                              0.0,      4.6051701859880909,],
        [     -69.077552551235968,     -18.420680726986294,
                              0.0,      4.6051701870484703,],
        [     -69.053699741007833,      -18.41898424072776,
                              0.0,      4.6052762255780619,],
        [                    -1.0,    -0.99999998999999995,
                              0.0,                    99.0,],
        [    -0.10000000000000001,    -0.10000000000000001,
                              0.0,   9.999999999999998e+18,],
    ])
    b = boxcox(xs, lams.reshape(-1, 1))
    assert_allclose(b, expected, rtol=1e-13)


def test_boxcox1p_basic():
    x = np.array([-0.25, -1e-20, 0, 1e-20, 0.25, 1, 3])
    y = boxcox1p(x, 0)
    yield assert_almost_equal, y, np.log1p(x)

    y = boxcox1p(x, 1)
    yield assert_almost_equal, y, x

    y = boxcox1p(x, 2)
    yield assert_almost_equal, y, 0.5*x*(2 + x)

    lam = np.array([0.5, 1, 2])
    y = boxcox1p(-1, lam)
    yield assert_almost_equal, y, -1.0 / lam


def test_boxcox1p_nonfinite():
    x = np.array([-2, -2, -1.5])
    y = boxcox1p(x, [0.5, 2.0, -1.5])
    yield assert_equal, y, np.array([np.nan, np.nan, np.nan])
    x = -1
    y = boxcox1p(x, [-2.5, 0])
    yield assert_equal, y, np.array([-np.inf, -np.inf])


def test_boxcox1p_vs_mp():
    lams = np.array([
        -10.0,
        -1.0,
        -9.9999999999999995e-07,
        9.9999999999999995e-07,
        1.0,
        10.0,
    ])
    xs = np.array([
        -0.10000000000000001,
        -1e-08,
        -9.9999999999999995e-21,
        9.9999999999999995e-21,
        1e-08,
        0.10000000000000001,
    ])
    # `expected` was computed with high precision using mpmath.
    expected = np.array([
        [    -0.18679719907924416, -1.0000000550000023e-08,
          -1.0000000000000001e-20,  9.9999999999999979e-21,
           9.9999994500000207e-09,    0.061445671057046826,],
        [    -0.11111111111111112, -1.0000000100000002e-08,
          -1.0000000000000001e-20,  9.9999999999999979e-21,
           9.9999999000000002e-09,    0.090909090909090912,],
        [    -0.10536052120824564, -1.0000000050000051e-08,
          -1.0000000000000001e-20,  9.9999999999999979e-21,
           9.9999999499999497e-09,     0.09531017526230981,],
        [    -0.10536051010740738, -1.0000000049999951e-08,
          -1.0000000000000001e-20,  9.9999999999999979e-21,
            9.999999950000049e-09,    0.095310184346340185,],
        [    -0.10000000000000001,                  -1e-08,
          -9.9999999999999995e-21,  9.9999999999999995e-21,
                            1e-08,     0.10000000000000001,],
        [    -0.06513215599000001, -9.9999995500000125e-09,
          -9.9999999999999995e-21,  9.9999999999999995e-21,
           1.0000000450000011e-08,           0.15937424601,],
    ])

    b = boxcox1p(xs, lams.reshape(-1, 1))
    assert_allclose(b, expected, rtol=1e-13)


if __name__ == '__main__':
    run_module_suite()
