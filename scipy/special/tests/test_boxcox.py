from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import run_module_suite, assert_equal, assert_almost_equal
from scipy.special import boxcox


def test_basic():
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
    x = np.array([-1.0, -0.5])
    y = boxcox(x, np.array([[1],[2]]))
    yield assert_almost_equal, y, np.array([[-2, -1.5], [0, -0.375]])


def test_nonfinite():
    x = np.array([-1, -0.5])
    y = boxcox(x, [0.5, -1.5])
    yield assert_equal, y, np.array([np.nan, np.nan])
    x = 0
    y = boxcox(x, [-2.5, 0])
    yield assert_equal, y, np.array([-np.inf, -np.inf])


if __name__ == '__main__':
    run_module_suite()
