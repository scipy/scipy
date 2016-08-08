from __future__ import division

from numpy.testing import (assert_almost_equal, assert_array_less, assert_)

from ..brent import find_minimum

def test_strictly_convex():
    def func(x, s):
        return (x - s)**2 - 0.8

    (x, fx, niters) = find_minimum(lambda x: func(x, 0), -10, 10)
    assert_almost_equal(x, 0)
    assert_almost_equal(fx, -0.8)
    assert_array_less(niters, 6)

    (x, fx, niters) = find_minimum(lambda x: func(x, 0), 0, 10)
    assert_almost_equal(x, 0)
    assert_array_less(niters, 41)

    (x, fx, niters) = find_minimum(lambda x: func(x, 5), -1000, -100, rtol=1e-10)
    assert_almost_equal(x, -100)
    assert_almost_equal(fx, func(-100, 5), decimal=4)
    assert_array_less(niters, 50)

    (x, fx, niters) = find_minimum(lambda x: func(x, 5), -100, -99, rtol=1e-10)
    assert_almost_equal(x, -99)
    assert_almost_equal(fx, func(-99, 5), decimal=4)

    (x, fx, niters) = find_minimum(lambda x: func(x, 5), 99, 100, rtol=1e-10)
    assert_almost_equal(x, 99)
    assert_almost_equal(fx, func(99, 5), decimal=4)
    assert_array_less(niters, 35)

    (x, fx, niters) = find_minimum(lambda x: func(x, 5), -1, -1/2)
    assert_almost_equal(x, -1/2)
    assert_almost_equal(fx, func(-1/2, 5), decimal=4)
    assert_array_less(niters, 34)

    (x, fx, niters) = find_minimum(lambda x: func(x, 5), 6.5, 6.5)
    assert_almost_equal(x, 6.5)
    assert_almost_equal(fx, func(6.5, 5), decimal=4)
    assert_array_less(niters, 1)

    (x, fx, niters) = find_minimum(lambda x: func(x, 5), 6, 7, rtol=1e-9)
    assert_almost_equal(x, 6)
    assert_almost_equal(fx, func(6, 5), decimal=4)
    assert_array_less(niters, 36)

def test_convex():
    def func(x, s):
        return abs(x - s) + 3

    (x, fx, niters) = find_minimum(lambda x: func(x, 0), -10, +10)
    assert_almost_equal(x, 0)
    assert_almost_equal(fx, 3)
    assert_array_less(niters, 36)

    (x, fx, niters) = find_minimum(lambda x: func(x, -9), -10, +10)
    assert_almost_equal(x, -9)
    assert_array_less(niters, 36)

def test_asymptotic():
    def func(x):
        return -3 + 1/x

    (x, _, niters) = find_minimum(func, 1e-6, +10, rtol=1e-9)
    assert_almost_equal(x, 10)
    assert_array_less(niters, 41)

    (x, _, niters) = find_minimum(func, -10, -1e-6)
    assert_almost_equal(x, -1e-06)
    assert_array_less(niters, 41)

def test_same_point():
    def func(x):
        return x**2
    (x, _, niters) = find_minimum(func, 1.2, 1.2)
    assert_array_less(niters, 1)
    assert_almost_equal(x, 1.2)

def test_piecewise_convex():
    def func(x):
        if abs(x) > 2:
            return x**2
        if abs(x) > 1:
            return abs(x)
        return 1

    (x, fx, niters) = find_minimum(func, -10, +10)
    assert_(-1-1e6 <= x <= 1+1e6)
    assert_almost_equal(fx, 1)
    assert_array_less(niters, 44)

    (x, fx, niters) = find_minimum(func, 1/2, +1)
    assert_(-1-1e6 <= x <= 1+1e6)
    assert_almost_equal(fx, 1)
    assert_array_less(niters, 34)

    (x, fx, niters) = find_minimum(func, -1, -1/2)
    assert_(-1-1e6 <= x <= 1+1e6)
    assert_almost_equal(fx, 1)
    assert_array_less(niters, 34)

    (x, fx, niters) = find_minimum(func, -1/2, +1/2)
    assert_(-1-1e6 <= x <= 1+1e6)
    assert_almost_equal(fx, 1)
    assert_array_less(niters, 36)
