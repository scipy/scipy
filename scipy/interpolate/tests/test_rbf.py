#!/usr/bin/env python
#
# Authors
# - John Travers, Robert Hetland (creators) 2007
# - J.J. Green (extra tests, tidying) 2015

""" Test functions for Rbf module """

from __future__ import division, print_function, absolute_import

import numpy as np
import numpy.testing as npt
from scipy.interpolate.rbf import Rbf

# standard and deprecated RBFs

FUNC_STD = ('multiquadric',
            'inverse_multiquadric',
            'gaussian',
            'thin_plate')


def check_rbf1d_interpolation(function):
    # Check that an Rbf instance interpolates through the nodes (1D)
    x = np.linspace(0, 10, 9)
    f = np.sin(x)
    rbf = Rbf(x, f, function=function)
    fi = rbf(x)
    npt.assert_array_almost_equal(f, fi)
    npt.assert_almost_equal(rbf(float(x[0])), f[0])


def test_rbf_interpolation_1d():
    for function in FUNC_STD:
        yield check_rbf1d_interpolation, function


def check_rbf2d_interpolation(function):
    # Check that an Rbf instance interpolates through the nodes (2D).
    x = np.random.rand(50, 1)*4 - 2
    y = np.random.rand(50, 1)*4 - 2
    f = x*np.exp(-x**2 - 1j*y**2)
    rbf = Rbf(x, y, f, function=function)
    fi = rbf(x, y)
    fi.shape = x.shape
    npt.assert_array_almost_equal(f, fi)


def test_rbf_interpolation_2d():
    for function in FUNC_STD:
        yield check_rbf2d_interpolation, function


def check_rbf3d_interpolation(function):
    # Check that an Rbf instance interpolates through the nodes (3D).
    x = np.random.rand(50, 1)*4 - 2
    y = np.random.rand(50, 1)*4 - 2
    z = np.random.rand(50, 1)*4 - 2
    f = x*np.exp(-x**2 - y**2)
    rbf = Rbf(x, y, z, f, function=function)
    fi = rbf(x, y, z)
    fi.shape = x.shape
    npt.assert_array_almost_equal(fi, f)


def test_rbf_interpolation_3d():
    for function in FUNC_STD:
        yield check_rbf3d_interpolation, function


def check_poly_order_raises(function, smooth, order):
    # Check that specified poly_order causes Rbf to raise ValueError
    x = f = np.linspace(0, 9, 10)
    npt.assert_raises(ValueError, Rbf, x, f, smooth=smooth, poly_order=order)


def check_poly_order_accepted(function, smooth, order):
    # Check that we can create an RBF with the specified poly_order
    # and that that is preserved in the RBF instance
    x = f = np.linspace(0, 9, 10)
    rbf = Rbf(x, f, function=function, smooth=smooth, poly_order=order)
    if order is None:
        npt.assert_(rbf.poly_order is None, 'poly_order lost')
    else:
        npt.assert_(rbf.poly_order == order, 'poly_order lost')


def check_poly_order(function, smooth):
    for order in (None, 0, 1):
        check_poly_order_accepted(function, smooth, order)
    check_poly_order_raises(function, smooth, -1)
    check_poly_order_raises(function, smooth, 2)


def test_rbf_poly_order():
    for function in FUNC_STD:
        for smooth in (0, 1e-3):
            yield check_poly_order, function, smooth


def check_shift_invariance(function, smooth):
    # Check that a vertical shift of interpolation points gives
    # rise to the same shift of the interpolation (or of the
    # approximant, when smooth != 0).
    n = 10
    x = np.linspace(0, n-1, n)
    f = np.random.randn(len(x))
    rbf0 = Rbf(x, f, smooth=smooth, function=function)
    rbf1 = Rbf(x, f + 1, smooth=smooth, function=function)
    x_fine = np.linspace(0, n-1, 5*n)
    f0 = rbf0(x_fine) + 1
    f1 = rbf1(x_fine)
    npt.assert_array_almost_equal(f0, f1)


def test_rbf_shift_invariance():
    for function in FUNC_STD:
        for smooth in [0, 1e-3, 1e-1]:
            yield check_shift_invariance, function, smooth


def check_affine_invariance(function, smooth):
    # Check that an affine transform of interpolation points
    # gives rise to the same for the interpolation (or of the
    # approximant, when smooth != 0).
    n = 10
    x = np.linspace(0, n-1, n)
    f = np.random.randn(len(x))
    rbf0 = Rbf(x, f, smooth=smooth, function=function)
    rbf1 = Rbf(x, f + 3*x + 2, smooth=smooth, function=function)
    x_fine = np.linspace(0, n-1, 5*n)
    f0 = rbf0(x_fine) + 3*x_fine + 2
    f1 = rbf1(x_fine)
    npt.assert_array_almost_equal(f0, f1)


def test_rbf_affine_invariance():
    for function in FUNC_STD:
        for smooth in [0, 1e-3, 1e-1]:
            yield check_affine_invariance, function, smooth


def check_positive_smoothing(function):
    # When using a small positive smoothing parameter on a coarsely
    # discretised step, we expect nodes away from the step to be
    # close (but not that close) to the data samples. (This test will
    # fail if the the `sign` in the smoothing setup is reversed)
    n = 10
    x = np.linspace(0, n-1, n)
    f = np.array([-1 if k < n/2 else 1 for k in x], dtype=np.float)
    rbf = Rbf(x, f, smooth=1e-3, function=function)
    offstep = range(0, n//2-1) + range(n//2+1, n)
    x_offstep = x[offstep]
    f0 = rbf(x_offstep)
    f1 = f[offstep]
    msg = "abs-diff: %f" % abs(f0 - f1).max()
    npt.assert_(np.allclose(f0, f1, atol=0.05), msg)


def test_rbf_positive_smoothing():
    for function in FUNC_STD:
        yield check_positive_smoothing, function


def check_rbf1d_regularity(function, atol):
    # Check that an Rbf instance approximates a smooth function well
    # away from the nodes.
    x = np.linspace(0, 10, 9)
    f = np.sin(x)
    rbf = Rbf(x, f, function=function)
    x_fine = np.linspace(0, 10, 100)
    f0 = rbf(x_fine)
    f1 = np.sin(x_fine)
    msg = "abs-diff: %f" % abs(f0 - f1).max()
    npt.assert_(np.allclose(f0, f1, atol=atol), msg)


def test_rbf_regularity():
    tol = {
        'multiquadric': 0.10,
        'inverse_multiquadric': 0.15,
        'gaussian': 0.15,
        'thin_plate': 0.10,
    }
    for function in FUNC_STD:
        yield check_rbf1d_regularity, function, tol.get(function, 1e-2)


def check_rbf1d_stability(function):
    # Check that an Rbf instance with default epsilon is not subject
    # to overshoot.  Regression for issue #4523.
    #
    # Generate some data (fixed random seed hence deterministic)
    np.random.seed(1234)
    x = np.linspace(0, 10, 50)
    f = x + 4.0 * np.random.randn(len(x))
    rbf = Rbf(x, f, function=function)
    x0 = np.linspace(0, 10, 1000)
    f0 = rbf(x0)
    # subtract the linear trend and make sure there no spikes
    npt.assert_(np.abs(f0-x0).max() / np.abs(f-x).max() < 1.1)


def test_rbf_stability():
    for function in FUNC_STD:
        yield check_rbf1d_stability, function


def test_default_construction():
    # Check that an Rbf instance can be constructed with the default
    # multiquadric basis function. Regression test for ticket #1228.
    x = np.linspace(0, 10, 9)
    f = np.sin(x)
    rbf = Rbf(x, f)
    f0 = rbf(x)
    npt.assert_array_almost_equal(f, f0)


def test_one_arg_function_is_callable():
    # Check that an Rbf instance can be constructed from a
    # one-argument callable function
    x = np.linspace(0,10,9)
    f = np.sin(x)
    rbf = Rbf(x, f, function=lambda x: x)
    f0 = rbf(x)
    npt.assert_array_almost_equal(f, f0)


def test_two_arg_function_is_callable():
    # Check that an Rbf instance can be constructed from a
    # two-argument callable function

    def _func(self, r):
        return self.epsilon + r

    x = np.linspace(0, 10, 9)
    f = np.sin(x)
    rbf = Rbf(x, f, function=_func)
    f0 = rbf(x)
    npt.assert_array_almost_equal(f, f0)


def test_rbf_epsilon_none():
    # Check that an Rbf instance can be constructed using
    # a None value of epsilon (and that results in epsilon
    # being chosen non-zero)
    x = np.linspace(0, 10, 9)
    f = np.sin(x)
    rbf = Rbf(x, f, epsilon=None)
    npt.assert_(rbf.epsilon > 0, 'non-positive epsilon')


if __name__ == "__main__":
    npt.run_module_suite()
