#!/usr/bin/env python
# Created by John Travers, Robert Hetland, 2007
""" Test functions for rbf module """

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_almost_equal
from numpy import linspace, sin, random, exp, log10, allclose
from scipy.interpolate.rbf import Rbf

FUNCTIONS = ('multiquadric', 'inverse multiquadric', 'gaussian',
             'cubic', 'quintic', 'thin-plate', 'linear')

def check_rbf1d_interpolation(function):
    """Check that the Rbf function interpolates throught the nodes (1D)"""
    x = linspace(0,10,9)
    y = sin(x)
    rbf = Rbf(x, y, function=function)
    yi = rbf(x)
    assert_array_almost_equal(y, yi)
    assert_almost_equal(rbf(float(x[0])), y[0])

def check_rbf2d_interpolation(function):
    """Check that the Rbf function interpolates throught the nodes (2D)"""
    x = random.rand(50,1)*4-2
    y = random.rand(50,1)*4-2
    z = x*exp(-x**2-1j*y**2)
    rbf = Rbf(x, y, z, epsilon=2, function=function)
    zi = rbf(x, y)
    zi.shape = x.shape
    assert_array_almost_equal(z, zi)

def check_rbf3d_interpolation(function):
    """Check that the Rbf function interpolates throught the nodes (3D)"""
    x = random.rand(50,1)*4-2
    y = random.rand(50,1)*4-2
    z = random.rand(50,1)*4-2
    d = x*exp(-x**2-y**2)
    rbf = Rbf(x, y, z, d, epsilon=2, function=function)
    di = rbf(x, y, z)
    di.shape = x.shape
    assert_array_almost_equal(di, d)

def test_rbf_interpolation():
    for function in FUNCTIONS:
        yield check_rbf1d_interpolation, function
        yield check_rbf2d_interpolation, function
        yield check_rbf3d_interpolation, function

def check_rbf1d_regularity(function, atol):
    """Check that the Rbf function approximates a smooth function well away
    from the nodes."""
    x = linspace(0, 10, 9)
    y = sin(x)
    rbf = Rbf(x, y, function=function)
    xi = linspace(0, 10, 100)
    yi = rbf(xi)
    #import matplotlib.pyplot as plt
    #plt.figure()
    #plt.plot(x, y, 'o', xi, sin(xi), ':', xi, yi, '-')
    #plt.title(function)
    #plt.show()
    msg = "abs-diff: %f" % abs(yi - sin(xi)).max()
    assert allclose(yi, sin(xi), atol=atol), msg

def test_rbf_regularity():
    tolerances = {
        'multiquadric': 0.05,
        'inverse multiquadric': 0.02,
        'gaussian': 0.01,
        'cubic': 0.15,
        'quintic': 0.1,
        'thin-plate': 0.1,
        'linear': 0.2
    }
    for function in FUNCTIONS:
        yield check_rbf1d_regularity, function, tolerances.get(function, 1e-2)

def test_default_construction():
    """Check that the Rbf class can be constructed with the default
    multiquadric basis function. Regression test for ticket #1228."""
    x = linspace(0,10,9)
    y = sin(x)
    rbf = Rbf(x, y)
    yi = rbf(x)
    assert_array_almost_equal(y, yi)


def test_function_is_callable():
    """Check that the Rbf class can be constructed with function=callable."""
    x = linspace(0,10,9)
    y = sin(x)
    linfunc = lambda x:x
    rbf = Rbf(x, y, function=linfunc)
    yi = rbf(x)
    assert_array_almost_equal(y, yi)
