#!/usr/bin/env python
# Created by John Travers, Robert Hetland, 2007
""" Test functions for rbf module """

from numpy.testing import assert_array_almost_equal, assert_almost_equal
from numpy import linspace, sin, random, exp
from scipy.interpolate.rbf import Rbf

FUNCTIONS = ('multiquadric', 'inverse multiquadric', 'gaussian',
             'cubic', 'quintic', 'thin-plate')

def check_rbf1d(function):
    x = linspace(0,10,9)
    y = sin(x)
    rbf = Rbf(x, y, function=function)
    yi = rbf(x)
    assert_array_almost_equal(y, yi)
    assert_almost_equal(rbf(float(x[0])), y[0])

def check_rbf2d(function):
    x = random.rand(50,1)*4-2
    y = random.rand(50,1)*4-2
    z = x*exp(-x**2-1j*y**2)
    rbf = Rbf(x, y, z, epsilon=2, function=function)
    zi = rbf(x, y)
    zi.shape = x.shape
    assert_array_almost_equal(z, zi)

def check_rbf3d(function):
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
        yield check_rbf1d, function
        yield check_rbf2d, function
        yield check_rbf3d, function
