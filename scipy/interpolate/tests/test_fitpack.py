#!/usr/bin/env python
# Created by Pearu Peterson, June 2003
""" Test functions for interpolate.fitpack2 module
"""
__usage__ = """
Build interpolate:
  python setup_interpolate.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.interpolate.test(<level>)'
Run tests if interpolate is not installed:
  python tests/test_fitpack.py [<level>]
"""
#import libwadpy

import sys
from numpy.testing import *
from numpy import array
set_package_path()
from interpolate.fitpack2 import UnivariateSpline,LSQUnivariateSpline,\
     InterpolatedUnivariateSpline
from interpolate.fitpack2 import LSQBivariateSpline, SmoothBivariateSpline,\
     RectBivariateSpline
restore_path()

class test_UnivariateSpline(NumpyTestCase):
    def check_linear_constant(self):
        x = [1,2,3]
        y = [3,3,3]
        lut = UnivariateSpline(x,y,k=1)
        assert_array_almost_equal(lut.get_knots(),[1,3])
        assert_array_almost_equal(lut.get_coeffs(),[3,3])
        assert_almost_equal(lut.get_residual(),0.0)
        assert_array_almost_equal(lut([1,1.5,2]),[3,3,3])

    def check_linear_1d(self):
        x = [1,2,3]
        y = [0,2,4]
        lut = UnivariateSpline(x,y,k=1)
        assert_array_almost_equal(lut.get_knots(),[1,3])
        assert_array_almost_equal(lut.get_coeffs(),[0,4])
        assert_almost_equal(lut.get_residual(),0.0)
        assert_array_almost_equal(lut([1,1.5,2]),[0,1,2])

class test_LSQBivariateSpline(NumpyTestCase):
    def check_linear_constant(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [3,3,3,3,3,3,3,3,3]
        s = 0.1
        tx = [1+s,3-s]
        ty = [1+s,3-s]
        lut = LSQBivariateSpline(x,y,z,tx,ty,kx=1,ky=1)
        #print lut.get_knots()
        #print lut.get_coeffs()
        #print lut.get_residual()

class test_SmoothBivariateSpline(NumpyTestCase):
    def check_linear_constant(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [3,3,3,3,3,3,3,3,3]
        lut = SmoothBivariateSpline(x,y,z,kx=1,ky=1)
        assert_array_almost_equal(lut.get_knots(),([1,1,3,3],[1,1,3,3]))
        assert_array_almost_equal(lut.get_coeffs(),[3,3,3,3])
        assert_almost_equal(lut.get_residual(),0.0)
        assert_array_almost_equal(lut([1,1.5,2],[1,1.5]),[[3,3],[3,3],[3,3]])

    def check_linear_1d(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [0,0,0,2,2,2,4,4,4]
        lut = SmoothBivariateSpline(x,y,z,kx=1,ky=1)
        assert_array_almost_equal(lut.get_knots(),([1,1,3,3],[1,1,3,3]))
        assert_array_almost_equal(lut.get_coeffs(),[0,0,4,4])
        assert_almost_equal(lut.get_residual(),0.0)
        assert_array_almost_equal(lut([1,1.5,2],[1,1.5]),[[0,0],[1,1],[2,2]])

class test_RectBivariateSpline(NumpyTestCase):
    def check_defaults(self):
        x = array([1,2,3,4,5])
        y = array([1,2,3,4,5])
        z = array([[1,2,1,2,1],[1,2,1,2,1],[1,2,3,2,1],[1,2,2,2,1],[1,2,1,2,1]])
        lut = RectBivariateSpline(x,y,z)
        assert_array_almost_equal(lut(x,y),z)

if __name__ == "__main__":
    NumpyTest().run()
