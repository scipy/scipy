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
import libwadpy

import sys
from scipy_test.testing import set_package_path
set_package_path()
from interpolate.fitpack2 import LSQBivariateSpline, SmoothBivariateSpline
del sys.path[0]

from scipy_test.testing import assert_array_almost_equal,assert_almost_equal
from scipy_test.testing import ScipyTestCase
import unittest

class test_LSQBivariateSpline(ScipyTestCase):
    def check_linear_constant(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [3,3,3,3,3,3,3,3,3]
        s = 0.1
        tx = [1+s,3-s]
        ty = [1+s,3-s]
        lut = LSQBivariateSpline(x,y,z,tx,ty,kx=1,ky=1)
        print lut.get_knots()
        print lut.get_coeffs()
        print lut.get_residual()

class test_SmoothBivariateSpline(ScipyTestCase):
    def check_linear_constant(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [3,3,3,3,3,3,3,3,3]
        lut = SmoothBivariateSpline(x,y,z,kx=1,ky=1)
        assert_array_almost_equal(lut.get_knots(),([1,1,3,3],[1,1,3,3]))
        assert_array_almost_equal(lut.get_coeffs(),[3,3,3,3])
        assert_almost_equal(lut.get_residual(),0.0)

#####################################

def test_suite(level=1):
    suites = []
    if level > 0:
        suites.append( unittest.makeSuite(test_LSQBivariateSpline,'check_') )
        suites.append( unittest.makeSuite(test_SmoothBivariateSpline,'check_') )

    total_suite = unittest.TestSuite(suites)
    return total_suite

def test(level=10):
    all_tests = test_suite(level)
    runner = unittest.TextTestRunner()
    runner.run(all_tests)
    return runner

if __name__ == "__main__":
    if len(sys.argv)>1:
        level = eval(sys.argv[1])
    else:
        level = 1
    test(level)
