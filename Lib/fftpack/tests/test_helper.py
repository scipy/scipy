#!/usr/bin/env python
# Created by Pearu Peterson, September 2002
""" Test functions for fftpack.helper module
"""
__usage__ = """
Build fftpack:
  python setup_fftpack.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.fftpack.test(<level>)'
Run tests if fftpack is not installed:
  python tests/test_helper.py [<level>]
"""

import sys
from scipy_test.testing import set_package_path
set_package_path()
from fftpack import fftshift,ifftshift,fftfreq,rfftfreq
del sys.path[0]


import Numeric
from Numeric import pi

from scipy_test.testing import rand
def random(size):
    return rand(*size)

from scipy_test.testing import assert_array_almost_equal, assert_equal
from scipy_test.testing import assert_almost_equal, assert_array_equal
from scipy_test.testing import ScipyTestCase
import unittest

class test_fftshift(ScipyTestCase):

    def check_definition(self):
        x = [0,1,2,3,4,-4,-3,-2,-1]
        y = [-4,-3,-2,-1,0,1,2,3,4]
        assert_array_almost_equal(fftshift(x),y)
        assert_array_almost_equal(ifftshift(y),x)
        x = [0,1,2,3,4,-5,-4,-3,-2,-1]
        y = [-5,-4,-3,-2,-1,0,1,2,3,4]
        assert_array_almost_equal(fftshift(x),y)
        assert_array_almost_equal(ifftshift(y),x)

    def check_inverse(self):
        for n in [1,4,9,100,211]:
            x = random((n,))
            assert_array_almost_equal(ifftshift(fftshift(x)),x)

class test_fftfreq(ScipyTestCase):

    def check_definition(self):
        x = [0,1,2,3,4,-4,-3,-2,-1]
        assert_array_almost_equal(9*fftfreq(9),x)
        assert_array_almost_equal(9*pi*fftfreq(9,pi),x)
        x = [0,1,2,3,4,-5,-4,-3,-2,-1]
        assert_array_almost_equal(10*fftfreq(10),x)
        assert_array_almost_equal(10*pi*fftfreq(10,pi),x)

class test_rfftfreq(ScipyTestCase):

    def check_definition(self):
        x = [0,1,1,2,2,3,3,4,4]
        assert_array_almost_equal(9*rfftfreq(9),x)
        assert_array_almost_equal(9*pi*rfftfreq(9,pi),x)
        x = [0,1,1,2,2,3,3,4,4,5]
        assert_array_almost_equal(10*rfftfreq(10),x)
        assert_array_almost_equal(10*pi*rfftfreq(10,pi),x)

#####################################
def test_suite(level=1):
    suites = []
    if level > 0:
        suites.append( unittest.makeSuite(test_fftshift,'check_') )
        suites.append( unittest.makeSuite(test_fftfreq,'check_') )
        suites.append( unittest.makeSuite(test_rfftfreq,'check_') )
    if level > 5:
        pass
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
