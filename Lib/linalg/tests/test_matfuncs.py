#!/usr/bin/env python
#
# Created by: Pearu Peterson, March 2002
#
""" Test functions for linalg.matfuncs module

"""
__usage__ = """
Build linalg:
  python setup_linalg.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.linalg.test(<level>)'
Run tests if linalg is not installed:
  python tests/test_matfuncs.py [<level>]
"""

from Numeric import array, identity

import sys
from scipy_test.testing import set_package_path
set_package_path()
import scipy_base
import linalg
from linalg import signm,logm
del sys.path[0]

from scipy_test.testing import assert_array_almost_equal
from scipy_test.testing import assert_almost_equal
from scipy_test.testing import ScipyTestCase
import unittest

class test_signm(ScipyTestCase):

    def check_nils(self):
        a = array([[ 29.2, -24.2,  69.5,  49.8,   7. ],
                   [ -9.2,   5.2, -18. , -16.8,  -2. ],
                   [-10. ,   6. , -20. , -18. ,  -2. ],
                   [ -9.6,   9.6, -25.5, -15.4,  -2. ],
                   [  9.8,  -4.8,  18. ,  18.2,   2. ]])
        cr = array([[ 11.94933333,-2.24533333,15.31733333,21.65333333,-2.24533333],
                    [ -3.84266667,0.49866667,-4.59066667,-7.18666667,0.49866667],
                    [ -4.08,0.56,-4.92,-7.6 ,0.56],
                    [ -4.03466667,1.04266667,-5.59866667,-7.02666667,1.04266667],
                    [4.15733333,-0.50133333,4.90933333,7.81333333,-0.50133333]])
        r = signm(a)
        assert_array_almost_equal(r,cr)

    def check_defective1(self):
        a = array([[0.0,1,0,0],[1,0,1,0],[0,0,0,1],[0,0,1,0]])
        r = signm(a)
        #XXX: what would be the correct result?

    def check_defective2(self):
        a = array((
            [29.2,-24.2,69.5,49.8,7.0],
            [-9.2,5.2,-18.0,-16.8,-2.0],
            [-10.0,6.0,-20.0,-18.0,-2.0],
            [-9.6,9.6,-25.5,-15.4,-2.0],
            [9.8,-4.8,18.0,18.2,2.0]))
        r = signm(a)
        #XXX: what would be the correct result?

    def check_defective3(self):
        a = array([[ -2.,  25.,   0.,   0.,   0.,   0.,   0.],
                   [  0.,  -3.,  10.,   3.,   3.,   3.,   0.],
                   [  0.,   0.,   2.,  15.,   3.,   3.,   0.],
                   [  0.,   0.,   0.,   0.,  15.,   3.,   0.],
                   [  0.,   0.,   0.,   0.,   3.,  10.,   0.],
                   [  0.,   0.,   0.,   0.,   0.,  -2.,  25.],
                   [  0.,   0.,   0.,   0.,   0.,   0.,  -3.]])
        r = signm(a)
        #XXX: what would be the correct result?

class test_logm(ScipyTestCase):

    def check_nils(self):
        a = array([[ -2.,  25.,   0.,   0.,   0.,   0.,   0.],
                   [  0.,  -3.,  10.,   3.,   3.,   3.,   0.],
                   [  0.,   0.,   2.,  15.,   3.,   3.,   0.],
                   [  0.,   0.,   0.,   0.,  15.,   3.,   0.],
                   [  0.,   0.,   0.,   0.,   3.,  10.,   0.],
                   [  0.,   0.,   0.,   0.,   0.,  -2.,  25.],
                   [  0.,   0.,   0.,   0.,   0.,   0.,  -3.]])
        logm((identity(7)*3.1+0j)-a)

#####################################
def test_suite(level=1):
    suites = []
    if level > 0:
        suites.append( unittest.makeSuite(test_signm,'check_') )
        suites.append( unittest.makeSuite(test_logm,'check_') )
    if level > 5:
        pass

    total_suite = unittest.TestSuite(suites)
    return total_suite

def test(level=10,verbosity=2):
    all_tests = test_suite(level)
    runner = unittest.TextTestRunner(verbosity=verbosity)
    runner.run(all_tests)
    return runner

if __name__ == "__main__":
    if len(sys.argv)>1:
        level = eval(sys.argv[1])
    else:
        level = 1
    test(level)
