#this program corresponds to special.py
import os
#import fpformat
import unittest
from unittest import TestCase
import scipy_base.limits as limits
from scipy_test.testing import assert_array_equal, assert_equal
from scipy_test.testing import assert_almost_equal, assert_array_almost_equal
import scipy.signal as signal

class test_convolve(TestCase):
    def check_basic(self):
        a = [3,4,5,6,5,4]
        b = [1,2,3]
        c = signal.convolve(a,b)
        assert_array_equal(c,array([3,10,22,28,32,32,23,12]))

class test_medfilt(TestCase):
    def check_basic(self):
        f = [[3,4,5],[2,3,4],[1,2,5]]
        d = signal.medfilt(f)
        assert_array_equal(d, [[0,3,0],[2,3,3],[0,2,0]])

class test_wiener(TestCase):
    def check_basic(self):
        g = Numeric.array([[5,6,4,3],[3,5,6,2],[2,3,5,6],[1,6,9,7]],'d')
        correct = Numeric.array([[2.16374269,3.2222222222, 2.8888888889, 1.6666666667],[2.666666667, 4.33333333333, 4.44444444444, 2.8888888888],[2.222222222, 4.4444444444, 5.4444444444, 4.801066874837],[1.33333333333, 3.92735042735, 6.0712560386, 5.0404040404]])
        h = wiener(g)
        assert_array_almost_equal(h,correct,decimal=6)


def test_suite(level=1):
    suites = []
    if level > 0:
        suites.append( unittest.makeSuite(test_convolve, 'check_') )
        suites.append( unittest.makeSuite(test_wiener, 'check_') )
        suites.append( unittest.makeSuite(test_medfilt, 'check_') )
        
    total_suite = unittest.TestSuite(suites)
    return total_suite

def test(level=10):
    all_tests = test_suite(level=level)
    runner = unittest.TextTestRunner()
    runner.run(all_tests)
    return runner


if __name__ == "__main__":
    test()







