"""* Test functions for misc module

*"""

from Numeric import *
from fastumath import *
import unittest
from scipy_test import assert_array_equal, assert_equal
from scipy_test import assert_almost_equal
from scipy_test import assert_array_almost_equal
import scipy
import scipy.stats stats


##################################################
### Test for sum

class test_gmean(unittest.TestCase):

    def check_1D_list(self):
        a = (1,2,3,4)
        actual= stats.gmean(a)
        desired = power(1*2*3*4,1./4.)
        assert_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.gmean(a,dimension=-1)
        assert_almost_equal(desired1,actual,decimal=14)
    def check_1D_array(self):
        a = array((1,2,3,4),Float32)
        actual= stats.gmean(a)
        desired = power(1*2*3*4,1./4.)
        assert_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.gmean(a,dimension=-1)
        assert_almost_equal(desired1,actual,decimal=14)

    def check_2D_array_default(self):
        a = array(((1,2,3,4),
                   (1,2,3,4),
                   (1,2,3,4)))
        actual= stats.geometric_mean(a)
        v = power(1*2*3*4,1./4.)
        desired = array((v,v,v))
        assert_array_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.ageometricmean(a,dimension=-1)
        assert_array_almost_equal(desired1,actual,decimal=14)

    def check_2D_array_dim0(self):
        a = array(((1,2,3,4),
                   (1,2,3,4),
                   (1,2,3,4)))
        actual= stats.geometric_mean(a,dimension=0)
        desired = array((1,2,3,4))
        assert_array_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.ageometricmean(a,dimension=0)
        assert_array_almost_equal(desired1,actual,decimal=14)

class test_harmonic_mean(unittest.TestCase):
    def check_1D_list(self):
        a = (1,2,3,4)
        actual= stats.harmonic_mean(a)
        desired =  4. / (1./1 + 1./2 + 1./3 + 1./4)
        assert_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.aharmonicmean(array(a),dimension=-1)
        assert_almost_equal(desired1,actual,decimal=14)
    def check_1D_array(self):
        a = array((1,2,3,4),Float32)
        actual= stats.harmonic_mean(a)
        desired =  4. / (1./1 + 1./2 + 1./3 + 1./4)
        assert_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.aharmonicmean(a,dimension=-1)
        assert_almost_equal(desired1,actual,decimal=14)

    def check_2D_array_default(self):
        a = array(((1,2,3,4),
                   (1,2,3,4),
                   (1,2,3,4)))
        actual= stats.harmonic_mean(a)
        v = 4. / (1./1 + 1./2 + 1./3 + 1./4)
        desired = array((v,v,v))      
        assert_array_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.aharmonicmean(a,dimension=-1)
        assert_array_almost_equal(desired1,actual,decimal=14)

    def check_2D_array_dim0(self):
        a = array(((1,2,3,4),
                   (1,2,3,4),
                   (1,2,3,4)))
        actual= stats.harmonic_mean(a,dimension=0)
        desired = array((1.,2.,3.,4.))        
        assert_array_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.aharmonicmean(a,dimension=0)
        assert_array_almost_equal(desired1,actual,decimal=14)

def test_suite(level=1):
    suites = []
    if level > 0:
        suites.append( unittest.makeSuite(test_geometric_mean,'check_') )
        suites.append( unittest.makeSuite(test_harmonic_mean,'check_') )
    total_suite = unittest.TestSuite(suites)
    return total_suite

def test(level=10):
    all_tests = test_suite(level=level)
    runner = unittest.TextTestRunner()
    runner.run(all_tests)
    return runner


if __name__ == "__main__":
    test()
