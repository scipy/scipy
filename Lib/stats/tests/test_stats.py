"""* Test functions for stats module

*"""

from Numeric import *
import Numeric
from scipy_base.fastumath import *
import unittest
from scipy_base.testing import assert_array_equal, assert_equal
from scipy_base.testing import assert_almost_equal
from scipy_base.testing import assert_array_almost_equal
import scipy
import scipy.stats as stats
TestCase = scipy_base.testing.ScipyTestCase


##################################################
### Test for sum

class test_gmean(TestCase):

    def check_1D_list(self):
        a = (1,2,3,4)
        actual= stats.gmean(a)
        desired = power(1*2*3*4,1./4.)
        assert_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.gmean(a,axis=-1)
        assert_almost_equal(desired1,actual,decimal=14)
    def check_1D_array(self):
        a = array((1,2,3,4),Float32)
        actual= stats.gmean(a)
        desired = power(1*2*3*4,1./4.)
        assert_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.gmean(a,axis=-1)
        assert_almost_equal(desired1,actual,decimal=14)

    def check_2D_array_default(self):
        a = array(((1,2,3,4),
                   (1,2,3,4),
                   (1,2,3,4)))
        actual= stats.gmean(a)
        v = power(1*2*3*4,1./4.)
        desired = array((v,v,v))
        assert_array_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.gmean(a,axis=-1)
        assert_array_almost_equal(desired1,actual,decimal=14)

    def check_2D_array_dim0(self):
        a = array(((1,2,3,4),
                   (1,2,3,4),
                   (1,2,3,4)))
        actual= stats.gmean(a,axis=0)
        desired = array((1,2,3,4))
        assert_array_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.gmean(a,axis=0)
        assert_array_almost_equal(desired1,actual,decimal=14)

class test_hmean(TestCase):
    def check_1D_list(self):
        a = (1,2,3,4)
        actual= stats.hmean(a)
        desired =  4. / (1./1 + 1./2 + 1./3 + 1./4)
        assert_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.hmean(array(a),axis=-1)
        assert_almost_equal(desired1,actual,decimal=14)
    def check_1D_array(self):
        a = array((1,2,3,4),Float32)
        actual= stats.hmean(a)
        desired =  4. / (1./1 + 1./2 + 1./3 + 1./4)
        assert_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.hmean(a,axis=-1)
        assert_almost_equal(desired1,actual,decimal=14)

    def check_2D_array_default(self):
        a = array(((1,2,3,4),
                   (1,2,3,4),
                   (1,2,3,4)))
        actual= stats.hmean(a)
        v = 4. / (1./1 + 1./2 + 1./3 + 1./4)
        desired = array((v,v,v))      
        assert_array_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.hmean(a,axis=-1)
        assert_array_almost_equal(desired1,actual,decimal=14)

    def check_2D_array_dim0(self):
        a = array(((1,2,3,4),
                   (1,2,3,4),
                   (1,2,3,4)))
        actual= stats.hmean(a,axis=0)
        desired = array((1.,2.,3.,4.))        
        assert_array_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.hmean(a,axis=0)
        assert_array_almost_equal(desired1,actual,decimal=14)

class test_mean(TestCase):
    def check_basic(self):
        a = [3,4,5,10,-3,-5,6]
        af = [3.,4,5,10,-3,-5,-6]
        Na = len(a)
        Naf = len(af)
        mn1 = 0.0
        for el in a:
            mn1 += el / float(Na)
        assert_almost_equal(stats.mean(a),mn1,11) 
        mn2 = 0.0
        for el in af:
            mn2 += el / float(Naf)
        assert_almost_equal(stats.mean(af),mn2,11)

    def check_2d(self):
        a = [[1.0, 2.0, 3.0],
             [2.0, 4.0, 6.0],
             [8.0, 12.0, 7.0]]
        A = array(a,'d')
        N1,N2 = (3,3)
        mn1 = zeros(N2,'d')
        for k in range(N1):
            mn1 += A[k,:] / N1
        Numeric.allclose(stats.mean(a),mn1,rtol=1e-13,atol=1e-13)
        mn2 = zeros(N1,'d')
        for k in range(N2):
            mn2 += A[:,k] / N2
        Numeric.allclose(stats.mean(a,axis=0),mn2,rtol=1e-13,atol=1e-13)            

class test_median(TestCase):
    def check_basic(self):
        a1 = [3,4,5,10,-3,-5,6]
        a2 = [3,-6,-2,8,7,4,2,1]
        a3 = [3.,4,5,10,-3,-5,-6,7.0]
        assert_equal(stats.median(a1),4)
        assert_equal(stats.median(a2),2.5)
        assert_equal(stats.median(a3),3.5)        

class test_std(TestCase):
    def check_basic(self):
        a = [3,4,5,10,-3,-5,6]
        af = [3.,4,5,10,-3,-5,-6]
        Na = len(a)
        Naf = len(af)
        assert_almost_equal(stats.std(a),5.2098807225172772,11)
        assert_almost_equal(stats.std(af),5.9281411203561225,11)

    def check_2d(self):
        a = [[1.0, 2.0, 3.0],
             [2.0, 4.0, 6.0],
             [8.0, 12.0, 7.0]]
        b1 = array((3.7859388972001824, 5.2915026221291814,
                    2.0816659994661335))
        b2 = array((4.5276925690687087, 2.1213203435596424,
                    7.2111025509279782))
        assert_array_almost_equal(stats.std(a),b2,11)
        assert_array_almost_equal(stats.std(a,axis=0),b1,11)


class test_cmedian(TestCase):
    def check_basic(self):
        

def test_suite(level=1):
    suites = []
    if level > 0:
        suites.append( unittest.makeSuite(test_gmean,'check_') )
        suites.append( unittest.makeSuite(test_hmean,'check_') )
        suites.append( unittest.makeSuite(test_mean,'check_') )
        suites.append( unittest.makeSuite(test_median,'check_') )
        suites.append( unittest.makeSuite(test_std,'check_') )
    total_suite = unittest.TestSuite(suites)
    return total_suite

def test(level=10):
    all_tests = test_suite(level=level)
    runner = unittest.TextTestRunner()
    runner.run(all_tests)
    return runner


if __name__ == "__main__":
    test()
