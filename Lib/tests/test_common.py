""" Test functions for basic module

"""

import unittest
import scipy_base.limits as limits
from scipy_test.testing import assert_array_equal, assert_equal
from scipy_test.testing import assert_almost_equal, assert_array_almost_equal
from scipy_base import sqrt, product, add,  ravel, mgrid
from scipy.common import rand,randn,comb,factorial


##################################################

### Test for rand

class test_rand(unittest.TestCase):
    def __init__(self,*args,**kwds):
        unittest.TestCase.__init__(self,*args,**kwds)
        self.z = rand(100,10,20)
        
    def check_shape(self):
        assert_equal(self.z.shape,(100,10,20))

    def check_mean(self):
        mn = add.reduce(ravel(self.z)) / 200 / 100.0
        assert_almost_equal(mn,0.5,1)

    def check_std(self):
        std = add.reduce((ravel(self.z)-0.5)**2) / 200.0 / 100.0
        assert_almost_equal(std,1.0/12.0,2)

class test_randn(unittest.TestCase):
    def __init__(self,*args,**kwds):
        unittest.TestCase.__init__(self,*args,**kwds)
        self.z = randn(100,10,20)
        
    def check_shape(self):
        assert_equal(self.z.shape,(100,10,20))

    def check_mean(self):
        mn = add.reduce(ravel(self.z)) / 200 / 100.0
        assert_almost_equal(mn,0.0,1)

    def check_std(self):
        std = sqrt(add.reduce(ravel(self.z)**2) / 200.0 / 100.0)
        assert_almost_equal(std,1.0,1)

class test_factorial(unittest.TestCase):
    def check_basic(self):
        for k in range(0,13):
            assert_equal(factorial(k),product(mgrid[1:k+1]))
        assert_equal(factorial(-10),0)

    def check_exact(self):
        resdict = {4:24L,10:3628800L,15:1307674368000L,
                   19:121645100408832000L, -10:0L,
                   100:93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000L}
        for key in resdict.keys():
            assert_equal(factorial(key,exact=1),resdict[key])

class test_comb(unittest.TestCase):
    def check_basic(self):
        for N in range(0,11):
            for k in range(0,N+1):
                ans = product(mgrid[N-k+1:N+1]) / product(mgrid[1:k+1])
                assert_almost_equal(comb(N,k),ans,9)
        assert_equal(comb(-10,1),0)
        assert_equal(comb(10,-1),0)
        assert_equal(comb(-10,-3),0)
        assert_equal(comb(10,11),0)

    def check_exact(self):
        resdict = {(10,2):45L, (10,5):252L,
                   (1000,20):339482811302457603895512614793686020778700L,
                   (1000,975):47641862536236518640933948075167736642053976275040L,
                   (-10,1):0L, (10,-1):0L, (-10,-3):0L,(10,11):0L}
        for key in resdict.keys():
            assert_equal(comb(key[0],key[1],exact=1),resdict[key])
        
#-----------------------------------------------------------------------------

def test_suite(level=1):
    suites = []
    if level > 0:
        suites.append( unittest.makeSuite(test_rand,'check_') )
        suites.append( unittest.makeSuite(test_randn,'check_') )
        suites.append( unittest.makeSuite(test_factorial,'check_') )
        suites.append( unittest.makeSuite(test_comb,'check_') )
    total_suite = unittest.TestSuite(suites)
    return total_suite

def test(level=10):
    all_tests = test_suite(level)
    runner = unittest.TextTestRunner()
    runner.run(all_tests)
    return runner


if __name__ == "__main__":
    test()
