""" Test functions for basic module

"""

import unittest
import scipy_base.limits as limits
from scipy_test.testing import assert_array_equal, assert_equal
from scipy_test.testing import assert_almost_equal, assert_array_almost_equal
from scipy import *

##################################################


class test_roots(unittest.TestCase):        
    def check_basic(self):
        a1 = [1,-4,4]
        a2 = [4,-16,16]
        a3 = [1,5,6]
        assert_array_almost_equal(roots(a1),[2,2],11)
        assert_array_almost_equal(roots(a2),[2,2],11)
        assert_array_almost_equal(sort(roots(a3)),sort([-3,-2]),7)

    def check_inverse(self):
        a = rand(5)
        assert_array_almost_equal(sort(roots(poly(a))),sort(a),5)

class test_factorial(unittest.TestCase):
    def check_basic(self):
        for k in range(0,13):
            assert_equal(factorial(k),product(grid[1:k+1]))
        self.failUnlessRaises(ValueError, factorial, -10)        

    def check_exact(self):
        resdict = {4:24L,10:3628800L,15:1307674368000L,
                   19:121645100408832000L,
                   100:93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000L}
        for key in resdict.keys():
            assert_equal(factorial(key,exact=1),resdict[key])

class test_comb(unittest.TestCase):
    def check_basic(self):
        for N in range(0,11):
            for k in range(0,N+1):
                ans = product(grid[N-k+1:N+1]) / product(grid[1:k+1])
                assert_almost_equal(comb(N,k),ans,9)
        self.failUnlessRaises(ValueError, comb, -10,1)
        self.failUnlessRaises(ValueError, comb, 10,-1)
        self.failUnlessRaises(ValueError, comb, -10,-3)
        self.failUnlessRaises(ValueError, comb, 10,11)

    def check_exact(self):
        resdict = {(10,2):45L, (10,5):252L,
                   (1000,20):339482811302457603895512614793686020778700L,
                   (1000,975):47641862536236518640933948075167736642053976275040L
                   }
        for key in resdict.keys():
            assert_equal(comb(key[0],key[1],exact=1),resdict[key])

                                        
##################################################

def test_suite(level=1):
    suites = []
    if level > 0:
        suites.append( unittest.makeSuite(test_roots,'check_') )
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
