# Author:  Travis Oliphant, 2002
# 


import unittest
import scipy
from scipy_base.testing import assert_array_equal, assert_equal, assert_approx_equal
from scipy_base.testing import assert_almost_equal, assert_array_almost_equal
import Numeric
N = Numeric

class test_shapiro(unittest.TestCase):
    def check_basic(self):
        x1 = [0.11,7.87,4.61,10.14,7.95,3.14,0.46,
              4.43,0.21,4.75,0.71,1.52,3.24,
              0.93,0.42,4.97,9.53,4.55,0.47,6.66]
        w,pw = scipy.stats.shapiro(x1)
        assert_almost_equal(w,0.90047299861907959,8)
        assert_almost_equal(pw,0.042089745402336121,8)
        x2 = [1.36,1.14,2.92,2.55,1.46,1.06,5.27,-1.11,
              3.48,1.10,0.88,-0.51,1.46,0.52,6.20,1.69,
              0.08,3.67,2.81,3.49]
        w,pw = scipy.stats.shapiro(x2)
        assert_almost_equal(w,0.9590269923210144,8)
        assert_almost_equal(pw,0.52459925413131714,8)
        
        
def test_suite(level=1):
    suites = []
    if level > 0:
        suites.append( unittest.makeSuite(test_shapiro,'check_') )
        
    total_suite = unittest.TestSuite(suites)
    return total_suite

def test(level=10):
    all_tests = test_suite(level=level)
    runner = unittest.TextTestRunner()
    runner.run(all_tests)
    return runner

if __name__ == "__main__":
    test()


        

