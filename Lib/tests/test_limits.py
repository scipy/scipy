"""* Test functions for limits module.

     Currently empty -- not sure how to test these values
     and routines as they are machine dependent.  Suggestions?
*"""

from Numeric import *
import unittest
import scipy.limits



##################################################
### Test for sum

class test_float(unittest.TestCase):
    def check_nothing(self):
        pass

class test_double(unittest.TestCase):
    def check_nothing(self):
        pass

##################################################


def test_suite():
    suites = []
    suites.append( unittest.makeSuite(test_float,'check_') )
    suites.append( unittest.makeSuite(test_double,'check_') )
    
    total_suite = unittest.TestSuite(suites)
    return total_suite

def test():
    all_tests = test_suite()
    runner = unittest.TextTestRunner()
    runner.run(all_tests)
    return runner


if __name__ == "__main__":
    test()