#
# optimize - Optimization Tools
#

from pre___init__ import __doc__

from optimize import *
from minpack import *
from zeros import *
from anneal import *

################## test functions #########################

def test(level=10):
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite(level=level))
    return runner

def test_suite(level=1):
    import scipy_test.testing
    import scipy.optimize
    this_mod = scipy.optimize
    return scipy_test.testing.harvest_test_suites(this_mod,level=level)

