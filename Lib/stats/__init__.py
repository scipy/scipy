""" Statistical Functions
"""

from stats import *
from rv import *
from pstat import *

#---- testing ----#
def test(level=10):
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite(level=level))
    return runner

def test_suite(level=1):
    import scipy_test
    import scipy.stats
    this_mod = scipy.stats
    return scipy_test.harvest_test_suites(this_mod,level=level)
