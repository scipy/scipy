""" Statistical Functions
"""
import pstats
from rv import *
from stats import *
from rv2 import *
from distributions import *

#---- testing ----#
def test(level=10):
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite(level=level))
    return runner

def test_suite(level=1):
    import scipy_base.testing
    import scipy.stats
    this_mod = scipy.stats
    return scipy_base.testing.harvest_test_suites(this_mod,level=level)
