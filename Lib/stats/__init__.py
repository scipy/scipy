""" Statistical Functions
"""

from stats import *
from rv import *
from pstat import *

#---- testing ----#
def test():
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())
    return runner

def test_suite():
    import scipy_test
    import scipy.stats
    this_mod = scipy.stats
    return scipy_test.harvest_test_suites(this_mod)
