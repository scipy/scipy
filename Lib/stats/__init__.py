#
# stats - Statistical Functions
#

import os as _os
execfile(_os.path.join(__path__[0],'pre___init__.py'),globals(),locals())

import pstats
from stats import *
from distributions import *
from rv import *
from morestats import *

try:  # use R functions if installed.
    import rpy
    from rfuncs import *
except ImportError:
    pass


#---- testing ----#
def test(level=10):
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite(level=level))
    return runner

def test_suite(level=1):
    import scipy_test.testing
    import scipy.stats
    this_mod = scipy.stats
    return scipy_test.testing.harvest_test_suites(this_mod,level=level)
