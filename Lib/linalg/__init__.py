#
# linalg - Linear algebra routines
#

from pre___init__ import __doc__
from linalg_version import linalg_version as __version__

from basic import *
from decomp import *
from matfuncs import *
from blas import *

################## test functions #########################

def test(level=10,verbosity=2):
    import unittest
    runner = unittest.TextTestRunner(verbosity=verbosity)
    runner.run(test_suite(level))
    return runner

def test_suite(level=1):
    import scipy_test.testing
    exec 'import %s as this_mod' % (__name__)
    return scipy_test.testing.harvest_test_suites(this_mod,level=level)
