#
# linalg - Linear algebra routines
#

import os as _os
execfile(_os.path.join(__path__[0],'pre___init__.py'),globals(),locals())

from linalg_version import linalg_version as __version__

from basic import *
from decomp import *
from matfuncs import *
from blas import *

################## test functions #########################

def test(level=10):
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite(level))
    return runner

def test_suite(level=1):
    import scipy_test.testing
    exec 'import %s as this_mod' % (__name__)
    return scipy_test.testing.harvest_test_suites(this_mod,level=level)
