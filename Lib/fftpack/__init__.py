#
# fftpack - Discrete Fourier Transform algorithms.
#
# Created: Pearu Peterson, August,September 2002

import os as _os
execfile(_os.path.join(__path__[0],'pre___init__.py'),globals(),locals())

from fftpack_version import fftpack_version as __version__

from basic import *
from pseudo_diffs import *
from helper import *

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
