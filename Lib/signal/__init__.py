#
# signal - Signal Processing Tools
#

import os as _os
execfile(_os.path.join(__path__[0],'pre___init__.py'),globals(),locals())

import sigtools
from waveforms import *
from signaltools import *
from bsplines import *
from filter_design import *
from ltisys import *

#---- testing ----#
def test(level=10):
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite(level=level))
    return runner

def test_suite(level=1):
    import scipy_test.testing
    import scipy.signal
    this_mod = scipy.signal
    return scipy_test.testing.harvest_test_suites(this_mod,level=level)

