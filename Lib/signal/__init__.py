#
# signal - Signal Processing Tools
#

from pre___init__ import __doc__

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

