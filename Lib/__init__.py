#try:
    #import gui_thread
#except:
    #print "Warning:  wxPython not loaded"

from scipy_version import scipy_version as __version__

# SciPy has levels
# Level 0 -- Numeric and core routines in scipy_base
#
# Level 1 -- Level 0 + fft, special, linalg, stats
# Level 1a -- Core routines which depend on Level 1.

# Level 2 -- plotting interface.
#
# Level 3
# Packages which define own functions plus depend on Levels 0-2.
#
# Level 0, 1, 2 should be imported in order and then other levels imported
#   as available.


import Numeric
import os,sys
import string

# Level 0
# modules to import under the scipy namespace
from scipy_base import *

from helpermod import *

# Level 1
# these modules will just be imported (not subsumed)
import special, io, linalg, stats, fftpack
from fftpack import fft, fftn, fft2, ifft, ifft2, ifftn
from stats import mean, median, std, cov, corrcoef
from special import isinf, isfinite, isnan

# Functions to be subsumed that need Level 0 and Level 1
from common import *

try:
    import fftw
except ImportError:
    pass  # Not default anymore


# Level 2
import optimize, integrate, signal, special, interpolate, cow, ga, cluster, weave

# Level 3
_plot = []
try:
    if sys.platform != 'win32':
        import xplt
        #__all__.append('xplt')
        _plot.append('xplt')
except ImportError:
    pass

try:
#   gplt on win32 messes up focus and takes up 99%
#   of processor -- works fine on *nix.
    if sys.platform != 'win32':
        import gplt
        #__all__.append('gplt')
        _plot.append('gplt')
except ImportError:
    pass

#if _plot == []:
#    print "Warning: No plotting available."
#else:
#    print "Plotting methods available: ", _plot



#---- testing ----#

def test(level=1):
    """ From this top level, there are possibly many many tests.
        Test only the quick tests by default.
    """
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite(level))
    return runner

def test_all(level=10):
    test(level)
    
def test_suite(level = 1):
    import scipy_base.testing
    import scipy
    ignore = ['xplt','plt','gplt','gui_thread','linalg','sparse','scipy_version']
    return scipy_base.testing.harvest_test_suites(scipy,ignore,level=level)












