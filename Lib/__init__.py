"""
SciPy --- A scientific computing package for Python

Available subpackages
======================
"""

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
import Matrix 
mat = Matrix.Matrix

from helpmod import *

try:
    id(help)
except NameError:
    help = info

# Level 1
# these modules will just be imported (not subsumed)
import special, io, linalg, stats, fftpack
from fftpack import fft, fftn, fft2, ifft, ifft2, ifftn, fftshift,\
     ifftshift, fftfreq #, zeropad, cont_ft
from stats import mean, median, std, cov, corrcoef
from special import isinf, isfinite, isnan

# Functions to be subsumed that need Level 0 and Level 1
from common import *

try:
    from pilutil import *
except ImportError:
    pass

# Level 2
import optimize, integrate, signal, special, interpolate, cow
import ga
#import cluster
import weave

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

# removed gplt because it screws up imports on win32
# removed plt cause it screws up unix when gtk isn't available.
_pkgs = ['cluster','cow','fftpack','ga','integrate',
         'interpolate', 'io', 'linalg', 'optimize', 
         'signal', 'sparse', 'special', 'stats', 'weave',
         'xplt'] #,'plt', 'gplt']

if sys.platform != 'win32':
    _pkgs.append('gplt')
    
_pkg_doc = ['Vector Quantization / Kmeans',
            'Cluster of Workstations',
            'FFT algorithms',
            'Genetic Algorithms',
            'Integration',
            'Interpolation',
            'Input/Output',
            'Linear Algebra',
            'Optimization',
            'Signal Processing',
            'Sparse Matrices',
            'Special Functions',
            'Statistics',
            'C/C++ integration',
            'Plotting with X',
	    'Plotting with gnuplot']
	    #'Plotting with wxPython']

def _extend_doc(st,p,pd):
    _lengths = [len(x) for x in p]
    _ml = max(_lengths)
    k = 0
    for name in p:
        try:
            exec "import %s" % name
            st += "%s%s --- %s\n" % (name, ' '*(_ml-_lengths[k]), pd[k])
        except ImportError:
            pass
        k += 1
    return st

__doc__ = _extend_doc(__doc__,_pkgs,_pkg_doc)

            
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
    import scipy_test.testing
    import scipy
    ignore = ['xplt','plt','gplt','gui_thread','sparse','scipy_version']
    suites = [scipy_test.testing.harvest_test_suites(scipy,ignore,level=level)]
    import unittest
    return unittest.TestSuite(suites)
