try:
    import gui_thread
except:
    print "Warning:  wxPython not loaded"

#get version
import scipy_version
__version__ = scipy_version.scipy_version
del scipy_version

# SciPy has levels
# Level 0 -- Numeric and core routines in basic.py, misc.py, and handy.py
#
# Level 1 -- Level 0 + fft, special, linalg (these can depend on Level 0)
# Level 1a -- Core routines which depend on Level 1.
# Level 2 -- plotting interface.
# Packages which define own functions plus depend on Levels 0-2.
#
# Level 0, 1, 2 should be imported in order and then other levels imported
#   as available.


import Numeric
import os,sys
from helpmod import help, source
from Matrix import Matrix as Mat
import fastumath

Inf = inf = Numeric.array(1e308)**10
NaN = nan = Numeric.array(0.0) / Numeric.array(0.0)

import string
def somenames2all(alist, namedict, gldict):
    for key in namedict.keys():
        exec("from %s import %s" % (key, string.join(namedict[key],',')), gldict)
        alist.extend(namedict[key])
    
def names2all(alist, spaces, gldict):
    for name in spaces:
        exec("import %s" % name, gldict)
        thelist = eval(name,gldict).__dict__.keys()
        exec("from %s import *" % name, gldict)
        exec("del %s" % name, gldict)
        for key in thelist:
            if key[0] == "_":
                thelist.remove(key)
        alist.extend(thelist)

def modules2all(alist, mods, gldict):
    for name in mods:
        exec("import %s" % name, gldict)
        alist.append(name)
    
def objects2all(alist, objlist):
    alist.extend(objlist)

# modules to import under the scipy namespace
_level0 = ["handy", "misc", "basic", "fastumath"]
_partials0 = {'Matrix' : ['Matrix']}
# these modules will just be imported (not subsumed)
_level0_importonly = []

_level1 = ["special", "io", "linalg"]  # fft is in this group.
_level1a = ["basic1a"] # functions to be subsumed into scipy namespace which
                      # require level 0 and level 1
# these modules will just be imported (not subsumed)                      
_level1_importonly = []

_level3 = ["optimize", "integrate", "signal", "special", "interpolate", "stats", "cow", "ga", "compiler", "cluster"]

__all__=[]
somenames2all(__all__, _partials0, globals())
names2all(__all__, _level0, globals())
modules2all(__all__, _level0, globals())
modules2all(__all__, _level0_importonly, globals())
objects2all(__all__, ['help', 'source', "Inf", "inf", "NaN", "nan", "Mat"])

# Level 1
try:
    import scipy.fftw
    __all__.append('fftw')
except ImportError:
    print "Warning: FFT package not found."

_partials1 = {'fftw' : ['fft', 'fftnd', 'fft2d', 'fft3d',
                        'ifft', 'ifft2d', 'ifft3d', 'ifftnd']}
modules2all(__all__, _level1, globals())
somenames2all(__all__, _partials1, globals())

# Level 1a
names2all(__all__, _level1a, globals())
modules2all(__all__, _level1a, globals())
modules2all(__all__, _level1a_importonly, globals())

# Level 2
_plot = []
try:
    import xplt
    __all__.append('xplt')
    _plot.append('xplt')
except ImportError:
    pass

try:
    import gplt
    __all__.append('gplt')
    _plot.append('gplt')
except ImportError:
    pass

if _plot == []:
    print "Warning: No plotting available."
else:
    print "Plotting methods available: ", _plot

# Level 3

modules2all(__all__, _level3, globals())


#---- testing ----#

def test(test_set = 'fast'):
    # The standard run test compiler which can take a while.
    # Set test_set = 'fast' to just do the basic tests
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite(test_set))
    return runner

def test_all():
    test('all')
    
def test_suite(test_set = 'all'):
    import scipy_test
    import scipy
    ignore = ['xplt','plt','gplt','gui_thread','linalg','sparse','scipy_version']
    if test_set != 'all':
        ignore += ['compiler']
    return scipy_test.harvest_test_suites(scipy,ignore)












