import gui_thread
#get version
import scipy_version
__version__ = scipy_version.version
del scipy_version

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
_modules = ["fastumath", "misc", "optimize", "integrate", "signal",
            "special", "io", "interpolate", "stats"]

# namespaces to subsume into the scipy namespace itself
_namespaces = ['MLab','handy', 'misc', 'fastumath'] # MLab includes Numeric

# partial list of namespaces to get
_partials = {'Matrix' : ['Matrix']}
import os,sys
from helpmod import help, source

__all__=[]

somenames2all(__all__, _partials, globals())
names2all(__all__, _namespaces, globals())
modules2all(__all__, _modules, globals())
objects2all(__all__, ['help', 'source'])


try:
    import scipy.fft
    __all__.append('fft')
except ImportError:
    pass

try:
    import xplt
    __all__.append('xplt')
except ImportError:
    pass

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












