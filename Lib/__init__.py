# modules to import under the scipy namespace
_modules = ["optimize", "integrate", "signal", "special", "io", 
            "interpolate", "stats"]

# namespaces to subsume into the scipy namespace itself
_namespaces = ['MLab', 'misc', 'fastumath'] # MLab includes Numeric
import os,sys
from helpmod import help, source

__all__=[]

for name in _namespaces:
    exec("import %s" % name)
    thelist = eval(name).__dict__.keys()
    exec("del %s" % name) # clean namespace
    exec("from %s import *" % name)
    for key in thelist:
        if key[0] == "_":
            thelist.remove(key)

    __all__.extend(thelist)


for name in _modules:
    exec("import %s" % name)
    __all__.append(name)

__all__.extend(['help', 'source'])


# add some directories to the path so we can import their
# modules.

d,f = os.path.split(__file__)
sys.path.append(os.path.join(d,'gui_thread'))
#import gui_thread

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

# for now, we'll aslo import misc
import misc
    
#---- testing ----#

def test(test_set = 'all'):
    # The standard run test compiler which can take a while.
    # Set test_set = 'fast' to just do the basic tests
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite(test_set))
    return runner

def test_suite(test_set = 'all'):
    import scipy_test
    import scipy
    ignore = ['xplt','plt','gui_thread','linalg','sparse']
    if test_set != 'all':
        ignore += ['compiler']
    return scipy_test.harvest_test_suites(scipy,ignore)

    
