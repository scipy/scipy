__all__ = ["optimize", "integrate", "signal", "special", "io", "fft",
           "interpolate", "stats",'fastumath']
           
from Numeric import *
import os,sys
from helpmod import help, source
from misc import *


for name in __all__:
    exec("import %s" % name)
# add some directories to the path so we can import their
# modules.

d,f = os.path.split(__file__)
sys.path.append(os.path.join(d,'gui_thread'))
#import gui_thread

sys.path.append(os.path.join(d,'pyunit-1.3.1'))

try:
    import scipy.xplt
    os.environ['GISTPATH'] = os.path.join('a','usr','local','scipy','xplt')[1:]    
except ImportError:
    pass

import unittest
   

    
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

    
