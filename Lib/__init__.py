
from Numeric import *
import os,sys

# add some directories to the path so we can import their
# modules.

d,f = os.path.split(__file__)
sys.path.append(os.path.join(d,'gui_thread'))
#import gui_thread

sys.path.append(os.path.join(d,'pyunit-1.3.1'))
import unittest
   

from misc import *

    
#---- testing ----#

def test():
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())
    return runner

def test_suite():
    import scipy_test
    import scipy
    return scipy_test.harvest_test_suites(scipy)

    
