
from Numeric import *
import os,sys

# add some directories to the path so we can import their
# modules.

d,f = os.path.split(__file__)
sys.path.append(os.path.join(d,'gui_thread'))
#import gui_thread

sys.path.append(os.path.join(d,'pyunit-1.1.0'))
import unittest

from common import *
    
#---- testing ----#
test_modules = [common]

def test():
    #XXX make this work on everything automatically
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())
    return runner

def test_suite():
    #XXX make this work on everything automatically
    suites=[]
    for module in test_modules:
        suites.append(module.test_suite())
    total_suite = unittest.TestSuite(suites)
    return total_suite            

    
