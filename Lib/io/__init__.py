#
# io - Data input and output
#

import os as _os
execfile(_os.path.join(__path__[0],'pre___init__.py'),globals(),locals())

from numpyio import packbits, unpackbits, bswap, fread, fwrite, \
     convert_objectarray
from mio import *
from array_import import *
from data_store import *
from pickler import *


#---- testing ----#
def test(level=10):
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite(level=level))
    return runner

def test_suite(level=1):
    import scipy_test.testing
    import scipy.io
    this_mod = scipy.io
    return scipy_test.testing.harvest_test_suites(this_mod,level=level)
