"""Data input and output.

  Classes

     fopen -- a class for easily reading and writing binary data.

  Functions

     read_array -- reading ascii streams into Numeric arrays
     write_array -- write an array to an ascii stream
     loadmat -- read a MATLAB (version <= 4) style mat file
     savemat -- write a MATLAB (version <= 4) style mat file

     fread -- low-level reading
     fwrite -- low-level writing
     bswap -- in-place byte-swapping
     packbits -- Pack a binary array of 1's and 0's into an array of bytes
     unpackbits -- Unpack an array packed by packbits.
"""
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
