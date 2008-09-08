#!/usr/bin/env python

# This python script tests the numpyio module.
# also check out numpyio.fread.__doc__ and other method docstrings.

import os
from numpy.testing import *
import scipy.io as io
from scipy.io import numpyio
from scipy.io import array_import

import numpy.oldnumeric as N
import tempfile

class TestNumpyio(TestCase):
    def test_basic(self):
        # Generate some data
        a = 255*rand(20)
        # Open a file
        fname = tempfile.mktemp('.dat')
        fid = open(fname,"wb")
        # Write the data as shorts
        numpyio.fwrite(fid,20,a,N.Int16)
        fid.close()
        # Reopen the file and read in data
        fid = open(fname,"rb")
        if verbose >= 3:
            print "\nDon't worry about a warning regarding the number of bytes read."
        b = numpyio.fread(fid,1000000,N.Int16,N.Int)
        fid.close()
        assert(N.product(a.astype(N.Int16) == b,axis=0))
        os.remove(fname)

class TestReadArray(TestCase):
    def test_complex(self):
        a = rand(13,4) + 1j*rand(13,4)
        fname = tempfile.mktemp('.dat')
        io.write_array(fname,a)
        b = io.read_array(fname,atype=N.Complex)
        assert_array_almost_equal(a,b,decimal=4)
        os.remove(fname)

    def test_float(self):
        a = rand(3,4)*30
        fname = tempfile.mktemp('.dat')
        io.write_array(fname,a)
        b = io.read_array(fname)
        assert_array_almost_equal(a,b,decimal=4)
        os.remove(fname)

    def test_integer(self):
        from scipy import stats
        a = stats.randint.rvs(1,20,size=(3,4))
        fname = tempfile.mktemp('.dat')
        io.write_array(fname,a)
        b = io.read_array(fname,atype=a.dtype.char)
        assert_array_equal(a,b)
        os.remove(fname)

class TestRegression(TestCase):
    def test_get_open_file_works_with_filelike_objects(self):
        f = tempfile.TemporaryFile()
        f2 = array_import.get_open_file(f)
        assert f2 is f
        f.close()

if __name__ == "__main__":
    run_module_suite()
