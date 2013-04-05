''' Tests for fortran sequential files '''

import os
import tempfile
import shutil
from os import path
import warnings

DATA_PATH = path.join(path.dirname(__file__), 'data')

from numpy.testing import assert_equal
import numpy as np
from glob import iglob
import re

from scipy.io.fortran import fortran_file

def test_fortranfiles_read():
    for filename in iglob(path.join(DATA_PATH, "fortran-*-*x*x*.dat")):
        m = re.search('fortran-([^-]+)-(\d+)x(\d+)x(\d+).dat', filename)
        if not m:
            raise RuntimeError("Couldn't match %s filename to regex" % filename)
        dims = (int(m.group(2)), int(m.group(3)), int(m.group(4)))

        f = fortran_file(filename, 'r', '<i4')
        data = f.readRecord(dtype=m.group(1)).reshape(dims)
        f.close()

        counter = 0
        for k in range(dims[2]):
            for j in range(dims[1]):
                for i in range(dims[0]):
                    assert_equal(counter, data[i,j,k])
                    counter += 1

def test_fortranfiles_write():
    for filename in iglob(path.join(DATA_PATH, "fortran-*-*x*x*.dat")):
        m = re.search('fortran-([^-]+)-(\d+)x(\d+)x(\d+).dat', filename)
        if not m:
            raise RuntimeError("Couldn't match %s filename to regex" % filename)
        dims = (int(m.group(2)), int(m.group(3)), int(m.group(4)))

        counter = 0
        data = np.zeros(dims, dtype=m.group(1))
        for k in range(dims[2]):
            for j in range(dims[1]):
                for i in range(dims[0]):
                    data[i,j,k] = counter
                    counter += 1
        try:
            tmpdir = tempfile.mkdtemp()
            testFile = path.join(tmpdir,path.basename(filename))
            f = fortran_file(testFile, 'w','<i4')
            f.writeRecord(data)
            f.close()
            originalfile = open(filename, 'rb')
            newfile = open(testFile, 'rb')
            assert_equal(originalfile.read(), newfile.read())
        except:
            shutil.rmtree(tmpdir)
            f.close()
            originalfile.close()
            newfile.close()
            raise
        shutil.rmtree(tmpdir)
        f.close()
        originalfile.close()
        newfile.close()
