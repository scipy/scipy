''' Tests for fortran sequential files '''

import tempfile
import shutil
from os import path
from glob import iglob
import re

from numpy.testing import assert_equal, assert_allclose, run_module_suite
import numpy as np

from scipy.io import FortranFile


DATA_PATH = path.join(path.dirname(__file__), 'data')


def test_fortranfiles_read():
    for filename in iglob(path.join(DATA_PATH, "fortran-*-*x*x*.dat")):
        m = re.search('fortran-([^-]+)-(\d+)x(\d+)x(\d+).dat', filename, re.I)
        if not m:
            raise RuntimeError("Couldn't match %s filename to regex" % filename)

        dims = (int(m.group(2)), int(m.group(3)), int(m.group(4)))

        f = FortranFile(filename, 'r', '<u4')
        data = f.read_record(dtype=m.group(1).replace('s', '<')).reshape(dims)
        f.close()

        counter = 0
        for k in range(dims[2]):
            for j in range(dims[1]):
                for i in range(dims[0]):
                    assert_equal(counter, data[i,j,k])
                    counter += 1


def test_fortranfiles_mixed_record():
    filename = path.join(DATA_PATH, "fortran-mixed.dat")
    with FortranFile(filename, 'r', '<u4') as f:
        record = f.read_record('<i4,<f4,<i8,(2)<f8')

    assert_equal(record['f0'][0], 1)
    assert_allclose(record['f1'][0], 2.3)
    assert_equal(record['f2'][0], 4)
    assert_allclose(record['f3'][0], [5.6, 7.8])


def test_fortranfiles_write():
    for filename in iglob(path.join(DATA_PATH, "fortran-*-*x*x*.dat")):
        m = re.search('fortran-([^-]+)-(\d+)x(\d+)x(\d+).dat', filename, re.I)
        if not m:
            raise RuntimeError("Couldn't match %s filename to regex" % filename)
        dims = (int(m.group(2)), int(m.group(3)), int(m.group(4)))

        counter = 0
        data = np.zeros(dims, dtype=m.group(1).replace('s', '<'))
        for k in range(dims[2]):
            for j in range(dims[1]):
                for i in range(dims[0]):
                    data[i,j,k] = counter
                    counter += 1
        tmpdir = tempfile.mkdtemp()
        try:
            testFile = path.join(tmpdir,path.basename(filename))
            f = FortranFile(testFile, 'w','<u4')
            f.write_record(data)
            f.close()
            originalfile = open(filename, 'rb')
            newfile = open(testFile, 'rb')
            assert_equal(originalfile.read(), newfile.read(),
                         err_msg=filename)
            originalfile.close()
            newfile.close()
        finally:
            shutil.rmtree(tmpdir)


if __name__ == "__main__":
    run_module_suite()
