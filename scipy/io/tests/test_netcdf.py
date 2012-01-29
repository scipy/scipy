''' Tests for netcdf '''

import os
from os.path import join as pjoin, dirname
import shutil
import tempfile
import time
import sys
if sys.version_info[0] >= 3:
    from io import BytesIO
else:
    from StringIO import StringIO as BytesIO
from glob import glob

import numpy as np
from numpy.compat import asbytes

from scipy.io.netcdf import netcdf_file

from nose.tools import assert_true, assert_false, assert_equal, assert_raises

TEST_DATA_PATH = pjoin(dirname(__file__), 'data')

N_EG_ELS = 11 # number of elements for example variable
VARTYPE_EG = 'b' # var type for example variable


def make_simple(*args, **kwargs):
    f = netcdf_file(*args, **kwargs)
    f.history = 'Created for a test'
    f.createDimension('time', N_EG_ELS)
    time = f.createVariable('time', VARTYPE_EG, ('time',))
    time[:] = np.arange(N_EG_ELS)
    time.units = 'days since 2008-01-01'
    f.flush()
    return f


def gen_for_simple(ncfileobj):
    ''' Generator for example fileobj tests '''
    yield assert_equal, ncfileobj.history, asbytes('Created for a test')
    time = ncfileobj.variables['time']
    yield assert_equal, time.units, asbytes('days since 2008-01-01')
    yield assert_equal, time.shape, (N_EG_ELS,)
    yield assert_equal, time[-1], N_EG_ELS-1


def test_read_write_files():
    # test round trip for example file
    cwd = os.getcwd()
    try:
        tmpdir = tempfile.mkdtemp()
        os.chdir(tmpdir)
        f = make_simple('simple.nc', 'w')
        f.close()
        # To read the NetCDF file we just created::
        f = netcdf_file('simple.nc')
        # Using mmap is the default
        yield assert_true, f.use_mmap
        for testargs in gen_for_simple(f):
            yield testargs
        f.close()
        # Now without mmap
        f = netcdf_file('simple.nc', mmap=False)
        # Using mmap is the default
        yield assert_false, f.use_mmap
        for testargs in gen_for_simple(f):
            yield testargs
        f.close()
        # To read the NetCDF file we just created, as file object, no
        # mmap.  When n * n_bytes(var_type) is not divisible by 4, this
        # raised an error in pupynere 1.0.12 and scipy rev 5893, because
        # calculated vsize was rounding up in units of 4 - see
        # http://www.unidata.ucar.edu/software/netcdf/docs/netcdf.html
        fobj = open('simple.nc', 'rb')
        f = netcdf_file(fobj)
        # by default, don't use mmap for file-like
        yield assert_false, f.use_mmap
        for testargs in gen_for_simple(f):
            yield testargs
        f.close()
    except:
        os.chdir(cwd)
        shutil.rmtree(tmpdir)
        raise
    os.chdir(cwd)
    shutil.rmtree(tmpdir)


def test_read_write_sio():
    eg_sio1 = BytesIO()
    f1 = make_simple(eg_sio1, 'w')
    str_val = eg_sio1.getvalue()
    f1.close()
    eg_sio2 = BytesIO(str_val)
    f2 = netcdf_file(eg_sio2)
    for testargs in gen_for_simple(f2):
        yield testargs
    f2.close()
    # Test that error is raised if attempting mmap for sio
    eg_sio3 = BytesIO(str_val)
    yield assert_raises, ValueError, netcdf_file, eg_sio3, 'r', True
    # Test 64-bit offset write / read
    eg_sio_64 = BytesIO()
    f_64 = make_simple(eg_sio_64, 'w', version=2)
    str_val = eg_sio_64.getvalue()
    f_64.close()
    eg_sio_64 = BytesIO(str_val)
    f_64 = netcdf_file(eg_sio_64)
    for testargs in gen_for_simple(f_64):
        yield testargs
    yield assert_equal, f_64.version_byte, 2
    # also when version 2 explicitly specified
    eg_sio_64 = BytesIO(str_val)
    f_64 = netcdf_file(eg_sio_64, version=2)
    for testargs in gen_for_simple(f_64):
        yield testargs
    yield assert_equal, f_64.version_byte, 2


def test_read_example_data():
    # read any example data files
    for fname in glob(pjoin(TEST_DATA_PATH, '*.nc')):
        f = netcdf_file(fname, 'r')
        f = netcdf_file(fname, 'r', mmap=False)

def test_itemset_no_segfault_on_readonly():
    # Regression test for ticket #1202.
    # Open the test file in read-only mode.
    filename = pjoin(TEST_DATA_PATH, 'example_1.nc')
    f = netcdf_file(filename, 'r')
    time_var = f.variables['time']
    # time_var.assignValue(42) should raise a RuntimeError--not seg. fault!
    assert_raises(RuntimeError, time_var.assignValue, 42)

def test_write_invalid_dtype():
    dtypes = ['int64', 'uint64']
    if np.dtype('int').itemsize == 8:   # 64-bit machines
        dtypes.append('int')
    if np.dtype('uint').itemsize == 8:   # 64-bit machines
        dtypes.append('uint')

    f = netcdf_file(BytesIO(), 'w')
    f.createDimension('time', N_EG_ELS)
    for dt in dtypes:
        yield assert_raises, ValueError, \
            f.createVariable, 'time', dt, ('time',)
    f.close()

