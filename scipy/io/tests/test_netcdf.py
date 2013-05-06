''' Tests for netcdf '''
from __future__ import division, print_function, absolute_import

import os
from os.path import join as pjoin, dirname
import shutil
import tempfile
import time
import sys
from io import BytesIO
from glob import glob
from contextlib import contextmanager

import numpy as np
from numpy.testing import dec, assert_, assert_allclose

from scipy.io.netcdf import netcdf_file

from nose.tools import assert_true, assert_false, assert_equal, assert_raises

TEST_DATA_PATH = pjoin(dirname(__file__), 'data')

N_EG_ELS = 11  # number of elements for example variable
VARTYPE_EG = 'b'  # var type for example variable


@contextmanager
def make_simple(*args, **kwargs):
    f = netcdf_file(*args, **kwargs)
    f.history = 'Created for a test'
    f.createDimension('time', N_EG_ELS)
    time = f.createVariable('time', VARTYPE_EG, ('time',))
    time[:] = np.arange(N_EG_ELS)
    time.units = 'days since 2008-01-01'
    f.flush()
    yield f
    f.close()


def gen_for_simple(ncfileobj):
    ''' Generator for example fileobj tests '''
    yield assert_equal, ncfileobj.history, b'Created for a test'
    time = ncfileobj.variables['time']
    yield assert_equal, time.units, b'days since 2008-01-01'
    yield assert_equal, time.shape, (N_EG_ELS,)
    yield assert_equal, time[-1], N_EG_ELS-1


def test_read_write_files():
    # test round trip for example file
    cwd = os.getcwd()
    try:
        tmpdir = tempfile.mkdtemp()
        os.chdir(tmpdir)
        with make_simple('simple.nc', 'w') as f:
            pass
        # To read the NetCDF file we just created::
        with netcdf_file('simple.nc') as f:
            # Using mmap is the default
            yield assert_true, f.use_mmap
            for testargs in gen_for_simple(f):
                yield testargs

        # Now without mmap
        with netcdf_file('simple.nc', mmap=False) as f:
            # Using mmap is the default
            yield assert_false, f.use_mmap
            for testargs in gen_for_simple(f):
                yield testargs

        # To read the NetCDF file we just created, as file object, no
        # mmap.  When n * n_bytes(var_type) is not divisible by 4, this
        # raised an error in pupynere 1.0.12 and scipy rev 5893, because
        # calculated vsize was rounding up in units of 4 - see
        # http://www.unidata.ucar.edu/software/netcdf/docs/netcdf.html
        fobj = open('simple.nc', 'rb')
        with netcdf_file(fobj) as f:
            # by default, don't use mmap for file-like
            yield assert_false, f.use_mmap
            for testargs in gen_for_simple(f):
                yield testargs
    except:
        os.chdir(cwd)
        shutil.rmtree(tmpdir)
        raise
    os.chdir(cwd)
    shutil.rmtree(tmpdir)


def test_read_write_sio():
    eg_sio1 = BytesIO()
    with make_simple(eg_sio1, 'w') as f1:
        str_val = eg_sio1.getvalue()

    eg_sio2 = BytesIO(str_val)
    with netcdf_file(eg_sio2) as f2:
        for testargs in gen_for_simple(f2):
            yield testargs

    # Test that error is raised if attempting mmap for sio
    eg_sio3 = BytesIO(str_val)
    yield assert_raises, ValueError, netcdf_file, eg_sio3, 'r', True
    # Test 64-bit offset write / read
    eg_sio_64 = BytesIO()
    with make_simple(eg_sio_64, 'w', version=2) as f_64:
        str_val = eg_sio_64.getvalue()

    eg_sio_64 = BytesIO(str_val)
    with netcdf_file(eg_sio_64) as f_64:
        for testargs in gen_for_simple(f_64):
            yield testargs
        yield assert_equal, f_64.version_byte, 2
    # also when version 2 explicitly specified
    eg_sio_64 = BytesIO(str_val)
    with netcdf_file(eg_sio_64, version=2) as f_64:
        for testargs in gen_for_simple(f_64):
            yield testargs
        yield assert_equal, f_64.version_byte, 2


def test_read_example_data():
    # read any example data files
    for fname in glob(pjoin(TEST_DATA_PATH, '*.nc')):
        with netcdf_file(fname, 'r') as f:
            pass
        with netcdf_file(fname, 'r', mmap=False) as f:
            pass


def test_itemset_no_segfault_on_readonly():
    # Regression test for ticket #1202.
    # Open the test file in read-only mode.
    filename = pjoin(TEST_DATA_PATH, 'example_1.nc')
    with netcdf_file(filename, 'r') as f:
        time_var = f.variables['time']

    # time_var.assignValue(42) should raise a RuntimeError--not seg. fault!
    assert_raises(RuntimeError, time_var.assignValue, 42)


def test_write_invalid_dtype():
    dtypes = ['int64', 'uint64']
    if np.dtype('int').itemsize == 8:   # 64-bit machines
        dtypes.append('int')
    if np.dtype('uint').itemsize == 8:   # 64-bit machines
        dtypes.append('uint')

    with netcdf_file(BytesIO(), 'w') as f:
        f.createDimension('time', N_EG_ELS)
        for dt in dtypes:
            yield assert_raises, ValueError, \
                f.createVariable, 'time', dt, ('time',)


def test_flush_rewind():
    stream = BytesIO()
    with make_simple(stream, mode='w') as f:
        x = f.createDimension('x',4)
        v = f.createVariable('v', 'i2', ['x'])
        v[:] = 1
        f.flush()
        len_single = len(stream.getvalue())
        f.flush()
        len_double = len(stream.getvalue())

    assert_(len_single == len_double)


def test_dtype_specifiers():
    # Numpy 1.7.0-dev had a bug where 'i2' wouldn't work.
    # Specifying np.int16 or similar only works from the same commit as this
    # comment was made.
    with make_simple(BytesIO(), mode='w') as f:
        f.createDimension('x',4)
        f.createVariable('v1', 'i2', ['x'])
        f.createVariable('v2', np.int16, ['x'])
        f.createVariable('v3', np.dtype(np.int16), ['x'])


def test_ticket_1720():
    io = BytesIO()

    items = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

    with netcdf_file(io, 'w') as f:
        f.history = 'Created for a test'
        f.createDimension('float_var', 10)
        float_var = f.createVariable('float_var', 'f', ('float_var',))
        float_var[:] = items
        float_var.units = 'metres'
        f.flush()
        contents = io.getvalue()

    io = BytesIO(contents)
    with netcdf_file(io, 'r') as f:
        assert_equal(f.history, b'Created for a test')
        float_var = f.variables['float_var']
        assert_equal(float_var.units, b'metres')
        assert_equal(float_var.shape, (10,))
        assert_allclose(float_var[:], items)
