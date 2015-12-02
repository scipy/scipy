''' Tests for netcdf '''
from __future__ import division, print_function, absolute_import

import os
from os.path import join as pjoin, dirname
import shutil
import tempfile
import warnings
from io import BytesIO
from glob import glob
from contextlib import contextmanager

import numpy as np
from numpy.testing import (assert_, assert_allclose, assert_raises,
    assert_equal, run_module_suite)

from scipy.io.netcdf import netcdf_file

from scipy._lib._tmpdirs import in_tempdir

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


def check_simple(ncfileobj):
    '''Example fileobj tests '''
    assert_equal(ncfileobj.history, b'Created for a test')
    time = ncfileobj.variables['time']
    assert_equal(time.units, b'days since 2008-01-01')
    assert_equal(time.shape, (N_EG_ELS,))
    assert_equal(time[-1], N_EG_ELS-1)


def test_read_write_files():
    # test round trip for example file
    cwd = os.getcwd()
    try:
        tmpdir = tempfile.mkdtemp()
        os.chdir(tmpdir)
        with make_simple('simple.nc', 'w') as f:
            pass
        # read the file we just created in 'a' mode
        with netcdf_file('simple.nc', 'a') as f:
            check_simple(f)
            # add something
            f._attributes['appendRan'] = 1

        # To read the NetCDF file we just created::
        with netcdf_file('simple.nc') as f:
            # Using mmap is the default
            assert_(f.use_mmap)
            check_simple(f)
            assert_equal(f._attributes['appendRan'], 1)

        # Read it in append (and check mmap is off)
        with netcdf_file('simple.nc', 'a') as f:
            assert_(not f.use_mmap)
            check_simple(f)
            assert_equal(f._attributes['appendRan'], 1)

        # Now without mmap
        with netcdf_file('simple.nc', mmap=False) as f:
            # Using mmap is the default
            assert_(not f.use_mmap)
            check_simple(f)

        # To read the NetCDF file we just created, as file object, no
        # mmap.  When n * n_bytes(var_type) is not divisible by 4, this
        # raised an error in pupynere 1.0.12 and scipy rev 5893, because
        # calculated vsize was rounding up in units of 4 - see
        # http://www.unidata.ucar.edu/software/netcdf/docs/netcdf.html
        with open('simple.nc', 'rb') as fobj:
            with netcdf_file(fobj) as f:
                # by default, don't use mmap for file-like
                assert_(not f.use_mmap)
                check_simple(f)

        # Read file from fileobj, with mmap
        with open('simple.nc', 'rb') as fobj:
            with netcdf_file(fobj, mmap=True) as f:
                assert_(f.use_mmap)
                check_simple(f)

        # Again read it in append mode (adding another att)
        with open('simple.nc', 'r+b') as fobj:
            with netcdf_file(fobj, 'a') as f:
                assert_(not f.use_mmap)
                check_simple(f)
                f.createDimension('app_dim', 1)
                var = f.createVariable('app_var', 'i', ('app_dim',))
                var[:] = 42

        # And... check that app_var made it in...
        with netcdf_file('simple.nc') as f:
            check_simple(f)
            assert_equal(f.variables['app_var'][:], 42)

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
        check_simple(f2)

    # Test that error is raised if attempting mmap for sio
    eg_sio3 = BytesIO(str_val)
    assert_raises(ValueError, netcdf_file, eg_sio3, 'r', True)
    # Test 64-bit offset write / read
    eg_sio_64 = BytesIO()
    with make_simple(eg_sio_64, 'w', version=2) as f_64:
        str_val = eg_sio_64.getvalue()

    eg_sio_64 = BytesIO(str_val)
    with netcdf_file(eg_sio_64) as f_64:
        check_simple(f_64)
        assert_equal(f_64.version_byte, 2)
    # also when version 2 explicitly specified
    eg_sio_64 = BytesIO(str_val)
    with netcdf_file(eg_sio_64, version=2) as f_64:
        check_simple(f_64)
        assert_equal(f_64.version_byte, 2)


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

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

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
            assert_raises(ValueError, f.createVariable, 'time', dt, ('time',))


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


def test_mmaps_segfault():
    filename = pjoin(TEST_DATA_PATH, 'example_1.nc')

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with netcdf_file(filename, mmap=True) as f:
            x = f.variables['lat'][:]
            # should not raise warnings
            del x

    def doit():
        with netcdf_file(filename, mmap=True) as f:
            return f.variables['lat'][:]

    # should not crash
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x = doit()
        x.sum()


def test_zero_dimensional_var():
    io = BytesIO()
    with make_simple(io, 'w') as f:
        v = f.createVariable('zerodim', 'i2', [])
        # This is checking that .isrec returns a boolean - don't simplify it
        # to 'assert not ...'
        assert v.isrec is False, v.isrec
        f.flush()


def test_byte_gatts():
    # Check that global "string" atts work like they did before py3k
    # unicode and general bytes confusion
    with in_tempdir():
        filename = 'g_byte_atts.nc'
        f = netcdf_file(filename, 'w')
        f._attributes['holy'] = b'grail'
        f._attributes['witch'] = 'floats'
        f.close()
        f = netcdf_file(filename, 'r')
        assert_equal(f._attributes['holy'], b'grail')
        assert_equal(f._attributes['witch'], b'floats')
        f.close()


def test_open_append():
    # open 'w' put one attr
    with in_tempdir():
        filename = 'append_dat.nc'
        f = netcdf_file(filename, 'w')
        f._attributes['Kilroy'] = 'was here'
        f.close()

        # open again in 'a', read the att and and a new one
        f = netcdf_file(filename, 'a')
        assert_equal(f._attributes['Kilroy'], b'was here')
        f._attributes['naughty'] = b'Zoot'
        f.close()

        # open yet again in 'r' and check both atts
        f = netcdf_file(filename, 'r')
        assert_equal(f._attributes['Kilroy'], b'was here')
        assert_equal(f._attributes['naughty'], b'Zoot')
        f.close()


def test_maskandscale():
    t = np.linspace(20, 30, 15)
    t[3] = 100
    tm = np.ma.masked_greater(t, 99)
    fname = pjoin(TEST_DATA_PATH, 'example_2.nc')
    with netcdf_file(fname, maskandscale=True) as f:
        Temp = f.variables['Temperature']
        assert_equal(Temp.missing_value, 9999)
        assert_equal(Temp.add_offset, 20)
        assert_equal(Temp.scale_factor, np.float32(0.01))
        found = Temp[:].compressed()
        del Temp  # Remove ref to mmap, so file can be closed.
        expected = np.round(tm.compressed(), 2)
        assert_allclose(found, expected)

    with in_tempdir():
        newfname = 'ms.nc'
        f = netcdf_file(newfname, 'w', maskandscale=True)
        f.createDimension('Temperature', len(tm))
        temp = f.createVariable('Temperature', 'i', ('Temperature',))
        temp.missing_value = 9999
        temp.scale_factor = 0.01
        temp.add_offset = 20
        temp[:] = tm
        f.close()

        with netcdf_file(newfname, maskandscale=True) as f:
            Temp = f.variables['Temperature']
            assert_equal(Temp.missing_value, 9999)
            assert_equal(Temp.add_offset, 20)
            assert_equal(Temp.scale_factor, np.float32(0.01))
            expected = np.round(tm.compressed(), 2)
            found = Temp[:].compressed()
            del Temp
            assert_allclose(found, expected)


if __name__ == "__main__":
    run_module_suite()
