#! /usr/bin/env python
# Last Change: Fri Dec 15 10:00 PM 2006 J

from numpy.testing import *
from numpy.random import randn, seed
from numpy import correlate, array, concatenate, require, corrcoef
from numpy.fft import fft, ifft

from numpy.ctypeslib import ndpointer, load_library
from ctypes import c_uint

set_package_path()
from cdavid.autocorr import _raw_autocorr_1d, _raw_autocorr_1d_noncontiguous
from cdavid.autocorr import autocorr_oneside_nofft as autocorr
from cdavid.autocorr import autocorr_fft , nextpow2
from cdavid.autocorr import _autocorr_oneside_nofft_py as autocorr_py
restore_path()

import numpy

# number of decimals to check
nd  = 20
# minimum number of correct decimals required
md  = 12

a   = array([1, 2, 3.])
b   = a + 3

x   = concatenate((a, b)).reshape(2, 3)
        
# float and double C order
xc      = require(x, dtype = numpy.float64, requirements = 'C')
xcf     = require(x, dtype = numpy.float32, requirements = 'C')
xc1     = xc[0]
xcf1    = xcf[0]

# float and double F order
xf  = require(x, dtype = numpy.float64, requirements = 'FORTRAN')
xff = require(x, dtype = numpy.float32, requirements = 'FORTRAN')
xf1     = xf[0]
xff1    = xff[0]

# This class tests the C functions directly. This is more a debugging tool
# that a test case, as the tested functions are not part of the public API
class test_ctype_1d(NumpyTestCase):
    def check_contiguous_double(self):
        # double test
        xt      = xc1
        yt      = _raw_autocorr_1d(xt, xt.size - 1)

        yr  = correlate(xt, xt, mode = 'full')
        yr  = yr[xt.size-1:]

        assert_array_equal(yt, yr)

    def check_contiguous_float(self):
        # float test
        xt      = xcf1

        yt      = _raw_autocorr_1d(xt, xt.size - 1)

        yr  = correlate(xt, xt, mode = 'full')
        yr  = yr[xt.size-1:]

        assert_array_equal(yt, yr)

    def check_non_contiguous_double(self):
        # double test
        xt      = xf1
        yt      = _raw_autocorr_1d_noncontiguous(xt, xt.size - 1)

        yr  = correlate(xt, xt, mode = 'full')
        yr  = yr[xt.size-1:]

        assert_array_equal(yt, yr)

    def check_non_contiguous_float(self):
        # float test
        xt      = xff1
        yt      = _raw_autocorr_1d_noncontiguous(xt, xt.size - 1)

        yr  = correlate(xt, xt, mode = 'full')
        yr  = yr[xt.size-1:]

        assert_array_equal(yt, yr)

# Test autocorrelation for rank 1 arrays
class test_autocorr_1d(NumpyTestCase):
    def check_contiguous_double(self):
        # double test
        xt      = xc1
        yt      = autocorr(xt, xt.size - 1)

        yr  = correlate(xt, xt, mode = 'full')
        yr  = yr[xt.size-1:]

        assert_array_equal(yt, yr)

    def check_contiguous_float(self):
        # float test
        xt      = xcf1

        yt      = autocorr(xt, xt.size - 1)

        yr  = correlate(xt, xt, mode = 'full')
        yr  = yr[xt.size-1:]

        assert_array_equal(yt, yr)

    def check_non_contiguous_double(self):
        # double test
        xt      = xf1
        yt      = autocorr(xt, xt.size - 1)

        yr  = correlate(xt, xt, mode = 'full')
        yr  = yr[xt.size-1:]

        assert_array_equal(yt, yr)

    def check_non_contiguous_float(self):
        # float test
        xt      = xff1
        yt      = autocorr(xt, xt.size - 1)

        yr  = correlate(xt, xt, mode = 'full')
        yr  = yr[xt.size-1:]

        assert_array_equal(yt, yr)

# This class is a pure python implementation of autocorrelation
# with rank 2 arrays. This will be used in the above test cases;
# this function implements the expected behaviour of the public 
# autocorr function.
class test_autocorr_py(NumpyTestCase):
    def check_full(self):
        xt      = xc
        axis    = -1
        lag     = xt.shape[axis] - 1
        yt      = autocorr_py(xt, lag, axis = axis)

        yr  = yt.copy()
        for i in range(xt.shape[(axis +1) % 2]):
            tmp     = correlate(xt[i], xt[i], 'full')
            center  = xt[i].size - 1
            assert_array_equal(tmp[center:center+1+lag], yt[i])

        xt      = xc
        axis    = 0
        lag     = xt.shape[axis] - 1
        yt      = autocorr_py(xt, lag, axis = axis)

        yr  = yt.copy()
        for i in range(xt.shape[(axis +1) % 2]):
            tmp = correlate(xt[:, i], xt[:, i], 'full')
            center  = xt[:,i].size - 1
            assert_array_equal(tmp[center:center+1+lag], yt[:, i])

    def check_partial(self):
        xt      = xc
        axis    = -1
        lag     = 1
        yt      = autocorr_py(xt, lag, axis = axis)

        yr  = yt.copy()
        for i in range(xt.shape[(axis +1) % 2]):
            tmp = correlate(xt[i], xt[i], 'full')
            center  = xt[i].size - 1
            assert_array_equal(tmp[center:center+1+lag], yt[i])

        xt      = xc
        axis    = 0
        lag     = 1
        yt      = autocorr_py(xt, lag, axis = axis)

        yr  = yt.copy()
        for i in range(xt.shape[(axis +1) % 2]):
            tmp = correlate(xt[:, i], xt[:, i], 'full')
            center  = xt[:,i].size - 1
            assert_array_equal(tmp[center:center+1+lag], yt[:, i])

# Test autocorrelation for rank 2 arrays
class test_autocorr_2d(NumpyTestCase):
    def check_double_full(self):
        # C, axis 1 test
        xt      = xc
        axis    = -1
        lag     = xt.shape[axis] - 1
        yt      = autocorr(xt, lag, axis = axis)

        yr      = autocorr_py(xt, lag, axis = axis)
        assert_array_equal(yt, yr)

        # C, axis 0 test
        xt      = xc
        axis    = 0
        lag     = xt.shape[axis] - 1
        yt      = autocorr(xt, lag, axis = axis)

        yr      = autocorr_py(xt, lag, axis = axis)
        assert_array_equal(yt, yr)

        # F, axis 0 test
        xt      = xf
        axis    = 0
        lag     = xt.shape[axis] - 1
        yt      = autocorr(xt, lag, axis = axis)

        yr      = autocorr_py(xt, lag, axis = axis)
        assert_array_equal(yt, yr)

        # F, axis 1 test
        xt      = xf
        axis    = -1
        lag     = xt.shape[axis] - 1
        yt      = autocorr(xt, lag, axis = axis)

        yr      = autocorr_py(xt, lag, axis = axis)
        assert_array_equal(yt, yr)

    def check_float(self):
        # C, axis 1 test
        xt      = xcf
        axis    = -1
        lag     = xt.shape[axis] - 1
        yt      = autocorr(xt, lag, axis = axis)

        yr      = autocorr_py(xt, lag, axis = axis)
        assert_array_equal(yt, yr)

        # C, axis 0 test
        xt      = xcf
        axis    = 0
        lag     = xt.shape[axis] - 1
        yt      = autocorr(xt, lag, axis = axis)

        yr      = autocorr_py(xt, lag, axis = axis)
        assert_array_equal(yt, yr)

        # F, axis 0 test
        xt      = xff
        axis    = 0
        lag     = xt.shape[axis] - 1
        yt      = autocorr(xt, lag, axis = axis)

        yr      = autocorr_py(xt, lag, axis = axis)
        assert_array_equal(yt, yr)

        # F, axis 1 test
        xt      = xff
        axis    = -1
        lag     = xt.shape[axis] - 1
        yt      = autocorr(xt, lag, axis = axis)

        yr      = autocorr_py(xt, lag, axis = axis)
        assert_array_equal(yt, yr)

    def check_double_partial(self):
        # C, axis 1 test
        xt      = xc
        axis    = -1
        lag     = 1
        yt      = autocorr(xt, lag, axis = axis)

        yr      = autocorr_py(xt, lag, axis = axis)
        assert_array_equal(yt, yr)

        # C, axis 0 test
        xt      = xc
        axis    = 0
        lag     = 0
        yt      = autocorr(xt, lag, axis = axis)

        yr      = autocorr_py(xt, lag, axis = axis)
        assert_array_equal(yt, yr)

        # F, axis 0 test
        xt      = xf
        axis    = 1
        lag     = xt.shape[axis] - 1
        yt      = autocorr(xt, lag, axis = axis)

        yr      = autocorr_py(xt, lag, axis = axis)
        assert_array_equal(yt, yr)

        # F, axis 1 test
        xt      = xf
        axis    = -1
        lag     = 1
        yt      = autocorr(xt, lag, axis = axis)

        yr      = autocorr_py(xt, lag, axis = axis)
        assert_array_equal(yt, yr)

class test_autocorr_fft(NumpyTestCase):
    n   = 5
    d   = 3
    def check_nextpow2(self):
        assert(nextpow2(255) == 8)
        assert(nextpow2(256) == 8)
        assert(nextpow2(257) == 9)

    def check_r1r(self):
        """real case, rank 1"""
        a   = randn(self.n)

        aref    = correlate(a, a, mode = 'full')
        atest   = autocorr_fft(a)
        assert_array_almost_equal(atest, aref, decimal = md)
        assert atest.dtype == a.dtype

    def check_r1c(self):
        """complex case, rank 1"""
        a   = randn(self.n) + 1.0j * randn(self.n)

        atest   = autocorr_fft(a)
        aref    = numpy.sum(a * numpy.conj(a))
        assert_array_almost_equal(atest[self.n - 1], aref, decimal = md)
        assert atest.dtype == a.dtype

    def check_r2c(self):
        """complex case, rank 2"""
        pass

    def check_r2r(self):
        """real case, rank 2"""

        # axis 0
        a       = randn(self.n, self.d)
        axis    = 0

        c       = [correlate(a[:, i], a[:, i], mode = 'full') for i in range(self.d)]
        aref    = array(c).T

        atest   = autocorr_fft(a, axis = axis)
        assert_array_almost_equal(atest, aref, decimal = md)

        # axis 1
        a       = randn(self.n, self.d)
        axis    = 1

        c       = [correlate(a[i], a[i], mode = 'full') for i in range(self.n)]
        aref    = array(c)

        atest   = autocorr_fft(a, axis = axis)
        assert_array_almost_equal(atest, aref, decimal = md)

if __name__ == "__main__":
    NumpyTest().run()

#class test_autocorr_2d(NumpyTestCase):
#    def check_double(self):
#        # C, axis 1 test
#        xt      = xc
#        axis    = -1
#        lag     = xt.shape[axis] - 1
#        yt      = autocorr(xt, lag, axis = axis)
#
#        yr  = yt.copy()
#        for i in range(xt.shape[(axis +1) % 2]):
#            tmp = correlate(xt[i], xt[i], 'full')
#            assert_array_equal(tmp[lag:], yt[i])
#
#        # C, axis 0 test
#        xt      = xc
#        axis    = 0
#        lag     = xt.shape[axis] - 1
#        yt      = autocorr(xt, lag, axis = axis)
#
#        yr  = yt.copy()
#        for i in range(xt.shape[(axis +1) % 2]):
#            tmp = correlate(xt[:, i], xt[:, i], 'full')
#            assert_array_equal(tmp[lag:], yt[:, i])
#
#        # F, axis 0 test
#        xt      = xf
#        axis    = 0
#        lag     = xt.shape[axis] - 1
#        yt      = autocorr(xt, lag, axis = axis)
#
#        yr  = yt.copy()
#        for i in range(xt.shape[(axis +1) % 2]):
#            tmp = correlate(xt[:, i], xt[:, i], 'full')
#            assert_array_equal(tmp[lag:], yt[:, i])
#
#        # F, axis 1 test
#        xt      = xf
#        axis    = -1
#        lag     = xt.shape[axis] - 1
#        yt      = autocorr(xt, lag, axis = axis)
#
#        yr  = yt.copy()
#        for i in range(xt.shape[(axis +1) % 2]):
#            tmp = correlate(xt[i], xt[i], 'full')
#            assert_array_equal(tmp[lag:], yt[i])
#
#    def check_float(self):
#        # C, axis 1 test
#        xt      = xcf
#        axis    = -1
#        lag     = xt.shape[axis] - 1
#        yt      = autocorr(xt, lag, axis = axis)
#
#        yr  = yt.copy()
#        for i in range(xt.shape[(axis +1) % 2]):
#            tmp = correlate(xt[i], xt[i], 'full')
#            assert_array_equal(tmp[lag:], yt[i])
#
#        # C, axis 0 test
#        xt      = xcf
#        axis    = 0
#        lag     = xt.shape[axis] - 1
#        yt      = autocorr(xt, lag, axis = axis)
#
#        yr  = yt.copy()
#        for i in range(xt.shape[(axis +1) % 2]):
#            tmp = correlate(xt[:, i], xt[:, i], 'full')
#            assert_array_equal(tmp[lag:], yt[:, i])
#
#        # F, axis 0 test
#        xt      = xff
#        axis    = 0
#        lag     = xt.shape[axis] - 1
#        yt      = autocorr(xt, lag, axis = axis)
#
#        yr  = yt.copy()
#        for i in range(xt.shape[(axis +1) % 2]):
#            tmp = correlate(xt[:, i], xt[:, i], 'full')
#            assert_array_equal(tmp[lag:], yt[:, i])
#
#        # F, axis 1 test
#        xt      = xff
#        axis    = -1
#        lag     = xt.shape[axis] - 1
#        yt      = autocorr(xt, lag, axis = axis)
#
#        yr  = yt.copy()
#        for i in range(xt.shape[(axis +1) % 2]):
#            tmp = correlate(xt[i], xt[i], 'full')
#            assert_array_equal(tmp[lag:], yt[i])
#
