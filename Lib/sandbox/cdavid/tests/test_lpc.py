#! /usr/bin/env python
# Last Change: Tue Nov 28 05:00 PM 2006 J

from numpy.testing import *
from numpy.random import randn, seed
from numpy import correlate, array, concatenate, require

from numpy.ctypeslib import ndpointer, load_library
from ctypes import c_uint

set_package_path()
from cdavid.lpc import _lpc2_py as lpc_py
from cdavid.lpc import lpc_ref, lpc2
from cdavid.autocorr import autocorr_oneside_nofft
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

# This class uses lpc in 1 dimension and loop on the axis. Is tested against
# a direct matrix inversion of the autocorrelation matrix (using matrix inverse 
# instead of levinson durbin)
class test_lpc_py(NumpyTestCase):
    def check_float(self):
        # Axis -1
        xt      = xcf
        axis    = -1
        order   = 1

        a, k, e = lpc_py(xt, order, axis)
        assert a.dtype == k.dtype == e.dtype == numpy.float32

        tmp     = numpy.zeros((xt.shape[0], order+1), xt.dtype)
        for i in range(xt.shape[0]):
            tmp[i]  = lpc_ref(xt[i], order)

        assert_array_almost_equal(tmp, a)

        # Axis 0
        xt      = xcf
        axis    = 0
        order   = 1

        a, k, e = lpc_py(xt, order, axis)
        assert a.dtype == k.dtype == e.dtype == numpy.float32

        tmp     = numpy.zeros((order + 1, xt.shape[1]), xt.dtype)
        for i in range(xt.shape[1]):
            tmp[:, i]   = lpc_ref(xt[:, i], order)

        assert_array_almost_equal(tmp, a)

    def check_double(self):
        # Axis -1 
        xt      = xc
        axis    = -1
        order   = 1

        a, e, k = lpc_py(xt, order, axis)
        assert a.dtype == k.dtype == e.dtype == numpy.float64

        tmp     = numpy.zeros((xt.shape[0], order+1), xt.dtype)
        for i in range(xt.shape[0]):
            tmp[i]  = lpc_ref(xt[i], order)

        assert_array_almost_equal(tmp, a)

        # Axis 0
        xt      = xc
        axis    = 0
        order   = 1

        a, e, k = lpc_py(xt, order, axis)
        assert a.dtype == k.dtype == e.dtype == numpy.float64

        tmp     = numpy.zeros((order + 1, xt.shape[1]), xt.dtype)
        for i in range(xt.shape[1]):
            tmp[:, i]   = lpc_ref(xt[:, i], order)

        assert_array_almost_equal(tmp, a)

class test_lpc(NumpyTestCase):
    def check_float(self):
        # Axis -1 
        xt      = xcf
        axis    = -1
        order   = 1

        a, e, k     = lpc2(xt, order, axis)
        at, et, kt  = lpc_py(xt, order, axis)

        assert a.dtype == e.dtype == k.dtype == numpy.float32

        assert_array_almost_equal(a, at)
        assert_array_almost_equal(e, et)
        assert_array_almost_equal(k, kt)

        # Axis 0 
        xt      = xcf
        axis    = 0
        order   = 1

        a, e, k     = lpc2(xt, order, axis)
        at, et, kt  = lpc_py(xt, order, axis)

        assert a.dtype == e.dtype == k.dtype == numpy.float32

        assert_array_almost_equal(a, at)
        assert_array_almost_equal(e, et)
        assert_array_almost_equal(k, kt)

    def check_float_rank1(self):
        # test rank 1
        xt      = xcf[0]
        axis    = 0
        order   = 1

        a, e, k     = lpc2(xt, order, axis)
        at, et, kt  = lpc_py(xt, order, axis)

        assert a.dtype == e.dtype == k.dtype == numpy.float32

        assert_array_almost_equal(a, at)
        assert_array_almost_equal(e, et)
        assert_array_almost_equal(k, kt)

    def check_double(self):
        # Axis -1 
        xt      = xc
        axis    = -1
        order   = 1

        a, e, k     = lpc2(xt, order, axis)
        at, et, kt  = lpc_py(xt, order, axis)

        assert_array_almost_equal(a, at)
        assert_array_almost_equal(e, et)
        assert_array_almost_equal(k, kt)

        # Axis 0 
        xt      = xc
        axis    = 0
        order   = 1

        a, e, k     = lpc2(xt, order, axis)
        at, et, kt  = lpc_py(xt, order, axis)

        assert_array_almost_equal(a, at)
        assert_array_almost_equal(e, et)
        assert_array_almost_equal(k, kt)

    def check_double_rank1(self):
        # test rank 1
        xt      = xc[0]
        axis    = 0
        order   = 1

        a, e, k     = lpc2(xt, order, axis)
        at, et, kt  = lpc_py(xt, order, axis)

        assert_array_almost_equal(a, at)
        assert_array_almost_equal(e, et)
        assert_array_almost_equal(k, kt)

if __name__ == "__main__":
    ScipyTest().run()
