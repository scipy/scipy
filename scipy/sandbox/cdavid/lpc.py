#! /usr/bin/env python

# Last Change: Tue Nov 28 03:00 PM 2006 J

# TODO: make lpc work for any rank
from warnings import warn

import numpy as _N
import scipy.signal as _sig

from scipy.linalg import toeplitz, inv

from autocorr import autocorr_oneside_nofft

# use ctype to have fast lpc computation
import ctypes
from ctypes import c_uint, c_int
from numpy.ctypeslib import ndpointer, load_library
ctypes_major    = int(ctypes.__version__.split('.')[0])
if ctypes_major < 1:
    msg =  "version of ctypes is %s, expected at least %s" \
            % (ctypes.__version__, '1.0.1')
    raise importerror(msg)

# load lpc
_lpc   = load_library('gabsig.so', __file__)

#===============================
# define the functions with args
#===============================
arg1    = ndpointer(dtype = _N.float64, flags = 'C')
arg2    = c_uint
arg3    = c_uint
arg4    = ndpointer(dtype = _N.float64, flags = 'C')
arg5    = ndpointer(dtype = _N.float64, flags = 'C')
arg6    = ndpointer(dtype = _N.float64, flags = 'C')

_lpc.dbl_lpc.argtypes   = [arg1, arg2, arg3, arg4, arg5, arg6]
_lpc.dbl_lpc.restype    = c_int

arg1    = ndpointer(dtype = _N.float32, flags = 'C')
arg2    = c_uint
arg3    = c_uint
arg4    = ndpointer(dtype = _N.float32, flags = 'C')
arg5    = ndpointer(dtype = _N.float32, flags = 'C')
arg6    = ndpointer(dtype = _N.float32, flags = 'C')

_lpc.flt_lpc.argtypes   = [arg1, arg2, arg3, arg4, arg5, arg6]
_lpc.flt_lpc.restype    = c_int

arg1    = ndpointer(dtype = _N.float64, flags = 'C')
arg2    = c_uint
arg3    = c_uint
arg4    = ndpointer(dtype = _N.float64, flags = 'C')
arg5    = ndpointer(dtype = _N.float64, flags = 'C')
arg6    = ndpointer(dtype = _N.float64, flags = 'C')

_lpc.dbl_levinson2d.argtypes   = [arg1, arg2, arg3, arg4, arg5, arg6]
_lpc.dbl_levinson2d.restype    = c_int

arg1    = ndpointer(dtype = _N.float32, flags = 'C')
arg2    = c_uint
arg3    = c_uint
arg4    = ndpointer(dtype = _N.float32, flags = 'C')
arg5    = ndpointer(dtype = _N.float32, flags = 'C')
arg6    = ndpointer(dtype = _N.float32, flags = 'C')

_lpc.flt_levinson2d.argtypes   = [arg1, arg2, arg3, arg4, arg5, arg6]
_lpc.flt_levinson2d.restype    = c_int

def lpc_ref(signal, order):
    """ Return the order + 1 LPC coefficients
    for the signal. This is just for reference, as it is using
    the direct inversion of the toeplitz matrix, which is really slow"""
    if signal.ndim > 1:
        print "Warning, not tested for rank > 1"
    p   = order + 1
    r   = autocorr_oneside_nofft(signal, order) / signal.shape[0]
    return _N.concatenate(([1.],  _N.dot(inv(toeplitz(r[:-1])), -r[1:])))

def lpc(signal, order):
    """ Return the order + 1 LPC coefficients
    for the signal using levinson durbin algorithm """
    if signal.ndim > 1:
        warn("Warning, not tested for rank > 1")
    if signal.size <= order:
        raise RuntimeError("size is smaller than order !")

    if signal.dtype == _N.float64:
        coeff   = _N.zeros((order+1), _N.float64)
        kcoeff  = _N.zeros((order), _N.float64)
        err     = _N.zeros((1), _N.float64)
        st  = _lpc.dbl_lpc(signal, signal.size, order, coeff, kcoeff, err)
    elif signal.dtype == _N.float32:
        coeff   = _N.zeros((order+1), _N.float32)
        kcoeff  = _N.zeros((order), _N.float32)
        err     = _N.zeros((1), _N.float32)
        st  = _lpc.flt_lpc(signal, signal.size, order, coeff, kcoeff, err)
    else:
        raise TypeError("Sorry, only float32 and float64 supported")
    if not (st == 0):
        raise RuntimeError("Error while using levinson algo, returns err is %d", st)
    return coeff, err, kcoeff

def lpcres_ref(signal, order, axis = 0):
    return _sig.lfilter(lpc_ref(signal, order), 1., signal, axis = 0)

def lpcres(signal, order, axis = 0):
    return _sig.lfilter(lpc(signal, order)[0], 1., signal, axis = 0)

def _lpc2_py(signal, order, axis = -1):
    """python implementation of lpc for rank 2., Do not use, for testing purpose only"""
    if signal.ndim > 2:
        raise NotImplemented("only for rank <=2")

    if signal.ndim < 2:
        return lpc(_N.require(signal, requirements = 'C'), order)

    # For each array of direction axis, compute levinson durbin
    if axis  % 2 == 0:
        # Prepare output arrays
        coeff   = _N.zeros((order+1, signal.shape[1]), signal.dtype)
        kcoeff  = _N.zeros((order, signal.shape[1]), signal.dtype)
        err     = _N.zeros(signal.shape[1], signal.dtype)
        for i in range(signal.shape[1]):
            coeff[:, i], err[i], kcoeff[:, i] = \
                    lpc(_N.require(signal[:, i], requirements = 'C'), order)
    elif axis % 2 == 1:
        # Prepare output arrays
        coeff   = _N.zeros((signal.shape[0], order+1), signal.dtype)
        kcoeff  = _N.zeros((signal.shape[0], order), signal.dtype)
        err     = _N.zeros(signal.shape[0], signal.dtype)
        for i in range(signal.shape[0]):
            coeff[i], err[i], kcoeff[i] = \
                lpc(_N.require(signal[i], requirements = 'C'), order)
    else:
        raise RuntimeError("this should not happen, please fill a bug")

    return coeff, err, kcoeff

def lpc2(signal, order, axis = -1):
    """ Returns ar coeff, err and k coeff"""
    sz      = signal.shape[axis]
    if order >= sz:
        raise RuntimeError("order should be strictly smaller than the length of signal")

    # rank 1
    if signal.ndim < 2:
        if signal.dtype == _N.float64:
            coeff   = _N.zeros((order+1), _N.float64)
            kcoeff  = _N.zeros((order), _N.float64)
            err     = _N.zeros((1), _N.float64)
            st  = _lpc.dbl_lpc(signal, signal.size, order, coeff, kcoeff, err)
        elif signal.dtype == _N.float32:
            coeff   = _N.zeros((order+1), _N.float32)
            kcoeff  = _N.zeros((order), _N.float32)
            err     = _N.zeros((1), _N.float32)
            st  = _lpc.flt_lpc(signal, signal.size, order, coeff, kcoeff, err)
        else:
            raise TypeError("Sorry, only float32 and float64 supported")

        if not (st == 0):
            raise RuntimeError("Error while using levinson algo, returns err is %d", st)

        return coeff, err, kcoeff
    # rank 2
    elif signal.ndim == 2:
        # Compute biased autocorrelation up to lag = order
        bias    = 1. / sz
        corr    = bias * autocorr_oneside_nofft(signal, order, axis)

        if axis  % 2 == 0:
            # we transpose to have a major row order
            icorr   = corr.T
            icorr   = _N.require(icorr, requirements = 'C')
            at      = _N.zeros((icorr.shape[0], order+1), icorr.dtype)
            kt      = _N.zeros((icorr.shape[0], order), icorr.dtype)
            et      = _N.zeros(icorr.shape[0], icorr.dtype)

            if icorr.dtype == _N.float64:
                _lpc.dbl_levinson2d(icorr, icorr.shape[0], icorr.shape[1], at, et, kt)
            elif icorr.dtype == _N.float32:
                _lpc.flt_levinson2d(icorr, icorr.shape[0], icorr.shape[1], at, et, kt)
            else:
                raise TypeError("Only float32 and float64 supported")

            return at.T, et.T, kt.T
        elif axis % 2 == 1:
            icorr   = _N.require(corr, requirements = 'C')
            at      = _N.zeros((icorr.shape[0], order+1), icorr.dtype)
            kt      = _N.zeros((icorr.shape[0], order), icorr.dtype)
            et      = _N.zeros(icorr.shape[0], icorr.dtype)

            if icorr.dtype == _N.float64:
                _lpc.dbl_levinson2d(icorr, icorr.shape[0], icorr.shape[1], at, et, kt)
            elif icorr.dtype == _N.float32:
                _lpc.flt_levinson2d(icorr, icorr.shape[0], icorr.shape[1], at, et, kt)
            else:
                raise TypeError("Only float32 and float64 supported")

            return at, et, kt
        else:
            raise RuntimeError("This should not happen, this is a bug")
    else:
        raise RuntimeError("Sorry, only rank <= 2 supported for now")

def bench():
    size    = 256
    nframes = 4000
    order   = 24

    X       = _N.random.randn(nframes, size)
    X       = _N.require(X, requirements = 'C')

    niter   = 10

    # Contiguous
    print "Running optimized with ctypes"
    for i in range(niter):
        at, et, kt  = lpc2(X, order, axis = 1)

    a, e, k = _lpc2_py(X, order, axis = 1)
    _N.testing.assert_array_almost_equal(a, at, 10)
    _N.testing.assert_array_almost_equal(e, et, 10)
    _N.testing.assert_array_almost_equal(k, kt, 10)

    ## Non contiguous
    #print "Running optimized with ctypes (non contiguous)"
    #def ncontig(*args, **kargs):
    #    return lpc2(*args, **kargs)
    #X       = _N.require(X, requirements = 'F')
    #for i in range(niter):
    #    at, et, kt  = ncontig(X, order, axis = 1)

    #a, e, k = _lpc2_py(X, order, axis = 1)
    #_N.testing.assert_array_almost_equal(a, at, 10)
    #_N.testing.assert_array_almost_equal(e, et, 10)
    #_N.testing.assert_array_almost_equal(k, kt, 10)

    print "Benchmark func done"

if __name__ == '__main__':
    import hotshot, hotshot.stats
    profile_file    = 'lpc.prof'
    prof    = hotshot.Profile(profile_file, lineevents=1)
    prof.runcall(bench)
    p = hotshot.stats.load(profile_file)
    print p.sort_stats('cumulative').print_stats(20)
    prof.close()
