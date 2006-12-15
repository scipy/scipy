#! /usr/bin/env python
# Last Change: Fri Dec 15 10:00 PM 2006 J

# TODO: - proper test
# TODO: - proper profiling

from numpy.fft import fft, ifft
from numpy import correlate, log2, floor, conj, real, \
        concatenate, sum, max

from warnings import warn

# use ctype to have one sided c imp of autocorr
import ctypes
from ctypes import c_uint, c_int
from numpy.ctypeslib import ndpointer, load_library

ctypes_major    = int(ctypes.__version__.split('.')[0])
if ctypes_major < 1:
    msg =  "version of ctypes is %s, expected at least %s" \
            % (ctypes.__version__, '1.0.1')
    raise importerror(msg)

import numpy as N

# load autocorr lib
_autocorr   = load_library('gabsig.so', __file__)

#===============================
# define the functions with args
#===============================

# contiguous 1d
arg1    = ndpointer(dtype = N.float64, flags='CONTIGUOUS,ALIGNED')
arg2    = c_uint
arg3    = ndpointer(dtype = N.float64, flags='CONTIGUOUS,ALIGNED')
arg4    = c_uint
_autocorr.dbl_xcorr_nofft_1d.argtypes  = [arg1, arg2, arg3, arg4]
_autocorr.dbl_xcorr_nofft_1d.restype   = c_int

arg1    = ndpointer(dtype = N.float32, flags='CONTIGUOUS,ALIGNED')
arg2    = c_uint
arg3    = ndpointer(dtype = N.float32, flags='CONTIGUOUS,ALIGNED')
arg4    = c_uint
_autocorr.flt_xcorr_nofft_1d.argtypes  = [arg1, arg2, arg3, arg4]
_autocorr.flt_xcorr_nofft_1d.restype   = c_int

# non contiguous 1d
arg1    = ndpointer(dtype = N.float64, flags = 'ALIGNED')
arg2    = c_uint
arg3    = c_uint
arg4    = ndpointer(dtype = N.float64, flags = 'ALIGNED')
arg5    = c_uint
arg6    = c_uint
_autocorr.dbl_xcorr_nofft_1d_noncontiguous.argtypes  = [arg1, \
        arg2, arg3, arg4, arg5, arg6]
_autocorr.dbl_xcorr_nofft_1d_noncontiguous.restype   = c_int

arg1    = ndpointer(dtype = N.float32, flags = 'ALIGNED')
arg2    = c_uint
arg3    = c_uint
arg4    = ndpointer(dtype = N.float32, flags = 'ALIGNED')
arg5    = c_uint
arg6    = c_uint
_autocorr.flt_xcorr_nofft_1d_noncontiguous.argtypes  = [arg1, \
        arg2, arg3, arg4, arg5, arg6]
_autocorr.flt_xcorr_nofft_1d_noncontiguous.restype   = c_int

# contiguous 2d
arg1    = ndpointer(dtype = N.float64, flags='ALIGNED')
arg2    = c_uint
arg3    = c_uint
arg4    = ndpointer(dtype = N.float64, flags='ALIGNED')
arg5    = c_uint
_autocorr.dbl_xcorr_nofft_2d.argtypes  = [arg1, arg2, arg3, arg4, arg5]
_autocorr.dbl_xcorr_nofft_2d.restype   = c_int

arg1    = ndpointer(dtype = N.float32, flags='ALIGNED')
arg2    = c_uint
arg3    = c_uint
arg4    = ndpointer(dtype = N.float32, flags='ALIGNED')
arg5    = c_uint
_autocorr.flt_xcorr_nofft_2d.argtypes  = [arg1, arg2, arg3, arg4, arg5]
_autocorr.flt_xcorr_nofft_2d.restype   = c_int

# non contiguous 2d
arg1    = ndpointer(dtype = N.float64, flags='ALIGNED')
arg2    = c_uint
arg3    = c_uint
arg4    = c_uint
arg5    = c_uint
arg6    = ndpointer(dtype = N.float64, flags='ALIGNED')
arg7    = c_uint
arg8    = c_uint
arg9    = c_uint
_autocorr.dbl_xcorr_nofft_2d_noncontiguous.argtypes  = [arg1, arg2, \
        arg3, arg4, arg5, arg6, arg7, arg8, arg9]
_autocorr.dbl_xcorr_nofft_2d_noncontiguous.restype   = c_int

arg1    = ndpointer(dtype = N.float32, flags='ALIGNED')
arg2    = c_uint
arg3    = c_uint
arg4    = c_uint
arg5    = c_uint
arg6    = ndpointer(dtype = N.float32, flags='ALIGNED')
arg7    = c_uint
arg8    = c_uint
arg9    = c_uint
_autocorr.flt_xcorr_nofft_2d_noncontiguous.argtypes  = [arg1, arg2, \
        arg3, arg4, arg5, arg6, arg7, arg8, arg9]
_autocorr.flt_xcorr_nofft_2d_noncontiguous.restype   = c_int

#======================================
# Fonctions to be used for testing only
#======================================
def _raw_autocorr_1d(signal, lag):
    assert signal.ndim == 1
    assert signal.flags['CONTIGUOUS']
    
    if lag >= signal.size:
        raise RuntimeError("lag should be < to input size")

    if signal.dtype == N.float64:
        res = N.zeros((lag+1), N.float64)
        _autocorr.dbl_xcorr_nofft_1d(signal, signal.size, res, lag) 
    elif signal.dtype == N.float32:
        res = N.zeros((lag+1), N.float32)
        _autocorr.flt_xcorr_nofft_1d(signal, signal.size, res, lag) 
    else:
        raise TypeError("only float 32 and 64 bits supported for now")

    return res

def _raw_autocorr_1d_noncontiguous(signal, lag):
    assert signal.ndim == 1
    
    if lag >= signal.size:
        raise RuntimeError("lag should be < to input size")

    if signal.dtype == N.float64:
        res = N.zeros((lag+1), N.float64)
        _autocorr.dbl_xcorr_nofft_1d_noncontiguous(signal, signal.size, 
                signal.strides[0], res, res.strides[0], lag) 
    elif signal.dtype == N.float32:
        res = N.zeros((lag+1), N.float32)
        _autocorr.flt_xcorr_nofft_1d_noncontiguous(signal, signal.size, 
                signal.strides[0], res, res.strides[0], lag) 
    else:
        raise TypeError("only float 32 and 64 bits supported for now")

    return res

# python implementation of autocorr for rank <= 2
def _autocorr_oneside_nofft_py(signal, lag, axis = -1):
    if signal.ndim > 2:
        raise NotImplemented("only for rank <=2")
    
    if axis  % 2 == 0:
        res     = N.zeros((lag+1, signal.shape[1]), signal.dtype)
        center  = signal.shape[0] - 1
        for i in range(signal.shape[1]):
            #print "compute corr of " + str(signal[:, i])
            res[:, i]   = correlate(signal[:, i], signal[:, i], \
                    'full')[center:center+lag+1]
    elif axis % 2 == 1:
        res     = N.zeros((signal.shape[0], lag+1), signal.dtype)
        center  = signal.shape[1] - 1
        for i in range(signal.shape[0]):
            #print "compute corr of " + str(signal[i])
            res[i]  = correlate(signal[i], signal[i], \
                    'full')[center:center+lag+1]
    else:
        raise RuntimeError("this should bnot happen, please fill a bug")

    return res

#=============
# Public API
#=============
def autocorr_oneside_nofft(signal, lag, axis = -1):
    """Compute the righ side of autocorrelation along the axis, for lags up to lag.
    
    This implementation does NOT use FFT."""
    # TODO  For rank < 2, the overhead of python code may be significant. Should
    # TODO not be difficult to do in C anyway (we can still use ctypes)

    # rank 0, 1
    if signal.ndim < 2:
        size    = signal.shape[-1]
        if lag >= size:
            raise RuntimeError("lag should be < to input size")

        res = N.zeros((lag+1), signal.dtype)

        if signal.flags['CONTIGUOUS']:
            if signal.dtype == N.float64:
                _autocorr.dbl_xcorr_nofft_1d(signal, size, res, lag) 
            elif signal.dtype == N.float32:
                _autocorr.flt_xcorr_nofft_1d(signal, size, res, lag) 
            else:
                raise TypeError("only float 32 and 64 bits supported for now")
        else:
            istride = signal.strides[0]
            ostride = signal.itemsize
            if signal.dtype == N.float64:
                _autocorr.dbl_xcorr_nofft_1d_noncontiguous(signal, size, istride, 
                        res, ostride,  lag) 
            elif signal.dtype == N.float32:
                _autocorr.flt_xcorr_nofft_1d_noncontiguous(signal, size, istride, 
                        res, ostride,  lag) 
            else:
                raise TypeError("only float 32 and 64 bits supported for now")

    # rank 2 case 
    elif signal.ndim == 2:
        size    = signal.shape[axis]
        if lag >= size:
            raise RuntimeError("lag should be < to input size")
            res = N.zeros((signal.shape[0], lag+1), signal.dtype)
        else:
            res = N.zeros((lag+1, signal.shape[1]), signal.dtype)

        if signal.dtype == N.float64:
            # contiguous case
            if signal.flags['C'] and axis % 2 == 1:
                res = N.zeros((signal.shape[0], lag+1), signal.dtype)
                _autocorr.dbl_xcorr_nofft_2d(signal, signal.shape[0], signal.shape[1], 
                        res, lag) 
            # contiguous case
            elif signal.flags['F'] and axis % 2 == 0:
                res = N.zeros((lag+1, signal.shape[1]), signal.dtype, order = 'F')
                _autocorr.dbl_xcorr_nofft_2d(signal, signal.shape[1], signal.shape[0], 
                        res, lag) 
            # non contiguous case
            elif axis % 2 == 0:
                res = N.zeros((lag+1, signal.shape[1]), signal.dtype)
                warn("non contiguous used, this will be slow")
                _autocorr.dbl_xcorr_nofft_2d_noncontiguous(signal, 
                        signal.shape[1], signal.shape[0], 
                        signal.strides[1], signal.strides[0],
                        res, res.strides[1], res.strides[0], lag) 
            elif axis % 2 == 1:
                res = N.zeros((signal.shape[0], lag+1), signal.dtype)
                warn("non contiguous used, this will be slow")
                _autocorr.dbl_xcorr_nofft_2d_noncontiguous(signal, 
                        signal.shape[0], signal.shape[1], 
                        signal.strides[0], signal.strides[1],
                        res, res.strides[0], res.strides[1], lag) 
        elif signal.dtype == N.float32:
            # contiguous case
            if signal.flags['C'] and axis % 2 == 1:
                res = N.zeros((signal.shape[0], lag+1), signal.dtype)
                _autocorr.flt_xcorr_nofft_2d(signal, signal.shape[0], signal.shape[1], 
                        res, lag) 
            # contiguous case
            elif signal.flags['F'] and axis % 2 == 0:
                res = N.zeros((lag+1, signal.shape[1]), signal.dtype, order = 'F')
                _autocorr.flt_xcorr_nofft_2d(signal, signal.shape[1], signal.shape[0], 
                        res, lag) 
            # non contiguous case
            elif axis % 2 == 0:
                res = N.zeros((lag+1, signal.shape[1]), signal.dtype)
                warn("non contiguous used, this will be slow")
                _autocorr.flt_xcorr_nofft_2d_noncontiguous(signal, 
                        signal.shape[1], signal.shape[0], 
                        signal.strides[1], signal.strides[0],
                        res, res.strides[1], res.strides[0], lag) 
            elif axis % 2 == 1:
                res = N.zeros((signal.shape[0], lag+1), signal.dtype)
                warn("non contiguous used, this will be slow")
                _autocorr.flt_xcorr_nofft_2d_noncontiguous(signal, 
                        signal.shape[0], signal.shape[1], 
                        signal.strides[0], signal.strides[1],
                        res, res.strides[0], res.strides[1], lag) 
        else:
            raise TypeError("only float 32 and 64 bits supported for now")
    else:
        raise RuntimeError("rank > 2 not supported yet")

    return res

def nextpow2(n):
    """Returns p such as 2 ** p >= n """
    p   = N.floor(N.log2(n))
    if 2 **  p ==  n:
        return p
    else:
        return p + 1

def autocorr_fft(signal, axis = -1):
    """Return full autocorrelation along specified axis. Use fft
    for computation."""
    if N.ndim(signal) == 0:
        return signal
    elif signal.ndim == 1:
        n       = signal.shape[0]
        nfft    = int(2 ** nextpow2(2 * n - 1))
        lag     = n - 1
        a       = fft(signal, n = nfft, axis = -1)
        au      = ifft(a * N.conj(a), n = nfft, axis = -1)
        return N.require(N.concatenate((au[-lag:], au[:lag+1])), dtype = signal.dtype)
    elif signal.ndim == 2:
        n       = signal.shape[axis]
        lag     = n - 1
        nfft    = int(2 ** nextpow2(2 * n - 1))
        a       = fft(signal, n = nfft, axis = axis)
        au      = ifft(a * N.conj(a), n = nfft, axis = axis)
        if axis == 0:
            return N.require(N.concatenate( (au[-lag:], au[:lag+1]), axis = axis), \
                    dtype = signal.dtype)
        else:
            return N.require(N.concatenate( (au[:, -lag:], au[:, :lag+1]), 
                        axis = axis), dtype = signal.dtype)
    else:
        raise RuntimeError("rank >2 not supported yet")
        
def bench():
    size    = 256
    nframes = 4000
    lag     = 24

    X       = N.random.randn(nframes, size)
    X       = N.require(X, requirements = 'C')

    niter   = 10

    # Contiguous
    print "Running optimized with ctypes"
    def contig(*args, **kargs):
        return autocorr_oneside_nofft(*args, **kargs)
    for i in range(niter):
        Yt  = contig(X, lag, axis = 1)

    Yr  = _autocorr_oneside_nofft_py(X, lag, axis = 1)
    N.testing.assert_array_almost_equal(Yt, Yr, 10)

    # Non contiguous
    print "Running optimized with ctypes (non contiguous)"
    def ncontig(*args, **kargs):
        return autocorr_oneside_nofft(*args, **kargs)
    X       = N.require(X, requirements = 'F')
    for i in range(niter):
        Yt  = ncontig(X, lag, axis = 1)

    Yr  = _autocorr_oneside_nofft_py(X, lag, axis = 1)
    N.testing.assert_array_almost_equal(Yt, Yr, 10)

    print "Benchmark func done"

if __name__ == '__main__':
    import hotshot, hotshot.stats
    profile_file    = 'autocorr.prof'
    prof    = hotshot.Profile(profile_file, lineevents=1)
    prof.runcall(bench)
    p = hotshot.stats.load(profile_file)
    print p.sort_stats('cumulative').print_stats(20)
    prof.close()
