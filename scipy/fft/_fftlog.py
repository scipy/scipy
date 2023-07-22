import numpy as np
from scipy._lib._array_api import array_namespace
from . import _fftlog_multimethods as npfft_log_mm
from . import _fftlog_np as npfft_log

__all__ = [
    'fht', 'ifht',
    'fhtoffset',
]


def fht(a, dln, mu, offset=0.0, bias=0.0):
    xp = array_namespace(a)
    a = np.asarray(a)
    b = npfft_log_mm.fht(a, dln, mu, offset=offset, bias=bias)
    return xp.asarray(b)


def ifht(A, dln, mu, offset=0.0, bias=0.0):
    xp = array_namespace(A)
    A = np.asarray(A)
    B = npfft_log_mm.ifht(A, dln, mu, offset=offset, bias=bias)
    return xp.asarray(B)


def fhtcoeff(n, dln, mu, offset=0.0, bias=0.0):
    xp = array_namespace(n)
    n = np.asarray(n)
    m = npfft_log.fhtcoeff(n, dln, mu, offset=0.0, bias=0.0)
    return xp.asarray(m)


def fhtoffset(dln, mu, initial=0.0, bias=0.0):
    return npfft_log.fhtoffset(dln, mu, initial=initial, bias=bias)


def _fhtq(a, u, inverse=False):
    xp = array_namespace(a)
    a = np.asarray(a)
    b = npfft_log._fhtq(a, u, inverse=inverse)
    return xp.asarray(b)
