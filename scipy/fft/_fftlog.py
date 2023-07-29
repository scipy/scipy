import numpy as np
from scipy._lib._array_api import array_namespace
from . import _fftlog_multimethods as npfft_log_mm

__all__ = [
    'fht', 'ifht'
]


def fht(a, dln, mu, offset=0.0, bias=0.0):
    xp = array_namespace(a)
    a = np.asarray(a)
    b = npfft_log_mm.fht(a, dln, mu, offset=offset, bias=bias)
    return xp.asarray(b.copy())


def ifht(A, dln, mu, offset=0.0, bias=0.0):
    xp = array_namespace(A)
    A = np.asarray(A)
    B = npfft_log_mm.ifht(A, dln, mu, offset=offset, bias=bias)
    return xp.asarray(B.copy())
