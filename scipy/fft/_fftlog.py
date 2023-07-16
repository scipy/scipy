import numpy as np
from warnings import warn
from ._basic import rfft, irfft
from ..special import loggamma, poch
from . import _fftlog_np as npfft_log

__all__ = [
    'fht', 'ifht',
    'fhtoffset',
]


def fht(a, dln, mu, offset=0.0, bias=0.0):
    return npfft_log.fht(a, dln, mu, offset=offset, bias=bias)


def ifht(A, dln, mu, offset=0.0, bias=0.0):
    return npfft_log.ifht(A, dln, mu, offset=offset, bias=bias)


def fhtcoeff(n, dln, mu, offset=0.0, bias=0.0):
    return npfft_log.fhtcoeff(n, dln, mu, offset=offset, bias=bias)


def fhtoffset(dln, mu, initial=0.0, bias=0.0):
    return npfft_log.fhtoffset(dln, mu, initial=initial, bias=bias)


def _fhtq(a, u, inverse=False):
    return npfft_log._fhtq(a, u, inverse=inverse)
