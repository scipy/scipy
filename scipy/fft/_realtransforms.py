import scipy.fftpack as _fftpack
from ._basic import _dispatch
from scipy._lib.uarray import Dispatchable
from scipy._lib.doccer import doc_replace
import numpy as np
import functools

__all__ = ['dct', 'idct', 'dst', 'idst', 'dctn', 'idctn', 'dstn', 'idstn']


@_dispatch
@doc_replace(_fftpack.dctn, 'fftpack', 'fft')
def dctn(x, type=2, shape=None, axes=None, norm=None, overwrite_x=False):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
@doc_replace(_fftpack.idctn, 'fftpack', 'fft')
def idctn(x, type=2, shape=None, axes=None, norm=None, overwrite_x=False):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
@doc_replace(_fftpack.dstn, 'fftpack', 'fft')
def dstn(x, type=2, shape=None, axes=None, norm=None, overwrite_x=False):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
@doc_replace(_fftpack.idstn, 'fftpack', 'fft')
def idstn(x, type=2, shape=None, axes=None, norm=None, overwrite_x=False):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
@doc_replace(_fftpack.dct, 'fftpack', 'fft')
def dct(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
@doc_replace(_fftpack.idct, 'fftpack', 'fft')
def idct(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
@doc_replace(_fftpack.dst, 'fftpack', 'fft')
def dst(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
@doc_replace(_fftpack.idst, 'fftpack', 'fft')
def idst(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False):
    return (Dispatchable(x, np.ndarray),)
