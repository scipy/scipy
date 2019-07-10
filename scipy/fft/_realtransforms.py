import scipy.fftpack as _fftpack
from ._basic import _dispatch
from scipy._lib.uarray import Dispatchable
import numpy as np
import functools

__all__ = ['dct', 'idct', 'dst', 'idst', 'dctn', 'idctn', 'dstn', 'idstn']

def _inherit_doc(fftpack_fun):
    def inner(fft_fun):
        doc = fftpack_fun.__doc__ or ''
        fft_fun.__doc__ = doc.replace('fftpack', 'fft')
        return fft_fun

    return inner


@_dispatch
@_inherit_doc(_fftpack.dctn)
def dctn(x, type=2, shape=None, axes=None, norm=None, overwrite_x=False):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
@_inherit_doc(_fftpack.idctn)
def idctn(x, type=2, shape=None, axes=None, norm=None, overwrite_x=False):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
@_inherit_doc(_fftpack.dstn)
def dstn(x, type=2, shape=None, axes=None, norm=None, overwrite_x=False):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
@_inherit_doc(_fftpack.idstn)
def idstn(x, type=2, shape=None, axes=None, norm=None, overwrite_x=False):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
@_inherit_doc(_fftpack.dct)
def dct(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
@_inherit_doc(_fftpack.idct)
def idct(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
@_inherit_doc(_fftpack.dst)
def dst(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
@_inherit_doc(_fftpack.idst)
def idst(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False):
    return (Dispatchable(x, np.ndarray),)
