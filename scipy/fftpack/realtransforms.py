"""
Real spectrum tranforms (DCT, DST, MDCT)
"""

__all__ = ['dct1', 'dct2']

import numpy as np
from scipy.fftpack import _fftpack

import atexit
atexit.register(_fftpack.destroy_dct1_cache)
atexit.register(_fftpack.destroy_dct2_cache)

def dct1(x, n=None):
    """
    Return Discrete Cosine Transform (type I) of arbitrary type sequence x.

    Parameters
    ----------
    x : array-like
        input array.
    n : int, optional
        Length of the transform.

    Returns
    -------
    y : real ndarray
    """
    return _dct(x, 1, n)

def dct2(x, n=None):
    """
    Return Discrete Cosine Transform (type II) of arbitrary type sequence x.

    Parameters
    ----------
    x : array-like
        input array.
    n : int, optional
        Length of the transform.

    Returns
    -------
    y : real ndarray
    """
    return _dct(x, 2, n)

def _dct(x, type, n=None, axis=-1, overwrite_x=0):
    """
    Return Discrete Cosine Transform of arbitrary type sequence x.

    Parameters
    ----------
    x : array-like
        input array.
    n : int, optional
        Length of the transform.
    axis : int, optional
        Axis along which the dct is computed. (default=-1)
    overwrite_x : bool, optional
        If True the contents of x can be destroyed. (default=False)

    Returns
    -------
    z : real ndarray

    """
    tmp = np.asarray(x)
    if not np.isrealobj(tmp):
        raise TypeError,"1st argument must be real sequence"

    if n is None:
        n = tmp.shape[axis]
    else:
        raise NotImplemented("Padding/truncating not yet implemented")

    if type == 1:
        f = _fftpack.dct1
    elif type == 2:
        f = _fftpack.dct2
    else:
        raise ValueError("Type %d not understood" % type)

    if axis == -1 or axis == len(tmp.shape) - 1:
        return f(tmp, n, 0, overwrite_x)
    else:
        raise NotImplementedError("Axis arg not yet implemented")

    #tmp = swapaxes(tmp, axis, -1)
    #tmp = work_function(tmp,n,1,0,overwrite_x)
    #return swapaxes(tmp, axis, -1)
