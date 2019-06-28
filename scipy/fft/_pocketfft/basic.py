"""
Discrete Fourier Transforms - basic.py
"""
from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.fft._pocketfft import pypocketfft as pfft
from scipy.fftpack.helper import _init_nd_shape_and_axes


# TODO: Build with OpenMp and add configuration support
_default_workers = 1

def _datacopied(arr, original):
    """
    Strict check for `arr` not sharing any data with `original`,
    under the assumption that arr = asarray(original)
    """
    if arr is original:
        return False
    if not isinstance(original, np.ndarray) and hasattr(original, '__array__'):
        return False
    return arr.base is None


def _fix_shape(x, shape, axes):
    """Internal auxiliary function for _raw_fft, _raw_fftnd."""
    must_copy = False

    # Build an nd slice with the dimensions to be read from x
    index = [slice(None)]*x.ndim
    for n, ax in zip(shape, axes):
        if x.shape[ax] >= n:
            index[ax] = slice(0, n)
        else:
            index[ax] = slice(0, x.shape[ax])
            must_copy = True

    index = tuple(index)

    if not must_copy:
        return x[index], False

    s = list(x.shape)
    for n, axis in zip(shape, axes):
        s[axis] = n

    z = np.zeros(s, x.dtype)
    z[index] = x[index]
    return z, True


def _fix_shape_1d(x, n, axis):
    if n < 1:
        raise ValueError(
            "invalid number of data points ({0}) specified".format(n))

    return _fix_shape(x, (n,), (axis,))


def _normalization(norm, forward):
    """Returns the pypocketfft normalization mode from the norm argument"""

    if norm is None:
        return pfft.norm_t.none if forward else pfft.norm_t.size

    if norm == 'ortho':
        return pfft.norm_t.ortho

    raise ValueError(
        "Invalid norm value {}, should be None or \"ortho\".".format(norm))


def fft(x, n=None, axis=-1, norm=None, overwrite_x=False):
    """ Return discrete Fourier transform of real or complex sequence. """
    tmp = np.asarray(x)
    overwrite_x = overwrite_x or _datacopied(tmp, x)
    norm = _normalization(norm, True)

    if n is not None:
        tmp, copied = _fix_shape_1d(tmp, n, axis)
        overwrite_x = overwrite_x or copied
    elif tmp.shape[axis] < 1:
        raise ValueError("invalid number of data points ({0}) specified"
                         .format(tmp.shape[axis]))

    return pfft.fftn(tmp, (axis,), norm, overwrite_x, _default_workers)


def ifft(x, n=None, axis=-1, norm=None, overwrite_x=False):
    """
    Return discrete inverse Fourier transform of real or complex sequence.
    """
    tmp = np.asarray(x)
    overwrite_x = overwrite_x or _datacopied(tmp, x)
    norm = _normalization(norm, False)

    if n is not None:
        tmp, copied = _fix_shape_1d(tmp, n, axis)
        overwrite_x = overwrite_x or copied
    elif tmp.shape[axis] < 1:
        raise ValueError("invalid number of data points ({0}) specified"
                         .format(tmp.shape[axis]))

    return pfft.ifftn(tmp, (axis,), norm, overwrite_x, _default_workers)


def rfft(x, n=None, axis=-1, norm=None, overwrite_x=False):
    """
    Discrete Fourier transform of a real sequence.
    """
    tmp = np.asarray(x)
    norm = _normalization(norm, True)

    if not np.isrealobj(tmp):
        raise TypeError("x must be a real sequence")

    if n is not None:
        tmp, _ = _fix_shape_1d(tmp, n, axis)
    elif tmp.shape[axis] < 1:
        raise ValueError("invalid number of data points ({0}) specified"
                         .format(tmp.shape[axis]))

    # Note: overwrite_x is not utilised
    return pfft.rfftn(tmp, (axis,), norm, _default_workers)


def irfft(x, n=None, axis=-1, norm=None, overwrite_x=False):
    """
    Return inverse discrete Fourier transform of real sequence x.
    """
    tmp = np.asarray(x)
    norm = _normalization(norm, False)

    # TODO: Optimize for hermitian and real?
    if np.isrealobj(tmp):
        tmp = tmp + 0.j

    # Last axis utilises hermitian symmetry
    if n is None:
        n = (tmp.shape[axis] - 1) * 2
        if n < 1:
            raise ValueError("Invalid number of data points ({0}) specified"
                             .format(n))
    else:
        tmp, _ = _fix_shape_1d(tmp, (n//2) + 1, axis)

    # Note: overwrite_x is not utilised
    return pfft.irfftn(tmp, (axis,), n, norm, _default_workers)


def fft2(x, shape=None, axes=(-2,-1), norm=None, overwrite_x=False):
    """
    2-D discrete Fourier transform.
    """
    return fftn(x, shape, axes, norm, overwrite_x)


def ifft2(x, shape=None, axes=(-2,-1), norm=None, overwrite_x=False):
    """
    2-D discrete inverse Fourier transform of real or complex sequence.
    """
    return ifftn(x, shape, axes, norm, overwrite_x)


def rfft2(x, shape=None, axes=(-2,-1), norm=None, overwrite_x=False):
    """
    2-D discrete Fourier transform of a real sequence
    """
    return rfftn(x, shape, axes, norm, overwrite_x)


def irfft2(x, shape=None, axes=(-2,-1), norm=None, overwrite_x=False):
    """
    2-D discrete inverse Fourier transform of a real sequence
    """
    return irfftn(x, shape, axes, norm, overwrite_x)


def fftn(x, shape=None, axes=None, norm=None, overwrite_x=False):
    """
    Return multidimensional discrete Fourier transform.
    """
    tmp = np.asarray(x)

    shape, axes = _init_nd_shape_and_axes(tmp, shape, axes)
    overwrite_x = overwrite_x or _datacopied(tmp, x)

    if len(axes) == 0:
        return x

    tmp, copied = _fix_shape(tmp, shape, axes)
    overwrite_x = overwrite_x or copied

    norm = _normalization(norm, True)

    return pfft.fftn(tmp, axes, norm, overwrite_x, _default_workers)


def ifftn(x, shape=None, axes=None, norm=None, overwrite_x=False):
    """
    Return inverse multi-dimensional discrete Fourier transform.
    """
    tmp = np.asarray(x)

    shape, axes = _init_nd_shape_and_axes(tmp, shape, axes)
    overwrite_x = overwrite_x or _datacopied(tmp, x)

    if len(axes) == 0:
        return x

    tmp, copied = _fix_shape(tmp, shape, axes)
    overwrite_x = overwrite_x or copied

    norm = _normalization(norm, False)

    return pfft.ifftn(tmp, axes, norm, overwrite_x, _default_workers)

def rfftn(x, shape=None, axes=None, norm=None, overwrite_x=False):
    """Return multi-dimensional discrete Fourier transform of real input"""
    tmp = np.asarray(x)

    if not np.isrealobj(tmp):
        raise TypeError("x must be a real sequence")

    shape, axes = _init_nd_shape_and_axes(tmp, shape, axes)
    tmp, _ = _fix_shape(tmp, shape, axes)
    norm = _normalization(norm, True)

    if len(axes) == 0:
        raise ValueError("at least 1 axis must be transformed")

    # Note: overwrite_x is not utilised
    return pfft.rfftn(tmp, axes, norm, _default_workers)

def irfftn(x, shape=None, axes=None, norm=None, overwrite_x=False):
    """Multi-dimensional inverse discrete fourier transform with real output"""
    tmp = np.asarray(x)

    # TODO: Optimize for hermitian and real?
    if np.isrealobj(tmp):
        tmp = tmp + 0.j

    noshape = shape is None
    shape, axes = _init_nd_shape_and_axes(tmp, shape, axes)

    if len(axes) == 0:
        raise ValueError("at least 1 axis must be transformed")

    if noshape:
        shape[-1] = (x.shape[axes[-1]] - 1) * 2

    norm = _normalization(norm, False)

    # Last axis utilises hermitian symmetry
    lastsize = shape[-1]
    shape[-1] = (shape[-1] // 2) + 1

    tmp, _ = _fix_shape(tmp, shape, axes)

    # Note: overwrite_x is not utilised
    return pfft.irfftn(tmp, axes, lastsize, norm, _default_workers)
