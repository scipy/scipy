"""
Discrete Fourier Transforms - basic.py
"""
from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.fft._pocketfft import pypocketfft as pfft
from scipy.fft._fftpack.helper import _init_nd_shape_and_axes


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


def _asfarray(x):
    """Like numpy asfarray, except that it does not modify x dtype if x is
    already an array with a float dtype, and do not cast complex types to
    real."""
    if not hasattr(x, "dtype"):
        x = np.asarray(x)

    if x.dtype.char in np.typecodes["AllFloat"]:
        # 'dtype' attribute does not ensure that the
        # object is an ndarray (e.g. Series class
        # from the pandas library)
        if x.dtype == np.half:
            # no half-precision routines, so convert to single precision
            return np.asarray(x, dtype=np.float32)
        return np.asarray(x, dtype=x.dtype)

    return np.asfarray(x)


def _fix_shape(x, shape, axes):
    """Internal auxiliary function for _raw_fft, _raw_fftnd."""
    s = list(x.shape)
    for n, axis in zip(shape, axes):
        s[axis] = n

    must_copy = False

    # Build an nd slice with the dimensions to be read from x
    index = [slice(None)]*len(s)
    for n, ax in zip(shape, axes):
        if x.shape[ax] > n:
            index[ax] = slice(0, n)
        else:
            index[ax] = slice(0, x.shape[ax])
            must_copy = True

    index = tuple(index)

    if not must_copy:
        return x[index], False

    z = np.zeros(s, x.dtype)
    z[index] = x[index]
    return z, True


def _normalization_factor(shape, axes, norm, forward):
    """Returns the normalisation factor for a given nd-fft"""

    # TODO: This isn't passed through to pypocketfft as long double
    factor = np.longfloat(1.)
    if norm is None:
        if forward:
            return factor
        else:
            return factor / np.prod(shape, dtype=np.longfloat)

    if norm == 'ortho':
        return factor / np.sqrt(np.prod(shape, dtype=np.longfloat))

    raise ValueError(
        "Invalid norm value {}, should be None or \"ortho\".".format(norm))


def _1d_to_nd_shape_and_axes(n, axis):
    """Convert 1d arguments into nd form"""
    if n is None:
        shape = None
    else:
        shape = (n,)

    axes = (axis,)

    return shape, axes


def fft(x, n=None, axis=-1, norm=None, overwrite_x=False):
    """ Return discrete Fourier transform of real or complex sequence. """
    shape, axes = _1d_to_nd_shape_and_axes(n, axis)
    return fftn(x, shape, axes, norm, overwrite_x)


def ifft(x, n=None, axis=-1, norm=None, overwrite_x=False):
    """
    Return discrete inverse Fourier transform of real or complex sequence.
    """
    shape, axes = _1d_to_nd_shape_and_axes(n, axis)
    return ifftn(x, shape, axes, norm, overwrite_x)


def rfft(x, n=None, axis=-1, norm=None, overwrite_x=False):
    """
    Discrete Fourier transform of a real sequence.
    """
    shape, axes = _1d_to_nd_shape_and_axes(n, axis)
    return rfftn(x, shape, axes, norm, overwrite_x)


def irfft(x, n=None, axis=-1, norm=None, overwrite_x=False):
    """
    Return inverse discrete Fourier transform of real sequence x.
    """
    shape, axes = _1d_to_nd_shape_and_axes(n, axis)
    return irfftn(x, shape, axes, norm, overwrite_x)


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
    2-D dicsrete Fourier transform of a real sequence
    """
    return rfftn(x, shape, axes, norm, overwrite_x)


def irfft2(x, shape=None, axes=(-2,-1), norm=None, overwrite_x=False):
    """
    2-D dicsrete inverse Fourier transform of a real sequence
    """
    return irfftn(x, shape, axes, norm, overwrite_x)


def fftn(x, shape=None, axes=None, norm=None, overwrite_x=False):
    """
    Return multidimensional discrete Fourier transform.
    """
    tmp = _asfarray(x)

    # TODO: Optimize for real input.
    # First perform rfftn, then use hermitian symmetry to fill remaining output
    if np.isrealobj(tmp):
        tmp = tmp + 0.j

    shape, axes = _init_nd_shape_and_axes(tmp, shape, axes)
    overwrite_x = overwrite_x or _datacopied(tmp, x)

    # TODO: pocketfft raises here, should we?
    if len(axes) == 0:
        return x

    tmp, copied = _fix_shape(tmp, shape, axes)
    overwrite_x = overwrite_x or copied

    fct = _normalization_factor(shape, axes, norm, True)

    return pfft.fftn(tmp, axes, fct, overwrite_x, _default_workers)


def ifftn(x, shape=None, axes=None, norm=None, overwrite_x=False):
    """
    Return inverse multi-dimensional discrete Fourier transform.
    """
    tmp = _asfarray(x)

    # TODO: Optimize for hermitian signal?
    if np.isrealobj(tmp):
        tmp = tmp + 0.j

    shape, axes = _init_nd_shape_and_axes(tmp, shape, axes)
    overwrite_x = overwrite_x or _datacopied(tmp, x)

    if len(axes) == 0:
        return x

    tmp, copied = _fix_shape(tmp, shape, axes)
    overwrite_x = overwrite_x or copied

    fct = _normalization_factor(shape, axes, norm, False)

    return pfft.ifftn(tmp, axes, fct, overwrite_x, _default_workers)

def rfftn(x, shape=None, axes=None, norm=None, overwrite_x=False):
    """Return multi-dimentional discrete Fourier transform of real input"""
    tmp = _asfarray(x)

    shape, axes = _init_nd_shape_and_axes(tmp, shape, axes)
    tmp, _ = _fix_shape(tmp, shape, axes)
    fct = _normalization_factor(shape, axes, norm, True)

    if len(axes) == 0:
        return x

    # Note: overwrite_x is not utilised
    return pfft.rfftn(tmp, axes, fct, _default_workers)

def irfftn(x, shape=None, axes=None, norm=None, overwrite_x=False):
    """Multi-dimensional inverse discrete fourier transform with real output"""
    tmp = _asfarray(x)

    # TODO: Optimize for hermitian and real?
    if np.isrealobj(tmp):
        tmp = tmp + 0.j

    noshape = shape is None
    shape, axes = _init_nd_shape_and_axes(tmp, shape, axes)

    if len(axes) == 0:
        return x

    # TODO: defaulting to 2n - 1 may be a better choice (numpy/numpy#13357)
    if noshape:
        shape[-1] = (x.shape[axes[-1]] - 1) * 2

    fct = _normalization_factor(shape, axes, norm, False)

    # Last axis utilises hermitian symmetry
    lastsize = shape[-1]
    shape[-1] = (shape[-1] // 2) + 1

    tmp, _ = _fix_shape(tmp, shape, axes)

    # Note: overwrite_x is not utilised
    return pfft.irfftn(tmp, axes, lastsize, fct, _default_workers)
