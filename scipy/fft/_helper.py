import numpy as np
from scipy._lib._array_api import array_namespace
from . import _helper_np as nphelper


def next_fast_len(target, real=False):
    """
    See the documentation in _helper_np.py.
    """
    return nphelper.next_fast_len(target, real=real)


def _init_nd_shape_and_axes(x, shape, axes):
    """
    See the documentation in _helper_np.py.
    """
    xp = array_namespace(x)
    x = np.asarray(x)
    y = nphelper._init_nd_shape_and_axes(x, shape, axes)
    return xp.asarray(y)


def fftfreq(n, d=1.0, *, xp=np, device=None):
    """
    Implements the array API specification of fftfreq.
    """
    # device is currently not supported in numpy, cupy or array-api-compat
    if device is not None:
        raise ValueError(
            "Providing the 'device' parameter is not yet supported.")
    if hasattr(xp, 'fft'):
        return xp.fft.fftfreq(n, d=d)
    return np.fft.fftfreq(n, d=d)


def rfftfreq(n, d=1.0, *, xp=np, device=None):
    """
    Implements the array API specification of rfftfreq.
    """
    # device is currently not supported in numpy, cupy or array-api-compat
    if device is not None:
        raise ValueError(
            "Providing the 'device' parameter is not yet supported.")
    if hasattr(xp, 'fft'):
        return xp.fft.rfftfreq(n, d=d)
    return np.fft.rfftfreq(n, d=d)


def fftshift(x, axes=None):
    """
    Implements the array API specification of fftshift.
    """
    xp = array_namespace(x)
    if hasattr(xp, 'fft'):
        return xp.fft.fftshift(x, axes=axes)
    x = np.asarray(x)
    y = np.fft.fftshift(x, axes=axes)
    return xp.asarray(y)


def ifftshift(x, axes=None):
    """
    Implements the array API specification of ifftshift.
    """
    xp = array_namespace(x)
    if hasattr(xp, 'fft'):
        return xp.fft.ifftshift(x, axes=axes)
    x = np.asarray(x)
    y = np.fft.ifftshift(x, axes=axes)
    return xp.asarray(y)
