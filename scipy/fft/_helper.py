import numpy as np
from scipy._lib._array_api import array_namespace
from ._pocketfft import helper as _helper
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
    return nphelper._init_nd_shape_and_axes(x, shape, axes)


def fftfreq(n, d=1.0, xp=np, device=None):
    """
    Implements the Array API specification of fftfreq.
    """
    if hasattr(xp, 'fft'):
        return xp.fft.fftfreq(n, d, device=device)
    return np.fft.fftfreq(n, d, device=device)


def rfftfreq(n, d=1.0, xp=np, device=None):
    """
    Implements the Array API specification of rfftfreq.
    """
    if hasattr(xp, 'fft'):
        return xp.fft.rfftfreq(n, d, device=device)
    return np.fft.rfftfreq(n, d, device=device)


def fftshift(x, axes=None):
    """
    Implements the Array API specification of fftshift.
    """
    xp = array_namespace(x)
    if hasattr(xp, 'fft'):
        return xp.fft.fftshift(x, axes)
    x = np.asarray(x)
    y = np.fft.fftshift(x, axes)
    return xp.asarray(y)


def ifftshift(x, axes=None):
    """
    Implements the Array API specification of ifftshift.
    """
    xp = array_namespace(x)
    if hasattr(xp, 'fft'):
        return xp.fft.ifftshift(x, axes)
    x = np.asarray(x)
    y = np.fft.ifftshift(x, axes)
    return xp.asarray(y)
