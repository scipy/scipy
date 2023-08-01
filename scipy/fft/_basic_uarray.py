from scipy._lib.uarray import generate_multimethod, Dispatchable
import numpy as np


def _x_replacer(args, kwargs, dispatchables):
    """
    uarray argument replacer to replace the transform input array (``x``)
    """
    if len(args) > 0:
        return (dispatchables[0],) + args[1:], kwargs
    kw = kwargs.copy()
    kw['x'] = dispatchables[0]
    return args, kw


def _dispatch(func):
    """
    Function annotation that creates a uarray multimethod from the function
    """
    return generate_multimethod(func, _x_replacer, domain="numpy.scipy.fft")


@_dispatch
def fft(x, n=None, axis=-1, norm=None, overwrite_x=False, workers=None, *,
        plan=None):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
def ifft(x, n=None, axis=-1, norm=None, overwrite_x=False, workers=None, *,
         plan=None):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
def rfft(x, n=None, axis=-1, norm=None, overwrite_x=False, workers=None, *,
         plan=None):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
def irfft(x, n=None, axis=-1, norm=None, overwrite_x=False, workers=None, *,
          plan=None):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
def hfft(x, n=None, axis=-1, norm=None, overwrite_x=False, workers=None, *,
         plan=None):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
def ihfft(x, n=None, axis=-1, norm=None, overwrite_x=False, workers=None, *,
          plan=None):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
def fftn(x, s=None, axes=None, norm=None, overwrite_x=False, workers=None, *,
         plan=None):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
def ifftn(x, s=None, axes=None, norm=None, overwrite_x=False, workers=None, *,
          plan=None):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
def fft2(x, s=None, axes=(-2, -1), norm=None, overwrite_x=False, workers=None, *,
         plan=None):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
def ifft2(x, s=None, axes=(-2, -1), norm=None, overwrite_x=False, workers=None, *,
          plan=None):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
def rfftn(x, s=None, axes=None, norm=None, overwrite_x=False, workers=None, *,
          plan=None):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
def rfft2(x, s=None, axes=(-2, -1), norm=None, overwrite_x=False, workers=None, *,
          plan=None):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
def irfftn(x, s=None, axes=None, norm=None, overwrite_x=False, workers=None, *,
           plan=None):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
def irfft2(x, s=None, axes=(-2, -1), norm=None, overwrite_x=False, workers=None, *,
           plan=None):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
def hfftn(x, s=None, axes=None, norm=None, overwrite_x=False, workers=None, *,
          plan=None):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
def hfft2(x, s=None, axes=(-2, -1), norm=None, overwrite_x=False, workers=None, *,
          plan=None):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
def ihfftn(x, s=None, axes=None, norm=None, overwrite_x=False, workers=None, *,
           plan=None):
    return (Dispatchable(x, np.ndarray),)


@_dispatch
def ihfft2(x, s=None, axes=(-2, -1), norm=None, overwrite_x=False, workers=None, *,
           plan=None):
    return (Dispatchable(x, np.ndarray),)
