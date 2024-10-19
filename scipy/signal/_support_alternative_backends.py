import sys
import functools
from scipy._lib._array_api import (
    array_namespace, is_cupy, is_jax, scipy_namespace_for, SCIPY_ARRAY_API
)
from ._signaltools import (convolve, fftconvolve, convolve2d, oaconvolve,
                           correlate, correlate2d, order_filter, medfilt, medfilt2d,
                           wiener, detrend, hilbert, hilbert2, lfilter, deconvolve,
                           sosfilt, sosfiltfilt, sosfilt_zi, lfilter_zi,
                           filtfilt,)

MODULE_NAME = 'signal'

# jax.scipy.signal has only partial coverage of scipy.signal, so we keep the list
# of functions we can delegate to JAX
# https://jax.readthedocs.io/en/latest/jax.scipy.html
JAX_SIGNAL_FUNCS = [
    'fftconvolve', 'convolve', 'convolve2d', 'correlate', 'correlate2d',
    'csd', 'detrend', 'istft', 'welch'
]

# some cupyx.scipy.signal functions are incompatible with their scipy counterparts
CUPY_BLACKLIST = ['lfilter_zi', 'sosfilt_zi']

def delegate_xp(delegator, module_name):
    def inner(func):
        @functools.wraps(func)
        def wrapper(*args, **kwds):
            try:
                xp = delegator(*args, **kwds)
            except TypeError:
                # object arrays
                import numpy as np
                xp = np

            # try delegating to a cupyx/jax namesake
            if is_cupy(xp) and func.__name__ not in CUPY_BLACKLIST:
                # https://github.com/cupy/cupy/issues/8336
                import importlib
                cupyx_module = importlib.import_module(f"cupyx.scipy.{module_name}")
                cupyx_func = getattr(cupyx_module, func.__name__)
                return cupyx_func(*args, **kwds)
            elif is_jax(xp) and func.__name__ in JAX_SIGNAL_FUNCS:
                spx = scipy_namespace_for(xp)
                jax_module = getattr(spx, module_name)
                jax_func = getattr(jax_module, func.__name__)
                return jax_func(*args, **kwds)
            else:
                # the original function
                return func(*args, **kwds)
        return wrapper
    return inner



# X_signature signature must match the signature of X

def convolve_signature(in1, in2, *args, **kwds):
    xp = array_namespace(in1, in2)
    return xp

fftconvolve_signature = convolve_signature
oaconvolve_signature = convolve_signature
correlate_signature = convolve_signature
correlate_signature = convolve_signature
convolve2d_signature = convolve_signature
correlate2d_signature = convolve_signature


def medfilt_signature(volume, kernel_size=None):
    xp = array_namespace(volume)
    return xp


def medfilt2d_signature(input, kernel_size=3):
    xp = array_namespace(input)
    return xp


def order_filter_signature(a, domain, rank):
    xp = array_namespace(a, domain)
    return xp


def wiener_signature(im, mysize=None, noise=None):
    xp = array_namespace(im)
    return xp


def detrend_signature(data, axis=-1, type='linear', bp=0, overwrite_data=False):
    xp = array_namespace(data, None if isinstance(bp, int) else bp)
    return xp


def hilbert_signature(x, *args, **kwds):
    xp = array_namespace(x)
    return xp

hilbert2_signature = hilbert_signature


def lfilter_signature(b, a, x, axis=-1, zi=None):
    return array_namespace(b, a, x, zi)


def lfilter_zi_signature(b, a):
    return array_namespace(b, a)


def filtfilt_signature(b, a, x, *args, **kwds):
    return array_namespace(b, a, x)


def sosfilt_signature(sos, x, axis=-1, zi=None):
    return array_namespace(sos, x, zi)


def sosfilt_zi_signature(sos):
    return array_namespace(sos)


def sosfiltfilt_signature(sos, x, axis=-1, padtype='odd', padlen=None):
    return array_namespace(sos, x)


def deconvolve_signature(signal, divisor):
    return array_namespace(signal, divisor)


# functions we patch for dispatch
_FUNC_MAP = {
    convolve: convolve_signature,
    fftconvolve: fftconvolve_signature,
    oaconvolve: oaconvolve_signature,
    correlate: correlate_signature,
    convolve2d: convolve2d_signature,
    correlate2d: correlate2d_signature,
    medfilt: medfilt_signature,
    medfilt2d: medfilt2d_signature,
    order_filter: order_filter_signature,
    wiener: wiener_signature,
    detrend: detrend_signature,
    hilbert: hilbert_signature,
    hilbert2: hilbert2_signature,
    lfilter: lfilter_signature,
    lfilter_zi: lfilter_zi_signature, 
    deconvolve: deconvolve_signature,
    sosfilt: sosfilt_signature,
    sosfiltfilt: sosfiltfilt_signature,
    sosfilt_zi : sosfilt_zi_signature,
    filtfilt: filtfilt_signature,
}


# ### decorate ###
for func in _FUNC_MAP:
    f = (delegate_xp(_FUNC_MAP[func], MODULE_NAME)(func)
         if SCIPY_ARRAY_API
         else func)
    sys.modules[__name__].__dict__[func.__name__] = f


__all__ = [f.__name__ for f in _FUNC_MAP]
