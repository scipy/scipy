import sys
import functools
from scipy._lib._array_api import (
    array_namespace, is_cupy, is_jax, scipy_namespace_for, SCIPY_ARRAY_API
)
from ._signaltools import (convolve, fftconvolve, convolve2d, oaconvolve,
                           medfilt, medfilt2d, wiener, correlate, correlate2d,
                           detrend)


MODULE_NAME = 'signal'

# https://jax.readthedocs.io/en/latest/jax.scipy.html
JAX_SIGNAL_FUNCS = [
    'fftconvolve', 'convolve', 'convolve2d', 'correlate', 'correlate2d',
    'csd', 'detrend', 'istft', 'welch'
]

def dispatch_xp(dispatcher, module_name):
    def inner(func):
        @functools.wraps(func)
        def wrapper(*args, **kwds):
            try:
                xp = dispatcher(*args, **kwds)
            except TypeError:
                # object arrays
                import numpy
                xp = numpy

            # try delegating to a cupyx/jax namesake
            if is_cupy(xp):
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



# X_dispatcher signature must match the signature of X

def convolve_dispatcher(in1, in2, mode='full', method='auto'):
    xp = array_namespace(in1, in2)
    return xp


def fftconvolve_dispatcher(in1, in2, mode="full", axes=None):
    xp = array_namespace(in1, in2)
    return xp

def oaconvolve_dispatcher(in1, in2, mode="full", axes=None):
    xp = array_namespace(in1, in2)
    return xp


def correlate_dispatcher(in1, in2, mode='full', method='auto'):
    xp = array_namespace(in1, in2)
    return xp


def convolve2d_dispatcher(in1, in2, mode='full', boundary='fill', fillvalue=0):
    xp = array_namespace(in1, in2)
    return xp

correlate2d_dispatcher = convolve2d_dispatcher


def medfilt_dispatcher(volume, kernel_size=None):
    xp = array_namespace(volume)
    return xp


def medfilt2d_dispatcher(input, kernel_size=3):
    xp = array_namespace(input)
    return xp


def wiener_dispatcher(im, mysize=None, noise=None):
    xp = array_namespace(im)
    return xp

def detrend_dispatcher(data, axis=-1, type='linear', bp=0, overwrite_data=False):
    xp = array_namespace(data, None if isinstance(bp, int) else bp)
    return xp


# functions we patch for dispatch
_FUNC_MAP = {
    convolve: convolve_dispatcher,
    fftconvolve: fftconvolve_dispatcher,
    oaconvolve: oaconvolve_dispatcher,
    correlate: correlate_dispatcher,
    convolve2d: convolve2d_dispatcher,
    correlate2d: correlate2d_dispatcher,
    medfilt: medfilt_dispatcher,
    medfilt2d: medfilt2d_dispatcher,
    wiener: wiener_dispatcher,
    detrend: detrend_dispatcher,
}


# ### decorate ###
for func in _FUNC_MAP:
    f = (dispatch_xp(_FUNC_MAP[func], MODULE_NAME)(func)
         if SCIPY_ARRAY_API
         else func)
    sys.modules[__name__].__dict__[func.__name__] = f


__all__ = [f.__name__ for f in _FUNC_MAP]
