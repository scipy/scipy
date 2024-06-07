import sys
import functools
from scipy._lib._array_api import (
    array_namespace, is_cupy, is_jax, scipy_namespace_for, SCIPY_ARRAY_API
)
from ._signaltools import convolve, fftconvolve, convolve2d


MODULE_NAME = 'signal'


def dispatch_xp(dispatcher, module_name):
    def inner(func):
        @functools.wraps(func)
        def wrapper(*args, **kwds):
            xp = dispatcher(*args, **kwds)

            # try delegating to a cupyx/jax namesake
            if is_cupy(xp):
                # https://github.com/cupy/cupy/issues/8336
                import importlib
                cupyx_module = importlib.import_module(f"cupyx.scipy.{module_name}")
                cupyx_func = getattr(cupyx_module, func.__name__)
                return cupyx_func(*args, **kwds)
            elif is_jax(xp):
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


def convolve2d_dispatcher(in1, in2, mode='full', boundary='fill', fillvalue=0):
    xp = array_namespace(in1, in2)
    return xp


# functions we patch for dispatch
_FUNC_MAP = {
    convolve: convolve_dispatcher,
    fftconvolve: fftconvolve_dispatcher,
    convolve2d: convolve2d_dispatcher,
}


# ### decorate ###
for func in _FUNC_MAP:
    f = (dispatch_xp(_FUNC_MAP[func], MODULE_NAME)(func)
         if SCIPY_ARRAY_API
         else func)
    sys.modules[__name__].__dict__[func.__name__] = f


__all__ = [f.__name__ for f in _FUNC_MAP]
