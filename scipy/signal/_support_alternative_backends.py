import sys
import functools
from scipy._lib._array_api import (
    array_namespace, is_cupy, is_jax, scipy_namespace_for, SCIPY_ARRAY_API
)

import numpy as np
from ._signal_api import *   # noqa: F403
from . import _signal_api
from . import _delegators
__all__ = _signal_api.__all__


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



# ### decorate ###
for func_name in _signal_api.__all__:
    bare_func = getattr(_signal_api, func_name)
    delegator = getattr(_delegators, func_name + "_signature", None)

    f = (delegate_xp(delegator, MODULE_NAME)(bare_func)
         if SCIPY_ARRAY_API
         else bare_func)

    # add the decorated function to the namespace, to be imported in __init__.py
    vars()[func_name] = f
