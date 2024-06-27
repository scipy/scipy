import sys
import functools
from scipy._lib._array_api import (
    array_namespace, is_cupy, is_jax, scipy_namespace_for, SCIPY_ARRAY_API
)


from ._ndimage_api import *
from . import _ndimage_api
from . import _dispatchers
__all__ = _ndimage_api.__all__


MODULE_NAME = 'ndimage'


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

# ### decorate ###
for func_name in _ndimage_api.__all__:
    bare_func = getattr(_ndimage_api, func_name)
    dispatcher = getattr(_dispatchers, func_name + "_dispatcher")

    f = (dispatch_xp(dispatcher, MODULE_NAME)(bare_func)
         if SCIPY_ARRAY_API
         else bare_func)

    sys.modules[__name__].__dict__[func_name] = f

