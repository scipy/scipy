from scipy._lib._array_api import array_namespace
from . import _ufuncs
from ._ufuncs import *

__all__ = _ufuncs.__all__

# functools.wraps doesn't work because:
# 'numpy.ufunc' object has no attribute '__module__'
def dispatcher(f_name, n_args):
    def wrapped(*args, **kwargs):
        xp = array_namespace(*args[:n_args])
        if xp.__name__.endswith('numpy'):
            f = getattr(_ufuncs, f_name)
        else:
            special = getattr(xp, "special")
            f = getattr(special, f_name)
        return f(*args, **kwargs)
    wrapped.__doc__ = getattr(_ufuncs, f_name).__doc__
    return wrapped

ndtr = dispatcher("ndtr", 1)
