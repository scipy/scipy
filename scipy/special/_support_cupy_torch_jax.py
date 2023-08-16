import sys
import numpy as np
from scipy._lib._array_api import array_namespace
from . import _ufuncs
from ._ufuncs import *

__all__ = _ufuncs.__all__

array_api_compat_prefix = "scipy._lib.array_api_compat.array_api_compat"


def get_array_special_func(f_name, xp, n_array_args):
    # TODO: replace these `__name__` checks with `is_numpy`, etc. after
    #       gh-19005 merges.
    if xp.__name__ == f"{array_api_compat_prefix}.numpy":
        f = getattr(_ufuncs, f_name, None)
    elif xp.__name__ == f"{array_api_compat_prefix}.torch":
        f = getattr(xp.special, f_name, None)
    elif xp.__name__ == f"{array_api_compat_prefix}.cupy":
        import cupyx
        f = getattr(cupyx.scipy.special, f_name, None)
    elif xp.__name__ == f"{array_api_compat_prefix}.jax":
        f = getattr(xp.scipy.special, f_name, None)
    else:
        f_scipy = getattr(_ufuncs, f_name, None)
        def f(*args, **kwargs):
            array_args = args[:n_array_args]
            other_args = args[n_array_args:]
            array_args = [np.asarray(arg) for arg in array_args]
            out = f_scipy(*array_args, *other_args, **kwargs)
            return xp.asarray(out)

    return f


# functools.wraps doesn't work because:
# 'numpy.ufunc' object has no attribute '__module__'
def support_cupy_torch_jax(f_name, n_array_args):
    def wrapped(*args, **kwargs):
        xp = array_namespace(*args[:n_array_args])
        f = get_array_special_func(f_name, xp, n_array_args)
        return f(*args, **kwargs)
    wrapped.__doc__ = getattr(_ufuncs, f_name).__doc__
    return wrapped


array_special_func_map = {
    'log_ndtr': 1,
    'ndtr': 1,
    'ndtri': 1
}

for f_name, n_array_args in array_special_func_map.items():
    sys.modules[__name__].__dict__[f_name] = (
        support_cupy_torch_jax(f_name, n_array_args)
    )
