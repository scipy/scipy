import sys
from scipy._lib.array_api_compat.array_api_compat import array_namespace
from . import _ufuncs
from ._ufuncs import *

__all__ = _ufuncs.__all__

array_api_comat_prefix = "scipy._lib.array_api_compat.array_api_compat"


def get_array_special_func(f_name, xp):
    if xp.__name__ == f"{array_api_comat_prefix}.numpy":
        f = getattr(_ufuncs, f_name, None)
    elif xp.__name__ == f"{array_api_comat_prefix}.torch":
        f = getattr(xp.special, f_name, None)
    elif xp.__name__ == f"{array_api_comat_prefix}.cupy":
        import cupyx
        f = getattr(cupyx.scipy.special, f_name, None)
    elif xp.__name__ == f"{array_api_comat_prefix}.jax":
        f = getattr(xp.scipy.special, f_name, None)
    else:
        f = None

    if f is None:
        raise ValueError("SciPy cannot dispatch to special function "
                         f"'{f_name}' in array library '{xp.__name__}'.")
    return f


# functools.wraps doesn't work because:
# 'numpy.ufunc' object has no attribute '__module__'
def support_cupy_torch_jax(f_name, n_args):
    def wrapped(*args, **kwargs):
        xp = array_namespace(*args[:n_args])
        f = get_array_special_func(f_name, xp)
        return f(*args, **kwargs)
    wrapped.__doc__ = getattr(_ufuncs, f_name).__doc__
    return wrapped


array_special_func_map = {
    'log_ndtr': 1,
    'ndtr': 1,
    'ndtri': 1
}

for f_name, n_args in array_special_func_map.items():
    sys.modules[__name__].__dict__[f_name] = (
        support_cupy_torch_jax(f_name, n_args)
    )
