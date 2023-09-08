import os
import sys
import functools

import numpy as np
from scipy._lib._docscrape import FunctionDoc
from scipy._lib._array_api import array_namespace, is_cupy, is_torch, is_numpy
from . import _ufuncs
# These don't really need to be imported, but otherwise IDEs might not realize
# that these are defined in this file / report an error in __init__.py
from ._ufuncs import (
    log_ndtr, ndtr, ndtri, erf, erfc, i0, i0e, i1, i1e, gammaln, gammainc,
    gammaincc, logit, expit,)  # noqa

_SCIPY_ARRAY_API = os.environ.get("SCIPY_ARRAY_API", False)
array_api_compat_prefix = "scipy._lib.array_api_compat.array_api_compat"


def get_array_special_func(f_name, xp, n_array_args):
    if is_numpy(xp):
        f = getattr(_ufuncs, f_name, None)
    elif is_torch(xp):
        f = getattr(xp.special, f_name, None)
    elif is_cupy(xp):
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
    func = getattr(_ufuncs, f_name)

    @functools.wraps(func)
    def wrapped(*args, **kwargs):
        xp = array_namespace(*args[:n_array_args])
        f = get_array_special_func(f_name, xp, n_array_args)
        return f(*args, **kwargs)

    return wrapped


array_special_func_map = {
    'log_ndtr': 1,
    'ndtr': 1,
    'ndtri': 1,
    'erf': 1,
    'erfc': 1,
    'i0': 1,
    'i0e': 1,
    'i1': 1,
    'i1e': 1,
    'gammaln': 1,
    'gammainc': 2,
    'gammaincc': 2,
    'logit': 1,
    'expit': 1,
}

for f_name, n_array_args in array_special_func_map.items():
    f = (support_cupy_torch_jax(f_name, n_array_args) if _SCIPY_ARRAY_API
         else getattr(_ufuncs, f_name))

    # Add standard note about CuPy/Torch/JAX support to documentation
    doc = FunctionDoc(f)
    first_args_note = ("first positional argument is" if n_array_args ==1
                       else f"first {n_array_args} positional arguments are")
    doc["Notes"] += [""] if len(doc["Notes"]) else []
    doc["Notes"] += [
        "This function has preliminary support for CuPy, PyTorch, and JAX \n"
        "arrays. When environment variable ``SCIPY_ARRAY_API=1``, this \n"
        "function automatically passes the argument array(s) to the \n"
        "appropriate function of their native library and returns the \n"
        f"result. In this case, only the {first_args_note} supported. \n"
        "Other arguments will be passed to the underlying function, but the \n"
        "behavior is not tested."]
    # There are consistently four lines at the top of these docstrings before
    # the top-line description starts.
    doc = str(doc).split("\n", 4)[4].strip()
    f.__doc__ = str(doc)

    sys.modules[__name__].__dict__[f_name] = f

__all__ = list(array_special_func_map)
