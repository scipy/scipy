import os
import sys
import functools

import numpy as np
import scipy
from scipy._lib._array_api import (
    array_namespace, scipy_namespace_for, is_numpy, is_torch
)
from . import _ufuncs
# These don't really need to be imported, but otherwise IDEs might not realize
# that these are defined in this file / report an error in __init__.py
from ._ufuncs import (
    log_ndtr, ndtr, ndtri, erf, erfc, i0, i0e, i1, i1e,  # noqa: F401
    gammaln, gammainc, gammaincc, logit, expit, entr, rel_entr)  # noqa: F401

_SCIPY_ARRAY_API = os.environ.get("SCIPY_ARRAY_API", False)
array_api_compat_prefix = "scipy._lib.array_api_compat"


def get_array_special_func(f_name, xp, n_array_args):
    spx = scipy_namespace_for(xp)
    f = None
    if is_numpy(xp):
        f = getattr(_ufuncs, f_name, None)
    elif is_torch(xp):
        f = getattr(xp.special, f_name, None)
    elif spx is not scipy:
        f = getattr(spx.special, f_name, None)

    if f is not None:
        return f

    # if generic array-API implementation is available, use that;
    # otherwise, fall back to NumPy/SciPy
    # Use of ` _f=_f, _xp=xp` is to avoid late-binding gotcha
    if f_name in _generic_implementations:
        _f = _generic_implementations[f_name]
        def f(*args, _f=_f, _xp=xp, **kwargs):
            return _f(*args, xp=_xp, **kwargs)
    else:
        _f = getattr(_ufuncs, f_name, None)
        def f(*args, _f=_f, _xp=xp, **kwargs):
            array_args = args[:n_array_args]
            other_args = args[n_array_args:]
            array_args = [np.asarray(arg) for arg in array_args]
            out = _f(*array_args, *other_args, **kwargs)
            return _xp.asarray(out)

    return f


def _get_shape_dtype(*args, xp):
    args = xp.broadcast_arrays(*args)
    shape = args[0].shape
    dtype = xp.result_type(*args)
    if xp.isdtype(dtype, 'integral'):
        dtype = xp.float64
        args = [xp.asarray(arg, dtype=dtype) for arg in args]
    return args, shape, dtype


def _rel_entr(x, y, *, xp):
    args, shape, dtype = _get_shape_dtype(x, y, xp=xp)
    x, y = args
    res = xp.full(x.shape, xp.inf, dtype=dtype)
    res[(x == 0) & (y >= 0)] = xp.asarray(0, dtype=dtype)
    i = (x > 0) & (y > 0)
    res[i] = x[i] * (xp.log(x[i]) - xp.log(y[i]))
    return res


_generic_implementations = {'rel_entr': _rel_entr}


# functools.wraps doesn't work because:
# 'numpy.ufunc' object has no attribute '__module__'
def support_alternative_backends(f_name, n_array_args):
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
    'entr': 1,
    'rel_entr': 2,
}

for f_name, n_array_args in array_special_func_map.items():
    f = (support_alternative_backends(f_name, n_array_args) if _SCIPY_ARRAY_API
         else getattr(_ufuncs, f_name))
    sys.modules[__name__].__dict__[f_name] = f

__all__ = list(array_special_func_map)
