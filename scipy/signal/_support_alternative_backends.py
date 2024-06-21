import sys

from scipy._lib._array_api import support_alternative_backends, SCIPY_ARRAY_API
from . import _signaltools
# These don't really need to be imported, but otherwise IDEs might not realize
# that these are defined in this file / report an error in __init__.py
from ._signaltools import (
    convolve, correlate, fftconvolve,  # noqa: F401
)

_generic_implementations = {}

array_special_func_map = {
    'convolve': 2,
    'correlate': 2,
    'fftconvolve': 2,
}

for f_name, n_array_args in array_special_func_map.items():
    _f = getattr(_signaltools, f_name)
    _f_generic = _generic_implementations.get(f_name, None)
    f = (support_alternative_backends(_f, n_array_args, 'signal', _f_generic)
         if SCIPY_ARRAY_API else getattr(_signaltools, f_name))
    sys.modules[__name__].__dict__[f_name] = f

__all__ = list(array_special_func_map)
