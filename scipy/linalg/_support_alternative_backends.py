import sys

from scipy._lib._array_api import support_alternative_backends, SCIPY_ARRAY_API
from . import _basic
# These don't really need to be imported, but otherwise IDEs might not realize
# that these are defined in this file / report an error in __init__.py
from ._basic import (
    solve, inv, det, pinv  # noqa: F401
)

_generic_implementations = {}

array_special_func_map = {
    'solve': 2,
    'inv': 1,
    'det': 1,
    'pinv': 1,
}

for f_name, n_array_args in array_special_func_map.items():
    _f = getattr(_basic, f_name)
    f = (support_alternative_backends(_f, n_array_args,'linalg',
                                      _generic_implementations) if SCIPY_ARRAY_API
         else _f)
    sys.modules[__name__].__dict__[f_name] = f

__all__ = list(array_special_func_map)
