# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.io` namespace for importing the functions
# included below.

import warnings
from . import _netcdf

__all__ = [  # noqa: F822
    'netcdf_file', 'netcdf_variable'
]


def __dir__():
    return __all__


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(
            "scipy.io.netcdf is deprecated and has no attribute "
            f"{name}. Try looking in scipy.io instead.")

    warnings.warn(f"Please use `{name}` from the `scipy.io` namespace, "
                  "the `scipy.io.netcdf` namespace is deprecated.",
                  category=DeprecationWarning, stacklevel=2)

    return getattr(_netcdf, name)
