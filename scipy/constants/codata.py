# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.constants` namespace for importing the functions
# included below.

import warnings
from . import _codata

__all__ = [  # noqa: F822
    'physical_constants', 'value', 'unit', 'precision', 'find',
    'ConstantWarning'
]


def __dir__():
    return __all__


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(
            "scipy.constants.codata is deprecated and has no attribute "
            f"{name}. Try looking in scipy.constants instead.")

    warnings.warn(f"Please use `{name}` from the `scipy.constants` namespace, "
                  "the `scipy.constants.codata` namespace is deprecated.",
                  category=DeprecationWarning, stacklevel=2)

    return getattr(_codata, name)
