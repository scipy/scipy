# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.stats` namespace for importing the functions
# included below.

import warnings
from . import _statlib  # type: ignore


__all__ = [  # noqa: F822
    'swilk',
    'gscale'
]


def __dir__():
    return __all__


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(
            "scipy.stats.statlib is deprecated and has no attribute "
            f"{name}. Try looking in scipy.stats instead.")

    warnings.warn(f"Please use `{name}` from the `scipy.stats` namespace, "
                  "the `scipy.stats.statlib` namespace is deprecated.",
                  category=DeprecationWarning, stacklevel=2)

    return getattr(_statlib, name)
