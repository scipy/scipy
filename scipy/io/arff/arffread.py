# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.io.arff` namespace for importing the functions
# included below.

import warnings
from . import _arffread

__all__ = [  # noqa: F822
    'MetaData', 'loadarff', 'ArffError', 'ParseArffError'
]


def __dir__():
    return __all__


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(
            "scipy.io.arff.arffread is deprecated and has no attribute "
            f"{name}. Try looking in scipy.io.arff instead.")

    warnings.warn(f"Please use `{name}` from the `scipy.io.arff` namespace, "
                  "the `scipy.io.arff.arffread` namespace is deprecated.",
                  category=DeprecationWarning, stacklevel=2)

    return getattr(_arffread, name)
