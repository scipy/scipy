"""
Module to read / write MATLAB .mat files.
"""

import warnings
from ._matlab import (loadmat, savemat, whosmat, byteordercodes,
                      matfile_version, MatReadError, MatReadWarning,
                      MatWriteError, MatlabOpaque, MatlabFunction, mat_struct)


__all__ = [  # noqa: F822
    'loadmat', 'savemat', 'whosmat', 'byteordercodes', 'matfile_version',
    'MatReadError', 'MatReadWarning', 'MatWriteError', 'MatlabOpaque',
    'MatlabFunction', 'mat_struct',
]


def __dir__():
    return __all__


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(
            "scipy.io.matlab is deprecated and has no attribute "
            f"{name}. Try looking in scipy.io instead.")

    warnings.warn(f"Please use `{name}` from the `scipy.io` namespace, "
                  "the `scipy.io.matlab` namespace is deprecated.",
                  category=DeprecationWarning, stacklevel=2)

    return getattr(_matlab, name)
