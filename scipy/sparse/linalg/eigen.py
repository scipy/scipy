# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.sparse.linalg` namespace for importing the functions
# included below.

import warnings
from . import _eigen


__all__ = [  # noqa: F822
    'ArpackError', 'ArpackNoConvergence',
    'eigs', 'eigsh', 'lobpcg', 'svds', 'arpack', 'test'
]

eigen_modules = ['arpack']

def __dir__():
    return __all__


def __getattr__(name):
    if name not in __all__ and name not in eigen_modules:
        raise AttributeError(
            "scipy.sparse.linalg.eigen is deprecated and has no attribute "
            f"{name}. Try looking in scipy.sparse.linalg instead.")

    warnings.warn(f"Please use `{name}` from the `scipy.sparse.linalg` namespace, "
                  "the `scipy.sparse.linalg.eigen` namespace is deprecated.",
                  category=DeprecationWarning, stacklevel=2)

    return getattr(_eigen, name)
