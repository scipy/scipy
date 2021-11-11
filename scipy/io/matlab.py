"""
MATLAB® file utilies (:mod:`scipy.io.matlab`)
=============================================

.. currentmodule:: scipy.io.matlab

This submodule is meant to provide lower-level file utilies related to reading
and writing MATLAB files.

.. autosummary::
   :toctree: generated/

   matfile_version - Get the MATLAB file version
   MatReadError - Exception indicating a read issue
   MatReadWarning - Warning class for read issues
   MatWriteError - Exception indicating a write issue
   mat_struct - Class used when ``struct_as_record=False``

.. autosummary::
   :toctree: generated/
   :template: autosummary/ndarray_subclass.rst
   :nosignatures:

   MatlabObject - Class for a MATLAB object
   MatlabOpaque - Class for a MATLAB opaque matrix
   MatlabFunction - Class for a MATLAB function object

The following utilities that live in the :mod:`scipy.io`
namespace also exist in this namespace:

.. autosummary::
   :toctree: generated/

   loadmat - Read a MATLAB style mat file (version 4 through 7.1)
   savemat - Write a MATLAB style mat file (version 4 through 7.1)
   whosmat - List contents of a MATLAB style mat file (version 4 through 7.1)
"""

import warnings
from ._matlab import (loadmat, savemat, whosmat, byteordercodes,
                      matfile_version, MatReadError, MatReadWarning,
                      MatWriteError, MatlabOpaque, MatlabFunction, mat_struct,
                      MatlabObject)
from . import _matlab


__all__ = [  # noqa: F822
    'loadmat', 'savemat', 'whosmat', 'byteordercodes', 'matfile_version',
    'MatReadError', 'MatReadWarning', 'MatWriteError', 'MatlabOpaque',
    'MatlabFunction', 'mat_struct', 'MatlabObject',
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
