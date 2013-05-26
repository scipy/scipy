"""
This module is deprecated -- use scipy.linalg.lapack instead
"""
from __future__ import division, print_function, absolute_import

try:
    from ._clapack import *
except ImportError:
    empty_module = True
import numpy as _np


@_np.deprecate(old_name="scipy.linalg.clapack", new_name="scipy.linalg.lapack")
def _deprecate():
    pass
_deprecate()
