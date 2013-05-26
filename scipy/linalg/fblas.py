"""
This module is deprecated -- use scipy.linalg.blas instead
"""
from __future__ import division, print_function, absolute_import

from ._fblas import *
import numpy as _np


@_np.deprecate(old_name="scipy.linalg.fblas", new_name="scipy.linalg.blas")
def _deprecate():
    pass
_deprecate()
