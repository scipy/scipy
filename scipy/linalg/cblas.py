"""
This module is deprecated -- use scipy.linalg.blas instead
"""
try:
    from _cblas import *
except ImportError:
    empty_module = True
import numpy as _np
@_np.deprecate(old_name="scipy.linalg.cblas", new_name="scipy.linalg.blas")
def _deprecate(): pass
_deprecate()
