"""
This module is deprecated -- use scipy.linalg.lapack instead
"""
from _flapack import *
import numpy as _np
@_np.deprecate(old_name="scipy.linalg.flapack", new_name="scipy.linalg.lapack")
def _deprecate(): pass
_deprecate()
