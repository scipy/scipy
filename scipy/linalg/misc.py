import numpy as np
from numpy.linalg import LinAlgError
import fblas

__all__ = ['LinAlgError', 'norm']

_nrm2_prefix = {'f' : 's', 'F': 'sc', 'D': 'dz'}

def norm(a, ord=None):
    # Differs from numpy only in non-finite handling and the use of
    # blas
    a = np.asarray_chkfinite(a)
    if ord in (None, 2) and (a.ndim == 1) and (a.dtype.char in 'fdFD'):
        # use blas for fast and stable euclidean norm
        func_name = _nrm2_prefix.get(a.dtype.char, 'd') + 'nrm2'
        nrm2 = getattr(fblas, func_name)
        return nrm2(a)
    return np.linalg.norm(a, ord=ord)

norm.__doc__ = np.linalg.norm.__doc__

def _datacopied(arr, original):
    """
    Strict check for `arr` not sharing any data with `original`,
    under the assumption that arr = asarray(original)

    """
    if arr is original:
        return False
    if not isinstance(original, np.ndarray) and hasattr(original, '__array__'):
        return False
    return arr.base is None
