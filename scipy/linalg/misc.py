import numpy as np
from numpy.linalg import LinAlgError
from blas import get_blas_funcs

__all__ = ['LinAlgError', 'norm']


def norm(a, ord=None):
    # Differs from numpy in non-finite handling and in implementation
    # of the 2-norm, where more stable BLAS routines are used.
    a = np.asarray_chkfinite(a)
    if a.ndim == 1 and ord in (2, None) and a.dtype.kind in 'fdFD':
        nrm2, = get_blas_funcs(('nrm2',), (a,))
        return nrm2(a)
    else:
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
