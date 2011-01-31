import numpy as np
from numpy.linalg import LinAlgError

__all__ = ['LinAlgError', 'norm']


def norm(a, ord=None):
    # Differs from numpy only in non-finite handling
    return np.linalg.norm(np.asarray_chkfinite(a), ord=ord)
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
