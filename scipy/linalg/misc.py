import numpy as np
from numpy.linalg import LinAlgError

__all__ = ['LinAlgError', 'norm']


def norm(a, ord=None):
    # Differs from numpy only in non-finite handling
    return np.linalg.norm(np.asarray_chkfinite(a), ord=ord)
norm.__doc__ = np.linalg.norm.__doc__


def _datanotshared(a1,a):
    if a1 is a:
        return False
    else:
        #try comparing data pointers
        try:
            return a1.__array_interface__['data'][0] != a.__array_interface__['data'][0]
        except:
            return True