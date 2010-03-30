import numpy as np
from numpy.linalg import LinAlgError

### Norm

def norm(a, ord=None):
    # Differs from numpy only in non-finite handling
    return np.linalg.norm(np.asarray_chkfinite(a), ord=ord)
norm.__doc__ = np.linalg.norm.__doc__

