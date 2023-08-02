'''Multimethods for fast Hankel transforms.
'''

import numpy as np

from ._basic_uarray import _dispatch
from scipy._lib.uarray import Dispatchable

# The functions in this file decorated with `@_dispatch` will only be called
# from the correspondingly named public functions with NumPy arrays or
# array-likes, not with CuPy, PyTorch and other array API standard supporting
# objects. See the docstrings in `_backend.py` for more details and examples.

__all__ = ['fht', 'ifht']


@_dispatch
def fht(a, dln, mu, offset=0.0, bias=0.0):
    """fht multimethod."""
    return (Dispatchable(a, np.ndarray),)


@_dispatch
def ifht(A, dln, mu, offset=0.0, bias=0.0):
    """ifht multimethod."""
    return (Dispatchable(A, np.ndarray),)
