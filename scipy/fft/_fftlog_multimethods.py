'''Multimethods for fast Hankel transforms.
'''

import numpy as np

from ._basic_uarray import _dispatch
from scipy._lib.uarray import Dispatchable


__all__ = ['fht', 'ifht']


@_dispatch
def fht(a, dln, mu, offset=0.0, bias=0.0):
    """fht multimethod."""
    return (Dispatchable(a, np.ndarray),)


@_dispatch
def ifht(A, dln, mu, offset=0.0, bias=0.0):
    """ifht multimethod."""
    return (Dispatchable(A, np.ndarray),)
