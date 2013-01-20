""" Utility functions for sparse matrix module
"""

from __future__ import division, print_function, absolute_import

__all__ = ['upcast','getdtype','isscalarlike','isintlike',
            'isshape','issequence','isdense']

import numpy as np

# keep this list syncronized with sparsetools
#supported_dtypes = ['int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32',
#        'int64', 'uint64', 'float32', 'float64',
#        'complex64', 'complex128']
supported_dtypes = ['int8','uint8','short','ushort','intc','uintc',
        'longlong','ulonglong','single','double','longdouble',
        'csingle','cdouble','clongdouble']
supported_dtypes = [np.typeDict[x] for x in supported_dtypes]

_upcast_memo = {}

def upcast(*args):
    """Returns the nearest supported sparse dtype for the
    combination of one or more types.

    upcast(t0, t1, ..., tn) -> T  where T is a supported dtype

    Examples
    --------

    >>> upcast('int32')
    <type 'numpy.int32'>
    >>> upcast('bool')
    <type 'numpy.int8'>
    >>> upcast('int32','float32')
    <type 'numpy.float64'>
    >>> upcast('bool',complex,float)
    <type 'numpy.complex128'>

    """

    t = _upcast_memo.get(hash(args))
    if t is not None:
        return t

    upcast = np.find_common_type(args, [])

    for t in supported_dtypes:
        if np.can_cast(upcast, t):
            _upcast_memo[hash(args)] = t
            return t

    raise TypeError('no supported conversion for types: %r' % (args,))

def upcast_char(*args):
    """Same as `upcast` but taking dtype.char as input (faster)."""
    t = _upcast_memo.get(args)
    if t is not None:
        return t
    t = upcast(*map(np.dtype, args))
    _upcast_memo[args] = t
    return t

def to_native(A):
    return np.asarray(A,dtype=A.dtype.newbyteorder('native'))


def getdtype(dtype, a=None, default=None):
    """Function used to simplify argument processing.  If 'dtype' is not
    specified (is None), returns a.dtype; otherwise returns a np.dtype
    object created from the specified dtype argument.  If 'dtype' and 'a'
    are both None, construct a data type out of the 'default' parameter.
    Furthermore, 'dtype' must be in 'allowed' set.
    """
    #TODO is this really what we want?
    canCast = True
    if dtype is None:
        try:
            newdtype = a.dtype
        except AttributeError:
            if default is not None:
                newdtype = np.dtype(default)
                canCast = False
            else:
                raise TypeError("could not interpret data type")
    else:
        newdtype = np.dtype(dtype)

    return newdtype

def isscalarlike(x):
    """Is x either a scalar, an array scalar, or a 0-dim array?"""
    return np.isscalar(x) or (isdense(x) and x.ndim == 0)

def isintlike(x):
    """Is x appropriate as an index into a sparse matrix? Returns True
    if it can be cast safely to a machine int.
    """
    if issequence(x):
        return False
    else:
        try:
            if int(x) == x:
                return True
            else:
                return False
        except TypeError:
            return False

def isshape(x):
    """Is x a valid 2-tuple of dimensions?
    """
    try:
        # Assume it's a tuple of matrix dimensions (M, N)
        (M, N) = x
    except:
        return False
    else:
        if isintlike(M) and isintlike(N):
            if np.rank(M) == 0 and np.rank(N) == 0:
                return True
        return False


def issequence(t):
    return isinstance(t, (list, tuple))\
           or (isinstance(t, np.ndarray) and (t.ndim == 1))

def isdense(x):
    return isinstance(x, np.ndarray)
