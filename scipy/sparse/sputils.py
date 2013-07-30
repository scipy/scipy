""" Utility functions for sparse matrix module
"""

from __future__ import division, print_function, absolute_import

__all__ = ['upcast','getdtype','isscalarlike','isintlike',
            'isshape','issequence','isdense','ismatrix']

import numpy as np

# keep this list syncronized with sparsetools
#supported_dtypes = ['bool', 'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32',
#        'int64', 'uint64', 'float32', 'float64',
#        'complex64', 'complex128']
supported_dtypes = ['bool', 'int8','uint8','short','ushort','intc','uintc',
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
    <type 'numpy.bool_'>
    >>> upcast('int32','float32')
    <type 'numpy.float64'>
    >>> upcast('bool',complex,float)
    <type 'numpy.complex128'>

    """

    t = _upcast_memo.get(hash(args))
    if t is not None:
        return t

    if np.all([np.issubdtype(np.bool, arg) for arg in args]):
        # numpy 1.5.x compat - it gives int8 for
        # np.find_common_type([np.bool, np.bool)
        upcast = np.bool
    else:
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
    return (isinstance(t, (list, tuple)) and np.isscalar(t[0]))\
           or (isinstance(t, np.ndarray) and (t.ndim == 1))


def ismatrix(t):
    return ((issequence(t) and issequence(t[0]) and np.isscalar(t[0][0]))
            or (isinstance(t, np.ndarray) and t.ndim == 2))


def isdense(x):
    return isinstance(x, np.ndarray)


class IndexMixin(object):
    """
    This class simply exists to hold the methods necessary for fancy indexing.
    """
    def _slicetoarange(self, j, shape):
        start, stop, step = j.indices(shape)
        return np.arange(start, stop, step)

    def _unpack_index(self, index):
        """ Parse index.
        """
        if isinstance(index, tuple):
            if len(index) == 2:
                return index
            elif len(index) == 1:
                return index[0], slice(None)
            else:
                raise IndexError('invalid number of indices')
        else:
            return index, slice(None)

    def _boolean_index_to_array(self, i):
        if i.ndim < 2:
            i = i.nonzero()
        elif (i.shape[0] > self.shape[0] or i.shape[1] > self.shape[1] or
              i.ndim > 2):
            raise IndexError('invalid index shape')
        else:
            i = i.nonzero()
        return i

    def _check_boolean(self, i, j):
        if isinstance(i, np.ndarray) and i.dtype.kind == 'b':
            i = self._boolean_index_to_array(i)
            if len(i) == 2:
                if isinstance(j, slice):
                    j = i[1]
                else:
                    raise ValueError('too many indices for array')
            i = i[0]
        if isinstance(j, np.ndarray) and j.dtype.kind == 'b':
            j = self._boolean_index_to_array(j)[0]
        return i, j

    def _index_to_arrays(self, i, j):
        i, j = self._check_boolean(i, j)

        i_slice = isinstance(i, slice)
        if i_slice:
            i = self._slicetoarange(i, self.shape[0])[:,None]
        else:
            i = np.atleast_1d(i)

        if isinstance(j, slice):
            j = self._slicetoarange(j, self.shape[1])[None,:]
            if i.ndim == 1:
                i = i[:,None]
            elif not i_slice:
                raise IndexError('index returns 3-dim structure')
        elif isscalarlike(j):
            # row vector special case
            j = np.atleast_1d(j)
            if i.ndim == 1:
                i, j = np.broadcast_arrays(i, j)
                i = i[:, None]
                j = j[:, None]
                return i, j
        else:
            j = np.atleast_1d(j)
            if i_slice and j.ndim > 1:
                raise IndexError('index returns 3-dim structure')

        i, j = np.broadcast_arrays(i, j)

        if i.ndim == 1:
            # return column vectors for 1-D indexing
            i = i[None,:]
            j = j[None,:]
        elif i.ndim > 2:
            raise IndexError("Index dimension must be <= 2")

        return i, j
