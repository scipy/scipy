"""Dictionary Of Keys based matrix"""

from __future__ import division, print_function, absolute_import

__docformat__ = "restructuredtext en"

__all__ = ['dok_matrix', 'isspmatrix_dok']


import numpy as np

from scipy.lib.six.moves import zip as izip, xrange
from scipy.lib.six import iteritems

from .base import spmatrix, isspmatrix
from .sputils import (isdense, getdtype, isshape, isintlike, isscalarlike,
                      upcast, upcast_scalar)

try:
    from operator import isSequenceType as _is_sequence
except ImportError:
    def _is_sequence(x):
        return (hasattr(x, '__len__') or hasattr(x, '__next__')
                or hasattr(x, 'next'))


class dok_matrix(spmatrix, dict):
    """
    Dictionary Of Keys based sparse matrix.

    This is an efficient structure for constructing sparse
    matrices incrementally.

    This can be instantiated in several ways:
        dok_matrix(D)
            with a dense matrix, D

        dok_matrix(S)
            with a sparse matrix, S

        dok_matrix((M,N), [dtype])
            create the matrix with initial shape (M,N)
            dtype is optional, defaulting to dtype='d'

    Attributes
    ----------
    dtype : dtype
        Data type of the matrix
    shape : 2-tuple
        Shape of the matrix
    ndim : int
        Number of dimensions (this is always 2)
    nnz
        Number of nonzero elements

    Notes
    -----

    Sparse matrices can be used in arithmetic operations: they support
    addition, subtraction, multiplication, division, and matrix power.

    Allows for efficient O(1) access of individual elements.
    Duplicates are not allowed.
    Can be efficiently converted to a coo_matrix once constructed.

    Examples
    --------
    >>> from scipy.sparse import *
    >>> from scipy import *
    >>> S = dok_matrix((5,5), dtype=float32)
    >>> for i in range(5):
    >>>     for j in range(5):
    >>>         S[i,j] = i+j # Update element

    """

    def __init__(self, arg1, shape=None, dtype=None, copy=False):
        dict.__init__(self)
        spmatrix.__init__(self)

        self.dtype = getdtype(dtype, default=float)
        if isinstance(arg1, tuple) and isshape(arg1):  # (M,N)
            M, N = arg1
            self.shape = (M, N)
        elif isspmatrix(arg1):  # Sparse ctor
            if isspmatrix_dok(arg1) and copy:
                arg1 = arg1.copy()
            else:
                arg1 = arg1.todok()

            if dtype is not None:
                arg1 = arg1.astype(dtype)

            self.update(arg1)
            self.shape = arg1.shape
            self.dtype = arg1.dtype
        else:  # Dense ctor
            try:
                arg1 = np.asarray(arg1)
            except:
                raise TypeError('invalid input format')

            if len(arg1.shape) != 2:
                raise TypeError('expected rank <=2 dense array or matrix')

            from .coo import coo_matrix
            d = coo_matrix(arg1, dtype=dtype).todok()
            self.update(d)
            self.shape = arg1.shape
            self.dtype = d.dtype

    def getnnz(self):
        return dict.__len__(self)
    nnz = property(fget=getnnz)

    def __len__(self):
        return dict.__len__(self)

    def get(self, key, default=0.):
        """This overrides the dict.get method, providing type checking
        but otherwise equivalent functionality.
        """
        try:
            i, j = key
            assert isintlike(i) and isintlike(j)
        except (AssertionError, TypeError, ValueError):
            raise IndexError('index must be a pair of integers')
        if (i < 0 or i >= self.shape[0] or j < 0 or j >= self.shape[1]):
            raise IndexError('index out of bounds')
        return dict.get(self, key, default)

    def __getitem__(self, key):
        """If key=(i,j) is a pair of integers, return the corresponding
        element.  If either i or j is a slice or sequence, return a new sparse
        matrix with just these elements.
        """
        try:
            i, j = key
        except (ValueError, TypeError):
            raise TypeError('index must be a pair of integers or slices')

        # Bounds checking
        if isintlike(i):
            if i < 0:
                i += self.shape[0]
            if i < 0 or i >= self.shape[0]:
                raise IndexError('index out of bounds')

        if isintlike(j):
            if j < 0:
                j += self.shape[1]
            if j < 0 or j >= self.shape[1]:
                raise IndexError('index out of bounds')

        # First deal with the case where both i and j are integers
        if isintlike(i) and isintlike(j):
            return dict.get(self, (i,j), 0.)
        else:
            # Either i or j is a slice, sequence, or invalid.  If i is a slice
            # or sequence, unfold it first and call __getitem__ recursively.

            if isinstance(i, slice):
                # Is there an easier way to do this?
                seq = xrange(i.start or 0, i.stop or self.shape[0], i.step or 1)
            elif _is_sequence(i):
                seq = i
            else:
                # Make sure i is an integer. (But allow it to be a subclass of int).
                if not isintlike(i):
                    raise TypeError('index must be a pair of integers or slices')
                seq = None
            if seq is not None:
                # i is a seq
                if isintlike(j):
                    # Create a new matrix of the correct dimensions
                    first = seq[0]
                    last = seq[-1]
                    if first < 0 or first >= self.shape[0] or last < 0 \
                                 or last >= self.shape[0]:
                        raise IndexError('index out of bounds')
                    newshape = (last-first+1, 1)
                    new = dok_matrix(newshape)
                    # ** This uses linear time in the size m of dimension 0:
                    # new[0:seq[-1]-seq[0]+1, 0] = \
                    #         [self.get((element, j), 0) for element in seq]
                    # ** Instead just add the non-zero elements.  This uses
                    # ** linear time in the number of non-zeros:
                    for (ii, jj) in self.keys():
                        if jj == j and ii >= first and ii <= last:
                            dict.__setitem__(new, (ii-first, 0),
                                             dict.__getitem__(self, (ii,jj)))
                else:
                    ###################################
                    # We should reshape the new matrix here!
                    ###################################
                    raise NotImplementedError("fancy indexing supported over"
                            " one axis only")
                return new

            # Below here, j is a sequence, but i is an integer
            if isinstance(j, slice):
                # Is there an easier way to do this?
                seq = xrange(j.start or 0, j.stop or self.shape[1], j.step or 1)
            elif _is_sequence(j):
                seq = j
            else:
                # j is not an integer
                raise TypeError("index must be a pair of integers or slices")

            # Create a new matrix of the correct dimensions
            first = seq[0]
            last = seq[-1]
            if first < 0 or first >= self.shape[1] or last < 0 \
                         or last >= self.shape[1]:
                raise IndexError("index out of bounds")
            newshape = (1, last-first+1)
            new = dok_matrix(newshape)
            # ** This uses linear time in the size n of dimension 1:
            # new[0, 0:seq[-1]-seq[0]+1] = \
            #         [self.get((i, element), 0) for element in seq]
            # ** Instead loop over the non-zero elements.  This is slower
            # ** if there are many non-zeros
            for (ii, jj) in self.keys():
                if ii == i and jj >= first and jj <= last:
                    dict.__setitem__(new, (0, jj-first),
                                     dict.__getitem__(self, (ii,jj)))
            return new

    def __setitem__(self, key, value):
        try:
            i, j = key
        except (ValueError, TypeError):
            raise TypeError("index must be a pair of integers or slices")

        # First deal with the case where both i and j are integers
        if isintlike(i) and isintlike(j):
            if i < 0:
                i += self.shape[0]
            if j < 0:
                j += self.shape[1]

            if i < 0 or i >= self.shape[0] or j < 0 or j >= self.shape[1]:
                raise IndexError("index out of bounds")

            if np.isscalar(value):
                if value == 0:
                    if (i,j) in self:
                        del self[(i,j)]
                else:
                    dict.__setitem__(self, (i,j), self.dtype.type(value))
            else:
                raise ValueError('setting an array element with a sequence')

        else:
            # Either i or j is a slice, sequence, or invalid.  If i is a slice
            # or sequence, unfold it first and call __setitem__ recursively.
            if isinstance(i, slice):
                # Is there an easier way to do this?
                seq = xrange(i.start or 0, i.stop or self.shape[0], i.step or 1)
            elif _is_sequence(i):
                seq = i
            else:
                # Make sure i is an integer. (But allow it to be a subclass of int).
                if not isintlike(i):
                    raise TypeError("index must be a pair of integers or slices")
                seq = None
            if seq is not None:
                # First see if 'value' is another dok_matrix of the appropriate
                # dimensions
                if isinstance(value, dok_matrix):
                    if value.shape[1] == 1:
                        for element in seq:
                            self[element, j] = value[element, 0]
                    else:
                        raise NotImplementedError("setting a 2-d slice of"
                                " a dok_matrix is not yet supported")
                elif np.isscalar(value):
                    for element in seq:
                        self[element, j] = value
                else:
                    # See if value is a sequence
                    try:
                        if len(seq) != len(value):
                            raise ValueError("index and value ranges must"
                                              " have the same length")
                    except TypeError:
                        # Not a sequence
                        raise TypeError("unsupported type for"
                                         " dok_matrix.__setitem__")

                    # Value is a sequence
                    for element, val in izip(seq, value):
                        self[element, j] = val   # don't use dict.__setitem__
                            # here, since we still want to be able to delete
                            # 0-valued keys, do type checking on 'val' (e.g. if
                            # it's a rank-1 dense array), etc.
            else:
                # Process j
                if isinstance(j, slice):
                    seq = xrange(j.start or 0, j.stop or self.shape[1], j.step or 1)
                elif _is_sequence(j):
                    seq = j
                else:
                    # j is not an integer
                    raise TypeError("index must be a pair of integers or slices")

                # First see if 'value' is another dok_matrix of the appropriate
                # dimensions
                if isinstance(value, dok_matrix):
                    if value.shape[0] == 1:
                        for element in seq:
                            self[i, element] = value[0, element]
                    else:
                        raise NotImplementedError("setting a 2-d slice of"
                                " a dok_matrix is not yet supported")
                elif np.isscalar(value):
                    for element in seq:
                        self[i, element] = value
                else:
                    # See if value is a sequence
                    try:
                        if len(seq) != len(value):
                            raise ValueError("index and value ranges must have"
                                              " the same length")
                    except TypeError:
                        # Not a sequence
                        raise TypeError("unsupported type for dok_matrix.__setitem__")
                    else:
                        for element, val in izip(seq, value):
                            self[i, element] = val

    def __add__(self, other):
        # First check if argument is a scalar
        if isscalarlike(other):
            new = dok_matrix(self.shape, dtype=self.dtype)
            # Add this scalar to every element.
            M, N = self.shape
            for i in xrange(M):
                for j in xrange(N):
                    aij = self.get((i, j), 0) + other
                    if aij != 0:
                        new[i, j] = aij
            # new.dtype.char = self.dtype.char
        elif isinstance(other, dok_matrix):
            if other.shape != self.shape:
                raise ValueError("matrix dimensions are not equal")
            # We could alternatively set the dimensions to the the largest of
            # the two matrices to be summed.  Would this be a good idea?
            new = dok_matrix(self.shape, dtype=self.dtype)
            new.update(self)
            for key in other.keys():
                new[key] += other[key]
        elif isspmatrix(other):
            csc = self.tocsc()
            new = csc + other
        elif isdense(other):
            new = self.todense() + other
        else:
            raise TypeError("data type not understood")
        return new

    def __radd__(self, other):
        # First check if argument is a scalar
        if isscalarlike(other):
            new = dok_matrix(self.shape, dtype=self.dtype)
            # Add this scalar to every element.
            M, N = self.shape
            for i in xrange(M):
                for j in xrange(N):
                    aij = self.get((i, j), 0) + other
                    if aij != 0:
                        new[i, j] = aij
        elif isinstance(other, dok_matrix):
            if other.shape != self.shape:
                raise ValueError("matrix dimensions are not equal")
            new = dok_matrix(self.shape, dtype=self.dtype)
            new.update(self)
            for key in other:
                new[key] += other[key]
        elif isspmatrix(other):
            csc = self.tocsc()
            new = csc + other
        elif isdense(other):
            new = other + self.todense()
        else:
            raise TypeError("data type not understood")
        return new

    def __neg__(self):
        new = dok_matrix(self.shape, dtype=self.dtype)
        for key in self.keys():
            new[key] = -self[key]
        return new

    def _mul_scalar(self, other):
        res_dtype = upcast_scalar(self.dtype, other)
        # Multiply this scalar by every element.
        new = dok_matrix(self.shape, dtype=res_dtype)
        for (key, val) in iteritems(self):
            new[key] = val * other
        return new

    def _mul_vector(self, other):
        # matrix * vector
        result = np.zeros(self.shape[0], dtype=upcast(self.dtype,other.dtype))
        for (i,j),v in iteritems(self):
            result[i] += v * other[j]
        return result

    def _mul_multivector(self, other):
        # matrix * multivector
        M,N = self.shape
        n_vecs = other.shape[1]  # number of column vectors
        result = np.zeros((M,n_vecs), dtype=upcast(self.dtype,other.dtype))
        for (i,j),v in iteritems(self):
            result[i,:] += v * other[j,:]
        return result

    def __imul__(self, other):
        if isscalarlike(other):
            # Multiply this scalar by every element.
            for (key, val) in iteritems(self):
                self[key] = val * other
            # new.dtype.char = self.dtype.char
            return self
        else:
            return NotImplementedError

    def __truediv__(self, other):
        if isscalarlike(other):
            new = dok_matrix(self.shape, dtype=self.dtype)
            # Multiply this scalar by every element.
            for (key, val) in iteritems(self):
                new[key] = val / other
            # new.dtype.char = self.dtype.char
            return new
        else:
            return self.tocsr() / other

    def __itruediv__(self, other):
        if isscalarlike(other):
            # Multiply this scalar by every element.
            for (key, val) in iteritems(self):
                self[key] = val / other
            return self
        else:
            return NotImplementedError

    # What should len(sparse) return? For consistency with dense matrices,
    # perhaps it should be the number of rows?  For now it returns the number
    # of non-zeros.

    def transpose(self):
        """ Return the transpose
        """
        M, N = self.shape
        new = dok_matrix((N, M), dtype=self.dtype)
        for key, value in iteritems(self):
            new[key[1], key[0]] = value
        return new

    def conjtransp(self):
        """ Return the conjugate transpose
        """
        M, N = self.shape
        new = dok_matrix((N, M), dtype=self.dtype)
        for key, value in iteritems(self):
            new[key[1], key[0]] = np.conj(value)
        return new

    def copy(self):
        new = dok_matrix(self.shape, dtype=self.dtype)
        new.update(self)
        return new

    def getrow(self, i):
        """Returns a copy of row i of the matrix as a (1 x n)
        DOK matrix.
        """
        out = self.__class__((1, self.shape[1]), dtype=self.dtype)
        for j in range(self.shape[1]):
            out[0, j] = self[i, j]
        return out

    def getcol(self, j):
        """Returns a copy of column j of the matrix as a (m x 1)
        DOK matrix.
        """
        out = self.__class__((self.shape[0], 1), dtype=self.dtype)
        for i in range(self.shape[0]):
            out[i, 0] = self[i, j]
        return out

    def tocoo(self):
        """ Return a copy of this matrix in COOrdinate format"""
        from .coo import coo_matrix
        if self.nnz == 0:
            return coo_matrix(self.shape, dtype=self.dtype)
        else:
            data = np.asarray(_list(self.values()), dtype=self.dtype)
            indices = np.asarray(_list(self.keys()), dtype=np.intc).T
            return coo_matrix((data,indices), shape=self.shape,
                              dtype=self.dtype)

    def todok(self,copy=False):
        if copy:
            return self.copy()
        else:
            return self

    def tocsr(self):
        """ Return a copy of this matrix in Compressed Sparse Row format"""
        return self.tocoo().tocsr()

    def tocsc(self):
        """ Return a copy of this matrix in Compressed Sparse Column format"""
        return self.tocoo().tocsc()

    def toarray(self, order=None, out=None):
        """See the docstring for `spmatrix.toarray`."""
        return self.tocoo().toarray(order=order, out=out)

    def resize(self, shape):
        """ Resize the matrix in-place to dimensions given by 'shape'.

        Any non-zero elements that lie outside the new shape are removed.
        """
        if not isshape(shape):
            raise TypeError("dimensions must be a 2-tuple of positive"
                             " integers")
        newM, newN = shape
        M, N = self.shape
        if newM < M or newN < N:
            # Remove all elements outside new dimensions
            for (i, j) in list(self.keys()):
                if i >= newM or j >= newN:
                    del self[i, j]
        self._shape = shape


def _list(x):
    """Force x to a list."""
    if not isinstance(x, list):
        x = list(x)
    return x


def isspmatrix_dok(x):
    return isinstance(x, dok_matrix)
