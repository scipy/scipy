"""Dictionary Of Keys based matrix"""

from __future__ import division, print_function, absolute_import

__docformat__ = "restructuredtext en"

__all__ = ['dok_matrix', 'isspmatrix_dok']

import functools
import operator

import numpy as np

from scipy._lib.six import zip as izip, xrange, iteritems, itervalues

from .base import spmatrix, isspmatrix
from .sputils import (isdense, getdtype, isshape, isintlike, isscalarlike,
                      upcast, upcast_scalar, IndexMixin, get_index_dtype)

try:
    from operator import isSequenceType as _is_sequence
except ImportError:
    def _is_sequence(x):
        return (hasattr(x, '__len__') or hasattr(x, '__next__')
                or hasattr(x, 'next'))


class dok_matrix(spmatrix, IndexMixin, dict):
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
    >>> import numpy as np
    >>> from scipy.sparse import dok_matrix
    >>> S = dok_matrix((5, 5), dtype=np.float32)
    >>> for i in range(5):
    ...     for j in range(5):
    ...         S[i, j] = i + j    # Update element

    """
    format = 'dok'

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

    def getnnz(self, axis=None):
        if axis is not None:
            raise NotImplementedError("getnnz over an axis is not implemented "
                                      "for DOK format")
        return dict.__len__(self)

    def count_nonzero(self):
        return sum(x != 0 for x in itervalues(self))

    getnnz.__doc__ = spmatrix.getnnz.__doc__
    count_nonzero.__doc__ = spmatrix.count_nonzero.__doc__

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

    def __getitem__(self, index):
        """If key=(i,j) is a pair of integers, return the corresponding
        element.  If either i or j is a slice or sequence, return a new sparse
        matrix with just these elements.
        """
        zero = self.dtype.type(0)
        i, j = self._unpack_index(index)

        i_intlike = isintlike(i)
        j_intlike = isintlike(j)

        if i_intlike and j_intlike:
            # Scalar index case
            i = int(i)
            j = int(j)
            if i < 0:
                i += self.shape[0]
            if i < 0 or i >= self.shape[0]:
                raise IndexError('index out of bounds')
            if j < 0:
                j += self.shape[1]
            if j < 0 or j >= self.shape[1]:
                raise IndexError('index out of bounds')
            return dict.get(self, (i,j), zero)
        elif ((i_intlike or isinstance(i, slice)) and
              (j_intlike or isinstance(j, slice))):
            # Fast path for slicing very sparse matrices
            i_slice = slice(i, i+1) if i_intlike else i
            j_slice = slice(j, j+1) if j_intlike else j
            i_indices = i_slice.indices(self.shape[0])
            j_indices = j_slice.indices(self.shape[1])
            i_seq = xrange(*i_indices)
            j_seq = xrange(*j_indices)
            newshape = (len(i_seq), len(j_seq))
            newsize = _prod(newshape)

            if len(self) < 2*newsize and newsize != 0:
                # Switch to the fast path only when advantageous
                # (count the iterations in the loops, adjust for complexity)
                #
                # We also don't handle newsize == 0 here (if
                # i/j_intlike, it can mean index i or j was out of
                # bounds)
                return self._getitem_ranges(i_indices, j_indices, newshape)

        i, j = self._index_to_arrays(i, j)

        if i.size == 0:
            return dok_matrix(i.shape, dtype=self.dtype)

        min_i = i.min()
        if min_i < -self.shape[0] or i.max() >= self.shape[0]:
            raise IndexError('index (%d) out of range -%d to %d)' %
                             (i.min(), self.shape[0], self.shape[0]-1))
        if min_i < 0:
            i = i.copy()
            i[i < 0] += self.shape[0]

        min_j = j.min()
        if min_j < -self.shape[1] or j.max() >= self.shape[1]:
            raise IndexError('index (%d) out of range -%d to %d)' %
                             (j.min(), self.shape[1], self.shape[1]-1))
        if min_j < 0:
            j = j.copy()
            j[j < 0] += self.shape[1]

        newdok = dok_matrix(i.shape, dtype=self.dtype)

        for a in xrange(i.shape[0]):
            for b in xrange(i.shape[1]):
                v = dict.get(self, (i[a,b], j[a,b]), zero)
                if v != 0:
                    dict.__setitem__(newdok, (a, b), v)

        return newdok

    def _getitem_ranges(self, i_indices, j_indices, shape):
        # performance golf: we don't want Numpy scalars here, they are slow
        i_start, i_stop, i_stride = map(int, i_indices)
        j_start, j_stop, j_stride = map(int, j_indices)

        newdok = dok_matrix(shape, dtype=self.dtype)

        for (ii, jj) in self.keys():
            # ditto for numpy scalars
            ii = int(ii)
            jj = int(jj)
            a, ra = divmod(ii - i_start, i_stride)
            if a < 0 or a >= shape[0] or ra != 0:
                continue
            b, rb = divmod(jj - j_start, j_stride)
            if b < 0 or b >= shape[1] or rb != 0:
                continue
            dict.__setitem__(newdok, (a, b),
                             dict.__getitem__(self, (ii, jj)))

        return newdok

    def __setitem__(self, index, x):
        if isinstance(index, tuple) and len(index) == 2:
            # Integer index fast path
            i, j = index
            if (isintlike(i) and isintlike(j) and 0 <= i < self.shape[0]
                    and 0 <= j < self.shape[1]):
                v = np.asarray(x, dtype=self.dtype)
                if v.ndim == 0 and v != 0:
                    dict.__setitem__(self, (int(i), int(j)), v[()])
                    return

        i, j = self._unpack_index(index)
        i, j = self._index_to_arrays(i, j)

        if isspmatrix(x):
            x = x.toarray()

        # Make x and i into the same shape
        x = np.asarray(x, dtype=self.dtype)
        x, _ = np.broadcast_arrays(x, i)

        if x.shape != i.shape:
            raise ValueError("shape mismatch in assignment")

        if np.size(x) == 0:
            return

        min_i = i.min()
        if min_i < -self.shape[0] or i.max() >= self.shape[0]:
            raise IndexError('index (%d) out of range -%d to %d)' %
                             (i.min(), self.shape[0], self.shape[0]-1))
        if min_i < 0:
            i = i.copy()
            i[i < 0] += self.shape[0]

        min_j = j.min()
        if min_j < -self.shape[1] or j.max() >= self.shape[1]:
            raise IndexError('index (%d) out of range -%d to %d)' %
                             (j.min(), self.shape[1], self.shape[1]-1))
        if min_j < 0:
            j = j.copy()
            j[j < 0] += self.shape[1]

        dict.update(self, izip(izip(i.flat, j.flat), x.flat))

        if 0 in x:
            zeroes = x == 0
            for key in izip(i[zeroes].flat, j[zeroes].flat):
                if dict.__getitem__(self, key) == 0:
                    # may have been superseded by later update
                    del self[key]

    def __add__(self, other):
        # First check if argument is a scalar
        if isscalarlike(other):
            res_dtype = upcast_scalar(self.dtype, other)
            new = dok_matrix(self.shape, dtype=res_dtype)
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
            # We could alternatively set the dimensions to the largest of
            # the two matrices to be summed.  Would this be a good idea?
            res_dtype = upcast(self.dtype, other.dtype)
            new = dok_matrix(self.shape, dtype=res_dtype)
            new.update(self)
            with np.errstate(over='ignore'):
                for key in other.keys():
                    new[key] += other[key]
        elif isspmatrix(other):
            csc = self.tocsc()
            new = csc + other
        elif isdense(other):
            new = self.todense() + other
        else:
            return NotImplemented
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
            return NotImplemented
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
            return NotImplemented

    def __truediv__(self, other):
        if isscalarlike(other):
            res_dtype = upcast_scalar(self.dtype, other)
            new = dok_matrix(self.shape, dtype=res_dtype)
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
            return NotImplemented

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
            idx_dtype = get_index_dtype(maxval=max(self.shape[0], self.shape[1]))
            data = np.asarray(_list(self.values()), dtype=self.dtype)
            indices = np.asarray(_list(self.keys()), dtype=idx_dtype).T
            return coo_matrix((data,indices), shape=self.shape, dtype=self.dtype)

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


def _prod(x):
    """Product of a list of numbers; ~40x faster vs np.prod for Python tuples"""
    if len(x) == 0:
        return 1
    return functools.reduce(operator.mul, x)
