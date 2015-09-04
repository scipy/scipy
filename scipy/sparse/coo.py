""" A sparse matrix in COOrdinate or 'triplet' format"""
from __future__ import division, print_function, absolute_import

__docformat__ = "restructuredtext en"

__all__ = ['coo_matrix', 'isspmatrix_coo']

from warnings import warn

import numpy as np

from scipy._lib.six import xrange, zip as izip

from ._sparsetools import coo_tocsr, coo_todense, coo_matvec
from .base import isspmatrix
from .data import _data_matrix, _minmax_mixin
from .sputils import (upcast, upcast_char, to_native, isshape, getdtype,
        isintlike, get_index_dtype, downcast_intp_index)


class coo_matrix(_data_matrix, _minmax_mixin):
    """
    A sparse matrix in COOrdinate format.

    Also known as the 'ijv' or 'triplet' format.

    This can be instantiated in several ways:
        coo_matrix(D)
            with a dense matrix D

        coo_matrix(S)
            with another sparse matrix S (equivalent to S.tocoo())

        coo_matrix((M, N), [dtype])
            to construct an empty matrix with shape (M, N)
            dtype is optional, defaulting to dtype='d'.

        coo_matrix((data, (i, j)), [shape=(M, N)])
            to construct from three arrays:
                1. data[:]   the entries of the matrix, in any order
                2. i[:]      the row indices of the matrix entries
                3. j[:]      the column indices of the matrix entries

            Where ``A[i[k], j[k]] = data[k]``.  When shape is not
            specified, it is inferred from the index arrays

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
    data
        COO format data array of the matrix
    row
        COO format row index array of the matrix
    col
        COO format column index array of the matrix

    Notes
    -----

    Sparse matrices can be used in arithmetic operations: they support
    addition, subtraction, multiplication, division, and matrix power.

    Advantages of the COO format
        - facilitates fast conversion among sparse formats
        - permits duplicate entries (see example)
        - very fast conversion to and from CSR/CSC formats

    Disadvantages of the COO format
        - does not directly support:
            + arithmetic operations
            + slicing

    Intended Usage
        - COO is a fast format for constructing sparse matrices
        - Once a matrix has been constructed, convert to CSR or
          CSC format for fast arithmetic and matrix vector operations
        - By default when converting to CSR or CSC format, duplicate (i,j)
          entries will be summed together.  This facilitates efficient
          construction of finite element matrices and the like. (see example)

    Examples
    --------
    >>> from scipy.sparse import coo_matrix
    >>> coo_matrix((3, 4), dtype=np.int8).toarray()
    array([[0, 0, 0, 0],
           [0, 0, 0, 0],
           [0, 0, 0, 0]], dtype=int8)

    >>> row  = np.array([0, 3, 1, 0])
    >>> col  = np.array([0, 3, 1, 2])
    >>> data = np.array([4, 5, 7, 9])
    >>> coo_matrix((data, (row, col)), shape=(4, 4)).toarray()
    array([[4, 0, 9, 0],
           [0, 7, 0, 0],
           [0, 0, 0, 0],
           [0, 0, 0, 5]])

    >>> # example with duplicates
    >>> row  = np.array([0, 0, 1, 3, 1, 0, 0])
    >>> col  = np.array([0, 2, 1, 3, 1, 0, 0])
    >>> data = np.array([1, 1, 1, 1, 1, 1, 1])
    >>> coo_matrix((data, (row, col)), shape=(4, 4)).toarray()
    array([[3, 0, 1, 0],
           [0, 2, 0, 0],
           [0, 0, 0, 0],
           [0, 0, 0, 1]])

    """
    def __init__(self, arg1, shape=None, dtype=None, copy=False):
        _data_matrix.__init__(self)

        if isinstance(arg1, tuple):
            if isshape(arg1):
                M, N = arg1
                self.shape = (M,N)
                idx_dtype = get_index_dtype(maxval=max(M, N))
                self.row = np.array([], dtype=idx_dtype)
                self.col = np.array([], dtype=idx_dtype)
                self.data = np.array([], getdtype(dtype, default=float))
                self.has_canonical_format = True
            else:
                try:
                    obj, ij = arg1
                except:
                    raise TypeError('invalid input format')

                try:
                    if len(ij) != 2:
                        raise TypeError
                except TypeError:
                    raise TypeError('invalid input format')

                row, col = ij
                if shape is None:
                    if len(row) == 0 or len(col) == 0:
                        raise ValueError('cannot infer dimensions from zero '
                                         'sized index arrays')
                    M = np.max(row) + 1
                    N = np.max(col) + 1
                    self.shape = (M, N)
                else:
                    # Use 2 steps to ensure shape has length 2.
                    M, N = shape
                    self.shape = (M, N)

                idx_dtype = get_index_dtype(maxval=max(self.shape))
                self.row = np.array(row, copy=copy, dtype=idx_dtype)
                self.col = np.array(col, copy=copy, dtype=idx_dtype)
                self.data = np.array(obj, copy=copy)
                self.has_canonical_format = False

        else:
            if isspmatrix(arg1):
                if isspmatrix_coo(arg1) and copy:
                    self.row = arg1.row.copy()
                    self.col = arg1.col.copy()
                    self.data = arg1.data.copy()
                    self.shape = arg1.shape
                else:
                    coo = arg1.tocoo()
                    self.row = coo.row
                    self.col = coo.col
                    self.data = coo.data
                    self.shape = coo.shape
                self.has_canonical_format = False
            else:
                #dense argument
                M = np.atleast_2d(np.asarray(arg1))

                if M.ndim != 2:
                    raise TypeError('expected dimension <= 2 array or matrix')
                else:
                    self.shape = M.shape

                self.row, self.col = M.nonzero()
                self.data = M[self.row, self.col]
                self.has_canonical_format = True

        if dtype is not None:
            self.data = self.data.astype(dtype)

        self._check()

    def getnnz(self, axis=None):
        """Get the count of explicitly-stored values (nonzeros)

        Parameters
        ----------
        axis : None, 0, or 1
            Select between the number of values across the whole matrix, in
            each column, or in each row.
        """
        if axis is None:
            nnz = len(self.data)
            if nnz != len(self.row) or nnz != len(self.col):
                raise ValueError('row, column, and data array must all be the '
                                 'same length')

            if self.data.ndim != 1 or self.row.ndim != 1 or \
                    self.col.ndim != 1:
                raise ValueError('row, column, and data arrays must be 1-D')

            return int(nnz)

        if axis < 0:
            axis += 2
        if axis == 0:
            return np.bincount(downcast_intp_index(self.col),
                               minlength=self.shape[1])
        elif axis == 1:
            return np.bincount(downcast_intp_index(self.row),
                               minlength=self.shape[0])
        else:
            raise ValueError('axis out of bounds')
    nnz = property(fget=getnnz)

    def _check(self):
        """ Checks data structure for consistency """
        nnz = self.nnz

        # index arrays should have integer data types
        if self.row.dtype.kind != 'i':
            warn("row index array has non-integer dtype (%s)  "
                    % self.row.dtype.name)
        if self.col.dtype.kind != 'i':
            warn("col index array has non-integer dtype (%s) "
                    % self.col.dtype.name)

        idx_dtype = get_index_dtype(maxval=max(self.shape))
        self.row = np.asarray(self.row, dtype=idx_dtype)
        self.col = np.asarray(self.col, dtype=idx_dtype)
        self.data = to_native(self.data)

        if nnz > 0:
            if self.row.max() >= self.shape[0]:
                raise ValueError('row index exceeds matrix dimensions')
            if self.col.max() >= self.shape[1]:
                raise ValueError('column index exceeds matrix dimensions')
            if self.row.min() < 0:
                raise ValueError('negative row index found')
            if self.col.min() < 0:
                raise ValueError('negative column index found')

    def transpose(self, copy=False):
        M,N = self.shape
        return coo_matrix((self.data, (self.col, self.row)), shape=(N,M), copy=copy)

    def toarray(self, order=None, out=None):
        """See the docstring for `spmatrix.toarray`."""
        B = self._process_toarray_args(order, out)
        fortran = int(B.flags.f_contiguous)
        if not fortran and not B.flags.c_contiguous:
            raise ValueError("Output array must be C or F contiguous")
        M,N = self.shape
        coo_todense(M, N, self.nnz, self.row, self.col, self.data,
                    B.ravel('A'), fortran)
        return B

    def tocsc(self):
        """Return a copy of this matrix in Compressed Sparse Column format

        Duplicate entries will be summed together.

        Examples
        --------
        >>> from numpy import array
        >>> from scipy.sparse import coo_matrix
        >>> row  = array([0, 0, 1, 3, 1, 0, 0])
        >>> col  = array([0, 2, 1, 3, 1, 0, 0])
        >>> data = array([1, 1, 1, 1, 1, 1, 1])
        >>> A = coo_matrix((data, (row, col)), shape=(4, 4)).tocsc()
        >>> A.toarray()
        array([[3, 0, 1, 0],
               [0, 2, 0, 0],
               [0, 0, 0, 0],
               [0, 0, 0, 1]])

        """
        from .csc import csc_matrix
        if self.nnz == 0:
            return csc_matrix(self.shape, dtype=self.dtype)
        else:
            M,N = self.shape
            idx_dtype = get_index_dtype((self.col, self.row),
                                        maxval=max(self.nnz, M))
            indptr = np.empty(N + 1, dtype=idx_dtype)
            indices = np.empty(self.nnz, dtype=idx_dtype)
            data = np.empty(self.nnz, dtype=upcast(self.dtype))

            coo_tocsr(N, M, self.nnz,
                      self.col.astype(idx_dtype),
                      self.row.astype(idx_dtype),
                      self.data,
                      indptr, indices, data)

            A = csc_matrix((data, indices, indptr), shape=self.shape)
            A.sum_duplicates()

            return A

    def tocsr(self):
        """Return a copy of this matrix in Compressed Sparse Row format

        Duplicate entries will be summed together.

        Examples
        --------
        >>> from numpy import array
        >>> from scipy.sparse import coo_matrix
        >>> row  = array([0, 0, 1, 3, 1, 0, 0])
        >>> col  = array([0, 2, 1, 3, 1, 0, 0])
        >>> data = array([1, 1, 1, 1, 1, 1, 1])
        >>> A = coo_matrix((data, (row, col)), shape=(4, 4)).tocsr()
        >>> A.toarray()
        array([[3, 0, 1, 0],
               [0, 2, 0, 0],
               [0, 0, 0, 0],
               [0, 0, 0, 1]])

        """
        from .csr import csr_matrix
        if self.nnz == 0:
            return csr_matrix(self.shape, dtype=self.dtype)
        else:
            M,N = self.shape
            idx_dtype = get_index_dtype((self.row, self.col),
                                        maxval=max(self.nnz, N))
            indptr = np.empty(M + 1, dtype=idx_dtype)
            indices = np.empty(self.nnz, dtype=idx_dtype)
            data = np.empty(self.nnz, dtype=upcast(self.dtype))

            coo_tocsr(M, N, self.nnz,
                      self.row.astype(idx_dtype),
                      self.col.astype(idx_dtype),
                      self.data,
                      indptr,
                      indices,
                      data)

            A = csr_matrix((data, indices, indptr), shape=self.shape)
            A.sum_duplicates()

            return A

    def tocoo(self, copy=False):
        if copy:
            return self.copy()
        else:
            return self

    def todia(self):
        from .dia import dia_matrix

        ks = self.col - self.row  # the diagonal for each nonzero
        diags = np.unique(ks)

        if len(diags) > 100:
            #probably undesired, should we do something?
            #should todia() have a maxdiags parameter?
            pass

        #initialize and fill in data array
        if self.data.size == 0:
            data = np.zeros((0, 0), dtype=self.dtype)
        else:
            data = np.zeros((len(diags), self.col.max()+1), dtype=self.dtype)
            data[np.searchsorted(diags,ks), self.col] = self.data

        return dia_matrix((data,diags), shape=self.shape)

    def todok(self):
        from .dok import dok_matrix

        self.sum_duplicates()
        dok = dok_matrix((self.shape), dtype=self.dtype)
        dok.update(izip(izip(self.row,self.col),self.data))

        return dok

    def diagonal(self):
        # Could be rewritten without the python loop.
        # Data entries at the same (row, col) are summed.
        n = min(self.shape)
        ndata = self.data.shape[0]
        d = np.zeros(n, dtype=self.dtype)
        for i in xrange(ndata):
            r = self.row[i]
            if r == self.col[i]:
                d[r] += self.data[i]
        return d
    diagonal.__doc__ = _data_matrix.diagonal.__doc__

    def _setdiag(self, values, k):
        M, N = self.shape
        if values.ndim and not len(values):
            return
        idx_dtype = self.row.dtype

        # Determine which triples to keep and where to put the new ones.
        full_keep = self.col - self.row != k
        if k < 0:
            max_index = min(M+k, N)
            if values.ndim:
                max_index = min(max_index, len(values))
            keep = np.logical_or(full_keep, self.col >= max_index)
            new_row = np.arange(-k, -k + max_index, dtype=idx_dtype)
            new_col = np.arange(max_index, dtype=idx_dtype)
        else:
            max_index = min(M, N-k)
            if values.ndim:
                max_index = min(max_index, len(values))
            keep = np.logical_or(full_keep, self.row >= max_index)
            new_row = np.arange(max_index, dtype=idx_dtype)
            new_col = np.arange(k, k + max_index, dtype=idx_dtype)

        # Define the array of data consisting of the entries to be added.
        if values.ndim:
            new_data = values[:max_index]
        else:
            new_data = np.empty(max_index, dtype=self.dtype)
            new_data[:] = values

        # Update the internal structure.
        self.row = np.concatenate((self.row[keep], new_row))
        self.col = np.concatenate((self.col[keep], new_col))
        self.data = np.concatenate((self.data[keep], new_data))
        self.has_canonical_format = False

    # needed by _data_matrix
    def _with_data(self,data,copy=True):
        """Returns a matrix with the same sparsity structure as self,
        but with different data.  By default the index arrays
        (i.e. .row and .col) are copied.
        """
        if copy:
            return coo_matrix((data, (self.row.copy(), self.col.copy())),
                                   shape=self.shape, dtype=data.dtype)
        else:
            return coo_matrix((data, (self.row, self.col)),
                                   shape=self.shape, dtype=data.dtype)

    def sum_duplicates(self):
        """Eliminate duplicate matrix entries by adding them together

        This is an *in place* operation
        """
        if self.has_canonical_format or len(self.data) == 0:
            return
        order = np.lexsort((self.row,self.col))
        self.row = self.row[order]
        self.col = self.col[order]
        self.data = self.data[order]
        unique_mask = ((self.row[1:] != self.row[:-1]) |
                       (self.col[1:] != self.col[:-1]))
        unique_mask = np.append(True, unique_mask)
        self.row = self.row[unique_mask]
        self.col = self.col[unique_mask]
        unique_inds, = np.nonzero(unique_mask)
        self.data = np.add.reduceat(self.data, unique_inds, dtype=self.dtype)
        self.has_canonical_format = True

    ###########################
    # Multiplication handlers #
    ###########################

    def _mul_vector(self, other):
        #output array
        result = np.zeros(self.shape[0], dtype=upcast_char(self.dtype.char,
                                                            other.dtype.char))
        coo_matvec(self.nnz, self.row, self.col, self.data, other, result)
        return result

    def _mul_multivector(self, other):
        return np.hstack([self._mul_vector(col).reshape(-1,1) for col in other.T])


def isspmatrix_coo(x):
    return isinstance(x, coo_matrix)
