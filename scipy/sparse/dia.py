"""Sparse DIAgonal format"""

from __future__ import division, print_function, absolute_import

__docformat__ = "restructuredtext en"

__all__ = ['dia_matrix', 'isspmatrix_dia']

import numpy as np

from .base import isspmatrix, _formats
from .data import _data_matrix
from .sputils import isshape, upcast, upcast_char, getdtype, get_index_dtype
from ._sparsetools import dia_matvec


class dia_matrix(_data_matrix):
    """Sparse matrix with DIAgonal storage

    This can be instantiated in several ways:
        dia_matrix(D)
            with a dense matrix

        dia_matrix(S)
            with another sparse matrix S (equivalent to S.todia())

        dia_matrix((M, N), [dtype])
            to construct an empty matrix with shape (M, N),
            dtype is optional, defaulting to dtype='d'.

        dia_matrix((data, offsets), shape=(M, N))
            where the ``data[k,:]`` stores the diagonal entries for
            diagonal ``offsets[k]`` (See example below)

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
        DIA format data array of the matrix
    offsets
        DIA format offset array of the matrix

    Notes
    -----

    Sparse matrices can be used in arithmetic operations: they support
    addition, subtraction, multiplication, division, and matrix power.

    Examples
    --------

    >>> import numpy as np
    >>> from scipy.sparse import dia_matrix
    >>> dia_matrix((3, 4), dtype=np.int8).toarray()
    array([[0, 0, 0, 0],
           [0, 0, 0, 0],
           [0, 0, 0, 0]], dtype=int8)

    >>> data = np.array([[1, 2, 3, 4]]).repeat(3, axis=0)
    >>> offsets = np.array([0, -1, 2])
    >>> dia_matrix((data, offsets), shape=(4, 4)).toarray()
    array([[1, 0, 3, 0],
           [1, 2, 0, 4],
           [0, 2, 3, 0],
           [0, 0, 3, 4]])

    """

    def __init__(self, arg1, shape=None, dtype=None, copy=False):
        _data_matrix.__init__(self)

        if isspmatrix_dia(arg1):
            if copy:
                arg1 = arg1.copy()
            self.data = arg1.data
            self.offsets = arg1.offsets
            self.shape = arg1.shape
        elif isspmatrix(arg1):
            if isspmatrix_dia(arg1) and copy:
                A = arg1.copy()
            else:
                A = arg1.todia()
            self.data = A.data
            self.offsets = A.offsets
            self.shape = A.shape
        elif isinstance(arg1, tuple):
            if isshape(arg1):
                # It's a tuple of matrix dimensions (M, N)
                # create empty matrix
                self.shape = arg1   # spmatrix checks for errors here
                self.data = np.zeros((0,0), getdtype(dtype, default=float))
                idx_dtype = get_index_dtype(maxval=max(self.shape))
                self.offsets = np.zeros((0), dtype=idx_dtype)
            else:
                try:
                    # Try interpreting it as (data, offsets)
                    data, offsets = arg1
                except:
                    raise ValueError('unrecognized form for dia_matrix constructor')
                else:
                    if shape is None:
                        raise ValueError('expected a shape argument')
                    self.data = np.atleast_2d(np.array(arg1[0], dtype=dtype, copy=copy))
                    self.offsets = np.atleast_1d(np.array(arg1[1],
                                                          dtype=get_index_dtype(maxval=max(shape)),
                                                          copy=copy))
                    self.shape = shape
        else:
            #must be dense, convert to COO first, then to DIA
            try:
                arg1 = np.asarray(arg1)
            except:
                raise ValueError("unrecognized form for"
                        " %s_matrix constructor" % self.format)
            from .coo import coo_matrix
            A = coo_matrix(arg1, dtype=dtype).todia()
            self.data = A.data
            self.offsets = A.offsets
            self.shape = A.shape

        if dtype is not None:
            self.data = self.data.astype(dtype)

        #check format
        if self.offsets.ndim != 1:
            raise ValueError('offsets array must have rank 1')

        if self.data.ndim != 2:
            raise ValueError('data array must have rank 2')

        if self.data.shape[0] != len(self.offsets):
            raise ValueError('number of diagonals (%d) '
                    'does not match the number of offsets (%d)'
                    % (self.data.shape[0], len(self.offsets)))

        if len(np.unique(self.offsets)) != len(self.offsets):
            raise ValueError('offset array contains duplicate values')

    def __repr__(self):
        nnz = self.getnnz()
        format = self.getformat()
        return "<%dx%d sparse matrix of type '%s'\n" \
               "\twith %d stored elements (%d diagonals) in %s format>" % \
               (self.shape + (self.dtype.type, nnz, self.data.shape[0],
                 _formats[format][1],))

    def getnnz(self):
        """number of nonzero values

        explicit zero values are included in this number
        """
        M,N = self.shape
        nnz = 0
        for k in self.offsets:
            if k > 0:
                nnz += min(M,N-k)
            else:
                nnz += min(M+k,N)
        return int(nnz)

    nnz = property(fget=getnnz)

    def _mul_vector(self, other):
        x = other

        y = np.zeros(self.shape[0], dtype=upcast_char(self.dtype.char,
                                                       x.dtype.char))

        L = self.data.shape[1]

        M,N = self.shape

        dia_matvec(M,N, len(self.offsets), L, self.offsets, self.data, x.ravel(), y.ravel())

        return y

    def _mul_multimatrix(self, other):
        return np.hstack([self._mul_vector(col).reshape(-1,1) for col in other.T])

    def _setdiag(self, values, k=0):
        M, N = self.shape

        if values.ndim == 0:
            # broadcast
            values_n = np.inf
        else:
            values_n = len(values)

        if k < 0:
            n = min(M + k, N, values_n)
            min_index = 0
            max_index = n
        else:
            n = min(M, N - k, values_n)
            min_index = k
            max_index = k + n

        if values.ndim != 0:
            # allow also longer sequences
            values = values[:n]

        if k in self.offsets:
            self.data[self.offsets == k, min_index:max_index] = values
        else:
            self.offsets = np.append(self.offsets, self.offsets.dtype.type(k))
            m = max(max_index, self.data.shape[1])
            data = np.zeros((self.data.shape[0]+1, m), dtype=self.data.dtype)
            data[:-1,:self.data.shape[1]] = self.data
            data[-1, min_index:max_index] = values
            self.data = data

    def todia(self,copy=False):
        if copy:
            return self.copy()
        else:
            return self

    def tocsr(self):
        #this could be faster
        return self.tocoo().tocsr()

    def tocsc(self):
        #this could be faster
        return self.tocoo().tocsc()

    def tocoo(self):
        num_data = len(self.data)
        len_data = self.data.shape[1]

        row = np.arange(len_data).reshape(1,-1).repeat(num_data,axis=0)
        col = row.copy()

        for i,k in enumerate(self.offsets):
            row[i,:] -= k

        row,col,data = row.ravel(),col.ravel(),self.data.ravel()

        mask = (row >= 0)
        mask &= (row < self.shape[0])
        mask &= (col < self.shape[1])
        mask &= data != 0
        row,col,data = row[mask],col[mask],data[mask]

        from .coo import coo_matrix
        return coo_matrix((data,(row,col)), shape=self.shape)

    # needed by _data_matrix
    def _with_data(self, data, copy=True):
        """Returns a matrix with the same sparsity structure as self,
        but with different data.  By default the structure arrays are copied.
        """
        if copy:
            return dia_matrix((data, self.offsets.copy()), shape=self.shape)
        else:
            return dia_matrix((data,self.offsets), shape=self.shape)


def isspmatrix_dia(x):
    return isinstance(x, dia_matrix)
