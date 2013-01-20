"""Compressed Sparse Column matrix format"""
from __future__ import division, print_function, absolute_import

__docformat__ = "restructuredtext en"

__all__ = ['csc_matrix', 'isspmatrix_csc']

from warnings import warn

import numpy as np
from scipy.lib.six.moves import xrange

from .sparsetools import csc_tocsr
from .sputils import upcast, isintlike

from .compressed import _cs_matrix


class csc_matrix(_cs_matrix):
    """
    Compressed Sparse Column matrix

    This can be instantiated in several ways:

        csc_matrix(D)
            with a dense matrix or rank-2 ndarray D

        csc_matrix(S)
            with another sparse matrix S (equivalent to S.tocsc())

        csc_matrix((M, N), [dtype])
            to construct an empty matrix with shape (M, N)
            dtype is optional, defaulting to dtype='d'.

        csc_matrix((data, ij), [shape=(M, N)])
            where ``data`` and ``ij`` satisfy the relationship
            ``a[ij[0, k], ij[1, k]] = data[k]``

        csc_matrix((data, indices, indptr), [shape=(M, N)])
            is the standard CSC representation where the row indices for
            column i are stored in ``indices[indptr[i]:indptr[i+1]]``
            and their corresponding values are stored in
            ``data[indptr[i]:indptr[i+1]]``.  If the shape parameter is
            not supplied, the matrix dimensions are inferred from
            the index arrays.

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
        Data array of the matrix
    indices
        CSC format index array
    indptr
        CSC format index pointer array
    has_sorted_indices
        Whether indices are sorted

    Notes
    -----

    Sparse matrices can be used in arithmetic operations: they support
    addition, subtraction, multiplication, division, and matrix power.

    Advantages of the CSC format
        - efficient arithmetic operations CSC + CSC, CSC * CSC, etc.
        - efficient column slicing
        - fast matrix vector products (CSR, BSR may be faster)

    Disadvantages of the CSC format
      - slow row slicing operations (consider CSR)
      - changes to the sparsity structure are expensive (consider LIL or DOK)


    Examples
    --------

    >>> from scipy.sparse import *
    >>> from scipy import *
    >>> csc_matrix( (3,4), dtype=int8 ).todense()
    matrix([[0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]], dtype=int8)

    >>> row = array([0,2,2,0,1,2])
    >>> col = array([0,0,1,2,2,2])
    >>> data = array([1,2,3,4,5,6])
    >>> csc_matrix( (data,(row,col)), shape=(3,3) ).todense()
    matrix([[1, 0, 4],
            [0, 0, 5],
            [2, 3, 6]])

    >>> indptr = array([0,2,3,6])
    >>> indices = array([0,2,2,0,1,2])
    >>> data = array([1,2,3,4,5,6])
    >>> csc_matrix( (data,indices,indptr), shape=(3,3) ).todense()
    matrix([[1, 0, 4],
            [0, 0, 5],
            [2, 3, 6]])

    """

    def transpose(self, copy=False):
        from .csr import csr_matrix
        M,N = self.shape
        return csr_matrix((self.data,self.indices,self.indptr),(N,M),copy=copy)

    def __iter__(self):
        csr = self.tocsr()
        for r in xrange(self.shape[0]):
            yield csr[r,:]

    def tocsc(self, copy=False):
        if copy:
            return self.copy()
        else:
            return self

    def tocsr(self):
        M,N = self.shape
        indptr  = np.empty(M + 1,    dtype=np.intc)
        indices = np.empty(self.nnz, dtype=np.intc)
        data    = np.empty(self.nnz, dtype=upcast(self.dtype))

        csc_tocsr(M, N, \
                 self.indptr, self.indices, self.data, \
                 indptr, indices, data)

        from .csr import csr_matrix
        A = csr_matrix((data, indices, indptr), shape=self.shape)
        A.has_sorted_indices = True
        return A


    def __getitem__(self, key):
        # use CSR to implement fancy indexing
        if isinstance(key, tuple):
            row = key[0]
            col = key[1]

            if isintlike(row) or isinstance(row, slice):
                return self.T[col,row].T
            else:
                #[[1,2],??] or [[[1],[2]],??]
                if isintlike(col) or isinstance(col,slice):
                    return self.T[col,row].T
                else:
                    row = np.asarray(row, dtype=np.intc)
                    col = np.asarray(col, dtype=np.intc)
                    if len(row.shape) == 1:
                        return self.T[col,row]
                    elif len(row.shape) == 2:
                        row = row.reshape(-1)
                        col = col.reshape(-1,1)
                        return self.T[col,row].T
                    else:
                        raise NotImplementedError('unsupported indexing')

            return self.T[col,row].T
        elif isintlike(key) or isinstance(key,slice):
            return self.T[:,key].T                              #[i] or [1:2]
        else:
            return self.T[:,key].T                              #[[1,2]]

    def getrow(self, i):
        """Returns a copy of row i of the matrix, as a (1 x n)
        CSR matrix (row vector).
        """
        # transpose to use CSR code
        # we convert to CSR to maintain compatibility with old impl.
        # in spmatrix.getrow()
        return self.T.getcol(i).T.tocsr()

    def getcol(self, i):
        """Returns a copy of column i of the matrix, as a (m x 1)
        CSC matrix (column vector).
        """
        # transpose to use CSR code
        return self.T.getrow(i).T

    # these functions are used by the parent class (_cs_matrix)
    # to remove redudancy between csc_matrix and csr_matrix
    def _swap(self,x):
        """swap the members of x if this is a column-oriented matrix
        """
        return (x[1],x[0])


def isspmatrix_csc(x):
    return isinstance(x, csc_matrix)
