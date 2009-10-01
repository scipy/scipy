"""Compressed Sparse Column matrix format"""

__docformat__ = "restructuredtext en"

__all__ = ['csc_matrix', 'isspmatrix_csc']

from warnings import warn

import numpy as np

from sparsetools import csc_tocsr
from sputils import upcast, isintlike

from compressed import _cs_matrix


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
            column i are stored in ``indices[indptr[i]:indices[i+1]]``
            and their corresponding values are stored in
            ``data[indptr[i]:indptr[i+1]]``.  If the shape parameter is
            not supplied, the matrix dimensions are inferred from
            the index arrays.

    Notes
    -----
    Advantages of the CSC format
        - efficient arithmetic operations CSC + CSC, CSC * CSC, etc.
        - efficient column slicing
        - fast matrix vector products (CSR, BSR may be faster)

    Disadvantages of the CSC format
      - slow row slicing operations (consider CSR)
      - changes to the sparsity structure are expensive (consider LIL or DOK)


    Examples
    ========

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

    def __getattr__(self, attr):
        if attr == 'rowind':
            warn("rowind attribute no longer in use. Use .indices instead",
                    DeprecationWarning)
            return self.indices
        else:
            return _cs_matrix.__getattr__(self, attr)

    def transpose(self, copy=False):
        from csr import csr_matrix
        M,N = self.shape
        return csr_matrix((self.data,self.indices,self.indptr),(N,M),copy=copy)

    def __iter__(self):
        csr = self.tocsr()
        for r in xrange(self.shape[0]):
            yield csr[r,:]

    @np.deprecate
    def rowcol(self, ind):
        #TODO remove after 0.7
        row = self.indices[ind]
        col = np.searchsorted(self.indptr, ind+1) - 1
        return (row, col)

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

        from csr import csr_matrix
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


    # these functions are used by the parent class (_cs_matrix)
    # to remove redudancy between csc_matrix and csr_matrix
    def _swap(self,x):
        """swap the members of x if this is a column-oriented matrix
        """
        return (x[1],x[0])


from sputils import _isinstance

def isspmatrix_csc(x):
    return _isinstance(x, csc_matrix)
