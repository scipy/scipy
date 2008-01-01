"""Compressed Sparse Row matrix format"""

__all__ = ['csr_matrix', 'isspmatrix_csr']


from warnings import warn

import numpy
from numpy import array, matrix, asarray, asmatrix, zeros, rank, intc, \
        empty, hstack, isscalar, ndarray, shape, searchsorted, where, \
        concatenate

from base import spmatrix, isspmatrix
from sparsetools import csr_tocsc
from sputils import upcast, to_native, isdense, isshape, getdtype, \
        isscalarlike, isintlike

from compressed import _cs_matrix

class csr_matrix(_cs_matrix):
    """Compressed Sparse Row matrix

    This can be instantiated in several ways:
      - csr_matrix(D)
        - with a dense matrix or rank-2 ndarray D

      - csr_matrix(S)
        - with another sparse matrix S (equivalent to S.tocsr())

      - csr_matrix((M, N), [dtype])
        - to construct an empty matrix with shape (M, N)
        - dtype is optional, defaulting to dtype='d'.

      - csr_matrix((data, ij), [shape=(M, N)])
        - where data, ij satisfy:
          - a[ij[0, k], ij[1, k]] = data[k]

      - csr_matrix((data, indices, indptr), [shape=(M, N)])
        - is the standard CSR representation where
          the column indices for row i are stored in
           - indices[ indptr[i]: indices[i+1] ] 
          and their corresponding values are stored in
           - data[ indptr[i]: indptr[i+1] ]
        - If the shape parameter is not supplied, the matrix dimensions
          are inferred from the index arrays.


    Notes
    =====
        Advantages of the CSR format
        ----------------------------
          - efficient arithmetic operations CSR + CSR, CSR * CSR, etc.
          - efficient row slicing
          - fast matrix vector products
        
        Disadvantages of the CSR format
        -------------------------------
          - slow column slicing operations (prefer CSC)
          - changes to the sparsity structure are expensive (prefer LIL, DOK)



    Examples
    ========    

    >>> from scipy.sparse import *
    >>> from scipy import *
    >>> csr_matrix( (3,4), dtype='i' ).todense()
    matrix([[0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]])

    >>> row = array([0,0,1,2,2,2])
    >>> col = array([0,2,2,0,1,2])
    >>> data = array([1,2,3,4,5,6])
    >>> csr_matrix( (data,(row,col)), shape=(3,3) ).todense()
    matrix([[1, 0, 2],
            [0, 0, 3],
            [4, 5, 6]])
    
    >>> indptr = array([0,2,3,6])
    >>> indices = array([0,2,2,0,1,2])
    >>> data = array([1,2,3,4,5,6])
    >>> csr_matrix( (data,indices,indptr), shape=(3,3) ).todense()
    matrix([[1, 0, 2],
            [0, 0, 3],
            [4, 5, 6]])

    """

    def __getattr__(self, attr):
        if attr == 'colind':
            warn("colind attribute no longer in use. Use .indices instead",
                    DeprecationWarning)
            return self.indices
        else:
            return _cs_matrix.__getattr__(self, attr)

    def transpose(self, copy=False):
        from csc import csc_matrix
        M,N = self.shape
        return csc_matrix((self.data,self.indices,self.indptr),(N,M),copy=copy)


    def rowcol(self, ind):
        #TODO remove after 0.7
        warn('rowcol() is deprecated',DeprecationWarning)
        col = self.indices[ind]
        row = searchsorted(self.indptr, ind+1)-1
        return (row, col)


    def tolil(self):
        from lil import lil_matrix
        lil = lil_matrix(self.shape,dtype=self.dtype)
     
        csr = self.sorted_indices() #lil_matrix needs sorted rows
        
        rows,data = lil.rows,lil.data
        ptr,ind,dat = csr.indptr,csr.indices,csr.data

        for n in xrange(self.shape[0]):
            start = ptr[n]
            end   = ptr[n+1]
            rows[n] = ind[start:end].tolist()
            data[n] = dat[start:end].tolist()

        return lil

    def tocsr(self, copy=False):
        if copy:
            return self.copy()
        else:
            return self

    def tocsc(self):
        indptr  = empty(self.shape[1] + 1, dtype=intc)
        indices = empty(self.nnz, dtype=intc)
        data    = empty(self.nnz, dtype=upcast(self.dtype))

        csr_tocsc(self.shape[0], self.shape[1], \
                  self.indptr, self.indices, self.data, \
                  indptr, indices, data)

        from csc import csc_matrix
        return csc_matrix((data, indices, indptr), self.shape)

    def tobsr(self,blocksize=None,copy=True):
        if blocksize == (1,1):
            from bsr import bsr_matrix
            arg1 = (self.data.reshape(-1,1,1),self.indices,self.indptr)  
            return bsr_matrix( arg1, shape=self.shape, copy=copy )
        else:
            #TODO make this more efficient
            return self.tocoo(copy=False).tobsr(blocksize=blocksize)
    
    
    def get_submatrix( self, slice0, slice1 ):
        """Return a submatrix of this matrix (new matrix is created).
        Contigous range of rows and columns can be selected using:
          1. a slice object
          2. a tuple (from, to)
          3. a scalar for single row/column selection.

        """

        aux = _cs_matrix._get_submatrix( self, self.shape[0], self.shape[1],
                                         slice0, slice1 )
        nr, nc = aux[3:]
        return self.__class__( aux[:3], shape = (nr, nc) )

    # these functions are used by the parent class (_cs_matrix)
    # to remove redudancy between csc_matrix and csr_matrix
    def _swap(self,x):
        """swap the members of x if this is a column-oriented matrix
        """
        return (x[0],x[1])


from sputils import _isinstance

def isspmatrix_csr(x):
    return _isinstance(x, csr_matrix)

