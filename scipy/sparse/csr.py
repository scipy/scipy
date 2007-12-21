"""Compressed Sparse Row matrix format
"""

__all__ = ['csr_matrix', 'isspmatrix_csr']


from warnings import warn

import numpy
from numpy import array, matrix, asarray, asmatrix, zeros, rank, intc, \
        empty, hstack, isscalar, ndarray, shape, searchsorted, where

from base import spmatrix,isspmatrix
from sparsetools import csr_tocsc
from sputils import upcast, to_native, isdense, isshape, getdtype, \
        isscalarlike

from compressed import _cs_matrix,resize1d

class csr_matrix(_cs_matrix):
    """Compressed Sparse Row matrix

    This can be instantiated in several ways:
      - csr_matrix(D)
        with a dense matrix or rank-2 ndarray D

      - csr_matrix(S)
        with another sparse matrix S (equivalent to S.tocsr())

      - csr_matrix((M, N), [dtype])
        to construct an empty matrix with shape (M, N)
        dtype is optional, defaulting to dtype='d'.

      - csr_matrix((data, ij), [shape=(M, N)])
        where data, ij satisfy:
            a[ij[0, k], ij[1, k]] = data[k]

      - csr_matrix((data, indices, indptr), [shape=(M, N)])
        is the native CSR representation where:
            the column indices for row i are stored in
                indices[ indptr[i]: indices[i+1] ] 
            and their corresponding values are stored in
                data[ indptr[i]: indptr[i+1] ]
        If the shape parameter is not supplied, the matrix dimensions
        are inferred from the index arrays.


    *Examples*
    ----------

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

    def __getitem__(self, key):
        #TODO unify these in _cs_matrix
        if isinstance(key, tuple):
            row = key[0]
            col = key[1]
            
            if isinstance(row, slice):
                # Returns a new matrix!
                return self.get_submatrix( row, col )
            elif isinstance(col, slice):
                return self._getslice(row, col)
            
            M, N = self.shape
            if (row < 0):
                row = M + row
            if (col < 0):
                col = N + col
            if not (0<=row<M) or not (0<=col<N):
                raise IndexError, "index out of bounds"
            
            major_index, minor_index = self._swap((row,col))

            start = self.indptr[major_index]
            end   = self.indptr[major_index+1]
            indxs = where(minor_index == self.indices[start:end])[0]

            num_matches = len(indxs)

            if num_matches == 0:
                # entry does not appear in the matrix
                return 0
            elif num_matches == 1:
                return self.data[start:end][indxs[0]]
            else:
                raise ValueError,'nonzero entry (%d,%d) occurs more than once' % (row,col)

        elif isintlike(key):
            return self[key, :]
        else:
            raise IndexError, "invalid index"

    def _getslice(self, i, myslice):
        return self._getrowslice(i, myslice)

    def _getrowslice(self, i, myslice):
        """Returns a view of the elements [i, myslice.start:myslice.stop].
        """
        start, stop, stride = myslice.indices(self.shape[1])
        return _cs_matrix._get_slice(self, i, start, stop, stride, (1, stop-start))

    def __setitem__(self, key, val):
        if isinstance(key, tuple):
            row = key[0]
            col = key[1]
            if not (isscalarlike(row) and isscalarlike(col)):
                raise NotImplementedError("Fancy indexing in assignment not "
                                          "supported for csr matrices.")
            M, N = self.shape
            if (row < 0):
                row = M + row
            if (col < 0):
                col = N + col
            if (row < 0) or (col < 0):
                raise IndexError, "index out of bounds"
            if (row >= M):
                self.indptr = resize1d(self.indptr, row+2)
                self.indptr[M+1:] = self.indptr[M]
                M = row+1
            if (col >= N):
                N = col+1
            self.shape = (M, N)

            indxs = numpy.where(col == self.indices[self.indptr[row]:self.indptr[row+1]])
            if len(indxs[0]) == 0:
                #value not present
                self.data    = resize1d(self.data, self.nnz + 1)
                self.indices = resize1d(self.indices, self.nnz + 1)

                newindex = self.indptr[row]
                self.data[newindex+1:]   = self.data[newindex:-1]
                self.indices[newindex+1:] = self.indices[newindex:-1]

                self.data[newindex]   = val
                self.indices[newindex] = col
                self.indptr[row+1:] += 1

            elif len(indxs[0]) == 1:
                #value already present
                self.data[self.indptr[row]:self.indptr[row+1]][indxs[0]] = val
            else:
                raise IndexError, "row index occurs more than once"

            self.check_format(full_check=False)
        else:
            # We should allow slices here!
            raise IndexError, "invalid index"

    def rowcol(self, ind):
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
        if blocksize in [None, (1,1)]:
            from bsr import bsr_matrix
            arg1 = (self.data.reshape(-1,1,1),self.indices,self.indptr)  
            return bsr_matrix( arg1, shape=self.shape, copy=copy )
        else:
            #TODO make this more efficient
            return self.tocoo(copy=False).tobsr(blocksize=blocksize)
    
#    def tobsc(self,blocksize=None):
#        if blocksize in [None, (1,1)]:
#            from bsc import bsc_matrix
#            csc = self.tocsc()
#            arg1 = (csc.data.reshape(-1,1,1),csc.indices,csc.indptr)  
#            return bsc_matrix( arg1, shape=self.shape )
#        else:
#            #TODO make this more efficient
#            return self.tocoo(copy=False).tobsc(blocksize=blocksize)
    
    def get_submatrix( self, slice0, slice1 ):
        """Return a submatrix of this matrix (new matrix is created).
        Contigous range of rows and columns can be selected using:
        1. a slice object
        2. a tuple (from, to)
        3. a scalar for single row/column selection."""
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

