"""Compressed Sparse Column sparse matrix format
"""

__all__ = ['csc_matrix', 'isspmatrix_csc']

from warnings import warn

import numpy
from numpy import array, matrix, asarray, asmatrix, zeros, rank, intc, \
        empty, hstack, isscalar, ndarray, shape, searchsorted, where

from base import spmatrix,isspmatrix
from sparsetools import csc_tocsr
from sputils import upcast, to_native, isdense, isshape, getdtype, \
        isscalarlike

from compressed import _cs_matrix,resize1d


class csc_matrix(_cs_matrix):
    """ Compressed sparse column matrix
        This can be instantiated in several ways:
          - csc_matrix(d)
            with a dense matrix d

          - csc_matrix(s)
            with another sparse matrix s (sugar for .tocsc())

          - csc_matrix((M, N), [dtype])
            to construct a container, where (M, N) are dimensions and
            dtype is optional, defaulting to dtype='d'.

          - csc_matrix((data, ij), [(M, N)])
            where data, ij satisfy:
                a[ij[0, k], ij[1, k]] = data[k]

          - csc_matrix((data, row, ptr), [(M, N)])
            standard CSC representation
    """
    
    def __getattr__(self, attr):
        if attr == 'rowind':
            warn("rowind attribute no longer in use. Use .indices instead",
                    DeprecationWarning)
            return self.indices
        else:
            return _cs_matrix.__getattr__(self, attr)

    def __iter__(self):
        csr = self.tocsr()
        for r in xrange(self.shape[0]):
            yield csr[r,:]

    def transpose(self, copy=False):
        from csr import csr_matrix
        return _cs_matrix._transpose(self, csr_matrix, copy)

    def __getitem__(self, key):
        #TODO unify these in _cs_matrix
        if isinstance(key, tuple):
            row = key[0]
            col = key[1]
            
            if isinstance(col, slice):
                # Returns a new matrix!
                return self.get_submatrix( row, col )
            elif isinstance(row, slice):
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
            # Was: return self.data[key]
            # If this is allowed, it should return the relevant row, as for
            # dense matrices (and __len__ should be supported again).
            raise IndexError, "integer index not supported for csc_matrix"
        else:
            raise IndexError, "invalid index"
    
    def _getslice(self, i, myslice):
        return self._getcolslice(i, myslice)

    def _getcolslice(self, myslice, j):
        """Returns a view of the elements [myslice.start:myslice.stop, j].
        """
        start, stop, stride = myslice.indices(self.shape[0])
        return _cs_matrix._get_slice(self, j, start, stop, stride, (stop - start, 1))



    def __setitem__(self, key, val):
        if isinstance(key, tuple):
            row = key[0]
            col = key[1]
            if not (isscalarlike(row) and isscalarlike(col)):
                raise NotImplementedError("Fancy indexing in assignments not"
                                          "supported for csc matrices.")
            M, N = self.shape
            if (row < 0):
                row = M + row
            if (col < 0):
                col = N + col
            if (row < 0) or (col < 0):
                raise IndexError, "index out of bounds"
            if (col >= N):
                self.indptr = resize1d(self.indptr, col+2)
                self.indptr[N+1:] = self.indptr[N]
                N = col+1
            if (row >= M):
                M = row+1
            self.shape = (M, N)

            indxs = numpy.where(row == self.indices[self.indptr[col]:self.indptr[col+1]])
    
            if len(indxs[0]) == 0:
                #value not present
                #TODO handle this with concatenation
                self.data    = resize1d(self.data,    self.nnz + 1)
                self.indices = resize1d(self.indices, self.nnz + 1)

                newindex = self.indptr[col]
                self.data[newindex+1:]    = self.data[newindex:-1]
                self.indices[newindex+1:] = self.indices[newindex:-1]

                self.data[newindex]   = val
                self.indices[newindex] = row
                self.indptr[col+1:] += 1

            elif len(indxs[0]) == 1:
                #value already present
                self.data[self.indptr[col]:self.indptr[col+1]][indxs[0]] = val
            else:
                raise IndexError, "row index occurs more than once"

            self.check_format(full_check=False)
        else:
            # We should allow slices here!
            raise IndexError, "invalid index"

    def rowcol(self, ind):
        row = self.indices[ind]
        col = searchsorted(self.indptr, ind+1)-1
        return (row, col)

    def tocsc(self, copy=False):
        if copy:
            return self.copy()
        else:
            return self
    
    def tocsr(self):
        indptr  = empty(self.shape[0] + 1, dtype=intc)
        indices = empty(self.nnz, dtype=intc)
        data    = empty(self.nnz, dtype=upcast(self.dtype))

        csc_tocsr(self.shape[0], self.shape[1], \
                 self.indptr, self.indices, self.data, \
                 indptr, indices, data)

        from csr import csr_matrix
        return csr_matrix((data, indices, indptr), self.shape)
    
    def toarray(self):
        A = self.tocoo(copy=False)
        M = zeros(self.shape, dtype=self.dtype)
        M[A.row, A.col] = A.data
        return M

    def get_submatrix( self, slice0, slice1 ):
        """Return a submatrix of this matrix (new matrix is created).
        Contigous range of rows and columns can be selected using:
        1. a slice object
        2. a tuple (from, to)
        3. a scalar for single row/column selection."""
        aux = _cs_matrix._get_submatrix( self, self.shape[1], self.shape[0],
                                         slice1, slice0 )
        nr, nc = aux[3:]
        return self.__class__( aux[:3], shape = (nc, nr) )
    
    # these functions are used by the parent class (_cs_matrix)
    # to remove redudancy between csc_matrix and csr_matrix
    def _swap(self,x):
        """swap the members of x if this is a column-oriented matrix
        """
        return (x[1],x[0])

    def _otherformat(self):
        return "csr"

    def _toother(self):
        return self.tocsr()

    def _tothis(self, other):
        return other.tocsc()

from sputils import _isinstance

def isspmatrix_csc(x):
    return _isinstance(x, csc_matrix)

