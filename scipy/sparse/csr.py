"""Compressed Sparse Row sparse matrix format
"""

__all__ = ['csr_matrix', 'isspmatrix_csr']


from warnings import warn

import numpy
from numpy import array, matrix, asarray, asmatrix, zeros, rank, intc, \
        empty, hstack, isscalar, ndarray, shape, searchsorted

from base import spmatrix,isspmatrix
from sparsetools import csr_tocsc
from sputils import upcast, to_native, isdense, isshape, getdtype, \
        isscalarlike

from compressed import _cs_matrix,resize1d

class csr_matrix(_cs_matrix):
    """ Compressed sparse row matrix
        This can be instantiated in several ways:
          - csr_matrix(d)
            with a dense matrix d

          - csr_matrix(s)
            with another sparse matrix s (sugar for .tocsr())

          - csr_matrix((M, N), [dtype])
            to construct a container, where (M, N) are dimensions and
            dtype is optional, defaulting to dtype='d'.

          - csr_matrix((data, ij), [shape=(M, N)])
            where data, ij satisfy:
                a[ij[0, k], ij[1, k]] = data[k]

          - csr_matrix((data, col, ptr), [shape=(M, N)])
            standard CSR representation
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
        return _cs_matrix._transpose(self, csc_matrix, copy)

    def __getitem__(self, key):
        if isinstance(key, tuple):
            #TODO use _swap() to unify this in _cs_matrix
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
            #this was implemented in fortran before - is there a noticable performance advantage?
            indxs = numpy.where(col == self.indices[self.indptr[row]:self.indptr[row+1]])
            if len(indxs[0]) == 0:
                return 0
            else:
                return self.data[self.indptr[row]:self.indptr[row+1]][indxs[0]]
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

    def _otherformat(self):
        return "csc"

    def _toother(self):
        return self.tocsc()

    def _tothis(self, other):
        return other.tocsr()


from sputils import _isinstance

def isspmatrix_csr(x):
    return _isinstance(x, csr_matrix)

