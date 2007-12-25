"""Compressed Block Sparse Row matrix format
"""

__all__ = ['bsr_matrix', 'isspmatrix_bsr']

from numpy import zeros, intc, array, asarray, arange, diff, tile, rank, \
        prod, ravel, empty, matrix, asmatrix, empty_like, hstack

from sparsetools import bsr_matvec, csr_matmat_pass1, bsr_matmat_pass2
from block import _block_matrix 
from base import isspmatrix
from sputils import isdense, upcast, isscalarlike

class bsr_matrix(_block_matrix):
    """Block Sparse Row matrix

    This can be instantiated in several ways:
      - bsr_matrix(D, [blocksize=(R,C)])
        with a dense matrix or rank-2 ndarray D

      - bsr_matrix(S, [blocksize=(R,C)])
        with another sparse matrix S (equivalent to S.tocsr())

      - bsr_matrix((M, N), [blocksize=(R,C), dtype])
        to construct an empty matrix with shape (M, N)
        dtype is optional, defaulting to dtype='d'.

      - bsr_matrix((data, ij), [blocksize=(R,C), shape=(M, N)])
        where data, ij satisfy:
            a[ij[0, k], ij[1, k]] = data[k]

      - bsr_matrix((data, indices, indptr), [shape=(M, N)])
        is the standard BSR representation where:
            the block column indices for row i are stored in
                indices[ indptr[i]: indices[i+1] ] 
            and their corresponding block values are stored in
                data[ indptr[i]: indptr[i+1] ]
        If the shape parameter is not supplied, the matrix dimensions
        are inferred from the index arrays.


    *Notes*
    -------
        The blocksize (R,C) must evenly divide the shape of 
        the matrix (M,N).  That is, R and C must satisfy the
        relationship M % R = 0 and N % C = 0.
    
        The Block Compressed Row (BSR) format is very similar to the
        Compressed Sparse Row (CSR) format.  BSR is appropriate for
        sparse matrices with dense sub matrices like the last example
        below.  Such matrices often arise, for instance, in finite
        element discretizations.


    *Examples*
    ----------

    >>> from scipy.sparse import *
    >>> from scipy import *
    >>> bsr_matrix( (3,4), dtype='i' ).todense()
    matrix([[0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]])

    >>> row = array([0,0,1,2,2,2])
    >>> col = array([0,2,2,0,1,2])
    >>> data = kron([1,2,3,4,5,6])
    >>> bsr_matrix( (data,(row,col)), shape=(3,3) ).todense()
    matrix([[1, 0, 2],
            [0, 0, 3],
            [4, 5, 6]])
    
    >>> indptr = array([0,2,3,6])
    >>> indices = array([0,2,2,0,1,2])
    >>> data = array([1,2,3,4,5,6]).repeat(4).reshape(6,2,2)
    >>> bsr_matrix( (data,indices,indptr), shape=(6,6) ).todense()
    matrix([[1, 1, 0, 0, 2, 2],
            [1, 1, 0, 0, 2, 2],
            [0, 0, 0, 0, 3, 3],
            [0, 0, 0, 0, 3, 3],
            [4, 4, 5, 5, 6, 6],
            [4, 4, 5, 5, 6, 6]])
    
    """


    def __mul__(self, other): # self * other
        """ Scalar, vector, or matrix multiplication
        """
        if isscalarlike(other):
            return self._with_data(self.data * other)
        else:
            return self.dot(other)

    def matvec(self, other, output=None):
        """Sparse matrix vector product (self * other)

        'other' may be a rank 1 array of length N or a rank 2 array 
        or matrix with shape (N,1).  
        
        If the optional 'output' parameter is defined, it will
        be used to store the result.  Otherwise, a new vector
        will be allocated.
             
        """
        if isdense(other):
            M,N = self.shape
            X,Y = self.blocksize

            if other.shape != (N,) and other.shape != (N,1):
                raise ValueError, "dimension mismatch"
    
            #output array
            if output is None:
                y = empty( self.shape[0], dtype=upcast(self.dtype,other.dtype) )
            else:
                if output.shape != (M,) and output.shape != (M,1):
                    raise ValueError, "output array has improper dimensions"
                if not output.flags.c_contiguous:
                    raise ValueError, "output array must be contiguous"
                if output.dtype != upcast(self.dtype,other.dtype):
                    raise ValueError, "output array has dtype=%s "\
                            "dtype=%s is required" % \
                            (output.dtype,upcast(self.dtype,other.dtype))
                y = output
            
            
            bsr_matvec(M/X, N/Y, X, Y, \
                self.indptr, self.indices, ravel(self.data), ravel(other), y)

            if isinstance(other, matrix):
                y = asmatrix(y)

            if other.ndim == 2 and other.shape[1] == 1:
                # If 'other' was an (nx1) column vector, reshape the result
                y = y.reshape(-1,1)

            return y

        elif isspmatrix(other):
            raise TypeError, "use matmat() for sparse * sparse"

        else:
            raise TypeError, "need a dense vector"


    def matmat(self, other):
        if isspmatrix(other):
            M, K1 = self.shape
            K2, N = other.shape
            if (K1 != K2):
                raise ValueError, "shape mismatch error"

            indptr = empty_like( self.indptr )
            
            R,n = self.blocksize

            if isspmatrix_bsr(other):
                C = other.blocksize[1]
            else:
                C = 1
            other = other.tobsr(blocksize=(n,C)) #convert to this format

            csr_matmat_pass1( M/R, N/C, \
                    self.indptr,  self.indices, \
                    other.indptr, other.indices, \
                    indptr)

            bnnz = indptr[-1]
            indices = empty( bnnz, dtype=intc)
            data    = empty( R*C*bnnz, dtype=upcast(self.dtype,other.dtype))

            bsr_matmat_pass2( M/R, N/C, R, C, n, \
                    self.indptr,  self.indices,  ravel(self.data), \
                    other.indptr, other.indices, ravel(other.data), \
                    indptr,       indices,       data)
            
            data = data.reshape(-1,R,C)
            #TODO eliminate zeros

            return bsr_matrix((data,indices,indptr),shape=(M,N),blocksize=(R,C))
        elif isdense(other):
            # TODO make sparse * dense matrix multiplication more efficient

            # matvec each column of other
            return hstack( [ self * col.reshape(-1,1) for col in other.T ] )
        else:
            raise TypeError, "need a dense or sparse matrix"

#    def transpose(self,copy=False):
#        M,N = self.shape
#            
#        data    = self.data.swapaxes(1,2)
#        indices = self.indices
#        indptr  = self.indptr
#
#        from bsc import bsc_matrix
#        return bsc_matrix( (data,indices,indptr), shape=(N,M), copy=copy)
   
    def tocoo(self,copy=True):
        """Convert this matrix to COOrdinate format.

        When copy=False the data array will be shared between
        this matrix and the resultant coo_matrix.
        """
        
        M,N = self.shape
        X,Y = self.blocksize

        row  = (X * arange(M/X)).repeat(diff(self.indptr))
        row  = row.repeat(X*Y).reshape(-1,X,Y)
        row += tile( arange(X).reshape(-1,1), (1,Y) )
        row  = row.reshape(-1) 

        col  = (Y * self.indices).repeat(X*Y).reshape(-1,X,Y)
        col += tile( arange(Y), (X,1) )
        col  = col.reshape(-1)

        data = self.data.reshape(-1)

        if copy:
            data = data.copy()

        from coo import coo_matrix
        return coo_matrix( (data,(row,col)), shape=self.shape )


    def transpose(self):
        
        X,Y = self.blocksize
        M,N = self.shape
        
        if self.nnz == 0:
            return bsr_matrix((N,M),blocksize=(Y,X))

        #use CSR.T to determine a permutation for BSR.T
        from csr import csr_matrix
        data = arange(len(self.indices), dtype=self.indices.dtype)
        proxy = csr_matrix((data,self.indices,self.indptr),shape=(M/X,N/Y))
        proxy = proxy.tocsc()

        data    = self.data.swapaxes(1,2)[proxy.data] #permute data

        indices = proxy.indices
        indptr  = proxy.indptr
       
        return bsr_matrix( (data,indices,indptr), shape=(N,M) )
    
    def tobsr(self,blocksize=None,copy=False):

        if blocksize not in [None, self.blocksize]:
            return self.tocoo(copy=False).tobsr(blocksize=blocksize)
        if copy:
            return self.copy()
        else:
            return self
    
    # these functions are used by the parent class
    # to remove redudancy between bsc_matrix and bsr_matrix
    def _swap(self,x):
        """swap the members of x if this is a column-oriented matrix
        """
        return (x[0],x[1])


from sputils import _isinstance

def isspmatrix_bsr(x):
    return _isinstance(x, bsr_matrix)

