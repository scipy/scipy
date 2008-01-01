"""Compressed Block Sparse Row matrix format"""

__all__ = ['bsr_matrix', 'isspmatrix_bsr']

from warnings import warn

from numpy import zeros, intc, array, asarray, arange, diff, tile, rank, \
        prod, ravel, empty, matrix, asmatrix, empty_like, hstack

import sparsetools
from sparsetools import bsr_matvec, csr_matmat_pass1, bsr_matmat_pass2
from data import _data_matrix
from compressed import _cs_matrix
from base import isspmatrix
from sputils import isshape, getdtype, to_native, isscalarlike, isdense, \
        upcast

class bsr_matrix(_cs_matrix):
    """Block Sparse Row matrix

    This can be instantiated in several ways:
      - bsr_matrix(D, [blocksize=(R,C)])
        - with a dense matrix or rank-2 ndarray D

      - bsr_matrix(S, [blocksize=(R,C)])
        - with another sparse matrix S (equivalent to S.tobsr())

      - bsr_matrix((M, N), [blocksize=(R,C), dtype])
        - to construct an empty matrix with shape (M, N)
        - dtype is optional, defaulting to dtype='d'.

      - bsr_matrix((data, ij), [blocksize=(R,C), shape=(M, N)])
        - where data, ij satisfy:
          - a[ij[0, k], ij[1, k]] = data[k]

      - bsr_matrix((data, indices, indptr), [shape=(M, N)])
        - is the standard BSR representation where:
          the block column indices for row i are stored in
           - indices[ indptr[i]: indices[i+1] ] 
          and their corresponding block values are stored in
           - data[ indptr[i]: indptr[i+1] ]
        - if the shape parameter is not supplied, the matrix dimensions
          are inferred from the index arrays.


    Notes
    =====
        
        - The blocksize (R,C) must evenly divide the shape of 
          the matrix (M,N).  That is, R and C must satisfy the
          relationship M % R = 0 and N % C = 0.
    
        - The Block Compressed Row (BSR) format is very similar to the
          Compressed Sparse Row (CSR) format.  BSR is appropriate for
          sparse matrices with dense sub matrices like the last example
          below.  Such matrices often arise, for instance, in finite
          element discretizations.


    Examples
    ========

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
    def __init__(self, arg1, shape=None, dtype=None, copy=False, blocksize=None):
        _data_matrix.__init__(self)

        
        if isspmatrix(arg1):
            if arg1.format == self.format and copy:
                arg1 = arg1.copy()
            else:
                arg1 = getattr(arg1,'to' + self.format)(blocksize=blocksize)
            self._set_self( arg1 )
            
        elif isinstance(arg1,tuple):
            if isshape(arg1):
                #it's a tuple of matrix dimensions (M,N)
                self.shape  = arg1
                M,N = self.shape
                #process blocksize
                if blocksize is None:
                    blocksize = (1,1)
                else:
                    if not isshape(blocksize):
                        raise ValueError,'invalid blocksize=%s',blocksize
                    blocksize = tuple(blocksize)
                self.data   = zeros( (0,) + blocksize, getdtype(dtype, default=float) )
                self.indices = zeros( 0, dtype=intc )
                
                X,Y = blocksize
                if (M % X) != 0 or (N % Y) != 0:
                    raise ValueError, 'shape must be multiple of blocksize'

                self.indptr  = zeros(M/X + 1, dtype=intc )
            
            elif len(arg1) == 2:
                # (data,(row,col)) format
                self._set_self( coo_matrix(arg1).tobsr(blocksize=blocksize) )

            elif len(arg1) == 3:
                # (data,indices,indptr) format
                (data, indices, indptr) = arg1
                self.indices = array(indices, copy=copy)
                self.indptr  = array(indptr,  copy=copy)
                self.data    = array(data,    copy=copy, \
                        dtype=getdtype(dtype, data))
            else:
                raise ValueError,'unrecognized bsr_matrix constructor usage'
        else:
            #must be dense
            try:
                arg1 = asarray(arg1)
            except:
                raise ValueError, "unrecognized form for" \
                        " %s_matrix constructor" % self.format
            from coo import coo_matrix
            arg1 = self.__class__( coo_matrix(arg1), blocksize=blocksize )
            self._set_self( arg1 )

        if shape is not None:
            self.shape = shape   # spmatrix will check for errors
        else:
            if self.shape is None:
                # shape not already set, try to infer dimensions
                try:
                    M = len(self.indptr) - 1
                    N = self.indices.max() + 1
                except:
                    raise ValueError,'unable to infer matrix dimensions'
                else:
                    R,C = self.blocksize
                    self.shape = (M*R,N*C)

        if self.shape is None:
            if shape is None:
                #infer shape here
                raise ValueError,'need to infer shape'
            else:
                self.shape = shape

        self.check_format()

    def check_format(self, full_check=True):
        """check whether the matrix format is valid

            *Parameters*:
                full_check:
                    True  - rigorous check, O(N) operations : default
                    False - basic check, O(1) operations

        """
        M,N = self.shape
        R,C = self.blocksize

        # index arrays should have integer data types
        if self.indptr.dtype.kind != 'i':
            warn("indptr array has non-integer dtype (%s)" \
                    % self.indptr.dtype.name )
        if self.indices.dtype.kind != 'i':
            warn("indices array has non-integer dtype (%s)" \
                    % self.indices.dtype.name )

        # only support 32-bit ints for now
        self.indptr  = asarray(self.indptr,intc)
        self.indices = asarray(self.indices,intc)
        self.data    = to_native(self.data)

        # check array shapes
        if (rank(self.indices) != 1) or (rank(self.indptr) != 1):
            raise ValueError,"indices, and indptr should be rank 1"
        if rank(self.data) != 3:
            raise ValueError,"data should be rank 3"

        # check index pointer
        if (len(self.indptr) != M/R + 1 ):
            raise ValueError, \
                "index pointer size (%d) should be (%d)" % \
                 (len(self.indptr), M/R + 1)
        if (self.indptr[0] != 0):
            raise ValueError,"index pointer should start with 0"

        # check index and data arrays
        if (len(self.indices) != len(self.data)):
            raise ValueError,"indices and data should have the same size"
        if (self.indptr[-1] > len(self.indices)):
            raise ValueError, \
                  "Last value of index pointer should be less than "\
                  "the size of index and data arrays"

        self.prune()

        if full_check:
            #check format validity (more expensive)
            if self.nnz > 0:
                if self.indices.max() >= N/C:
                    print "max index",self.indices.max()
                    raise ValueError, "column index values must be < %d" % N/C
                if self.indices.min() < 0:
                    raise ValueError, "column index values must be >= 0"
                if diff(self.indptr).min() < 0:
                    raise ValueError,'index pointer values must form a " \
                                        "non-decreasing sequence'

        #if not self.has_sorted_indices():
        #    warn('Indices were not in sorted order. Sorting indices.')
        #    self.sort_indices(check_first=False)

    def _get_blocksize(self):
        return self.data.shape[1:]
    blocksize = property(fget=_get_blocksize)
    
    def getnnz(self):
        R,C = self.blocksize
        return self.indptr[-1] * R * C
    nnz = property(fget=getnnz)
    
    def __repr__(self):
        nnz = self.getnnz()
        format = self.getformat()
        return "<%dx%d sparse matrix of type '%s'\n" \
               "\twith %d stored elements (blocksize = %dx%d) in %s format>" % \
               ( self.shape + (self.dtype.type, nnz) + self.blocksize + \
                 (_formats[format][1],) )


    ##########################
    # NotImplemented methods #
    ##########################

    def getdata(self,ind):
        raise NotImplementedError

    def __getitem__(self,key):
        raise NotImplementedError
    
    def __setitem__(self,key,val):
        raise NotImplementedError

    ######################
    # Arithmetic methods #
    ######################

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

  
    ######################
    # Conversion methods #
    ######################
    
    def tobsr(self,blocksize=None,copy=False):
        if blocksize not in [None, self.blocksize]:
            return self.tocoo(copy=False).tobsr(blocksize=blocksize)
        if copy:
            return self.copy()
        else:
            return self

    def tocsr(self):
        return self.tocoo(copy=False).tocsr()
        #TODO make this more efficient

    def tocsc(self):
        return self.tocoo(copy=False).tocsc()

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
    
    
    ############################################################## 
    # methods that examine or modify the internal data structure #
    ##############################################################
   
    def sum_duplicates(self):
        raise NotImplementedError

    def sort_indices(self):
        """Sort the indices of this matrix *in place*
        """
        if self.has_sorted_indices():
            return

        from csr import csr_matrix

        R,C = self.blocksize
        M,N = self.shape

        if self.nnz == 0:
            return

        #use CSR.sort_indices to determine a permutation for BSR blocks
        data = arange(len(self.indices), dtype=self.indices.dtype)
        proxy = csr_matrix((data,self.indices,self.indptr),shape=(M/R,N/C))
        proxy.sort_indices()

        self.data[:] = self.data[proxy.data]
        self.indices[:] = proxy.indices

        self._has_sorted_indices = True

    def prune(self):
        """ Remove empty space after all non-zero elements.
        """

        R,C = self.blocksize
        M,N = self.shape

        if len(self.indptr) != M/R + 1:
            raise ValueError, "index pointer has invalid length"
        
        bnnz = self.indptr[-1]

        if len(self.indices) < bnnz: 
            raise ValueError, "indices array has too few elements"
        if len(self.data) < bnnz:
            raise ValueError, "data array has too few elements"

        self.data    = self.data[:bnnz]
        self.indices = self.indices[:bnnz]
    
    # utility functions
    def _binopt(self, other, op, in_shape=None, out_shape=None):
        """apply the binary operation fn to two sparse matrices"""
        other = self.__class__(other,blocksize=self.blocksize)

        if in_shape is None:
            in_shape = self.shape
        if out_shape is None:
            out_shape = self.shape
        
        self.sort_indices()
        other.sort_indices()

        # e.g. bsr_plus_bsr, etc.
        fn = getattr(sparsetools, self.format + op + self.format)
        
        R,C = self.blocksize 

        max_bnnz = len(self.data) + len(other.data)
        indptr  = empty_like(self.indptr)
        indices = empty( max_bnnz, dtype=intc )
        data    = empty( R*C*max_bnnz, dtype=upcast(self.dtype,other.dtype) )

        fn(in_shape[0]/R, in_shape[1]/C, R, C, \
                self.indptr,  self.indices,  ravel(self.data),
                other.indptr, other.indices, ravel(other.data),
                indptr,       indices,       data)
        
        actual_bnnz = indptr[-1]
        indices = indices[:actual_bnnz]
        data    = data[:R*C*actual_bnnz]

        if actual_bnnz < max_bnnz/2:
            indices = indices.copy()
            data    = data.copy()

        data = data.reshape(-1,R,C)

        return self.__class__((data, indices, indptr), shape=out_shape)

    # needed by _data_matrix
    def _with_data(self,data,copy=True):
        """Returns a matrix with the same sparsity structure as self,
        but with different data.  By default the structure arrays
        (i.e. .indptr and .indices) are copied.
        """
        if copy:
            return self.__class__((data,self.indices.copy(),self.indptr.copy()), \
                                   shape=self.shape,dtype=data.dtype)
        else:
            return self.__class__((data,self.indices,self.indptr), \
                                   shape=self.shape,dtype=data.dtype)


    
#    # these functions are used by the parent class
#    # to remove redudancy between bsc_matrix and bsr_matrix
#    def _swap(self,x):
#        """swap the members of x if this is a column-oriented matrix
#        """
#        return (x[0],x[1])


from sputils import _isinstance

def isspmatrix_bsr(x):
    return _isinstance(x, bsr_matrix)

