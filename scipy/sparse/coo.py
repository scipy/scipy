""" A sparse matrix in COOrdinate or 'triplet' format"""

__all__ = ['coo_matrix', 'isspmatrix_coo']

from itertools import izip
from warnings import warn 

from numpy import array, asarray, empty, intc, zeros, bincount, \
        unique, searchsorted, atleast_2d, lexsort, cumsum, concatenate, \
        empty_like, arange

from sparsetools import coo_tocsr, coo_tocsc
from base import isspmatrix
from data import _data_matrix
from sputils import upcast, to_native, isshape, getdtype
from spfuncs import estimate_blocksize

class coo_matrix(_data_matrix):
    """A sparse matrix in COOrdinate format.
    Also known as the 'ijv' or 'triplet' format.

    This can be instantiated in several ways:
      - coo_matrix(D)
        - with a dense matrix D

      - coo_matrix(S)
        - with another sparse matrix S (equivalent to S.tocoo())

      - coo_matrix((M, N), [dtype])
        - to construct an empty matrix with shape (M, N)
          dtype is optional, defaulting to dtype='d'.

      - coo_matrix((data, ij), [shape=(M, N)])
        - When shape is not specified, it is inferred from the index arrays:
          - ij[0][:] and ij[1][:]

        - The arguments 'data' and 'ij' represent three arrays:
             1. data[:]   the entries of the matrix, in any order
             2. ij[0][:]  the row indices of the matrix entries
             3. ij[1][:]  the column indices of the matrix entries
          So the following holds:
           - A[ij[0][k], ij[1][k] = data[k]

    Notes
    =====
        Advantages of the COO format
        ----------------------------
          - facilitates fast conversion among sparse formats
          - permits duplicate entries (see example)
          - faster conversion to CSR/CSC than LIL
        
        Disadvantages of the COO format
        -------------------------------
          - does not currently support (forces COO->CSR conversion) 
            - arithmetic operations
            - slicing
            - matrix vector products
        
        Usage
        -----
          - COO is a fast format for constructing sparse matrices
          - once a matrix has been constructed, convert to CSR or 
            CSC format for fast arithmetic and matrix vector operations
          - By default when converting to CSR or CSC format, duplicate (i,j) 
            entries will be summed together.  This facilitates efficient 
            construction of finite element matrices and the like. (see example)

    Examples
    ========

    >>> from scipy.sparse import *
    >>> from scipy import *
    >>> coo_matrix( (3,4), dtype=int8 ).todense()
    matrix([[0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]], dtype=int8)

    >>> row  = array([0,3,1,0])
    >>> col  = array([0,3,1,2])
    >>> data = array([4,5,7,9])
    >>> coo_matrix( (data,(row,col)), shape=(4,4) ).todense()
    matrix([[4, 0, 9, 0],
            [0, 7, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 5]])

    >>> # example with duplicates
    >>> row  = array([0,0,1,3,1,0,0])
    >>> col  = array([0,2,1,3,1,0,0])
    >>> data = array([1,1,1,1,1,1,1])
    >>> coo_matrix( (data,(row,col)), shape=(4,4)).todense()
    matrix([[3, 0, 1, 0],
            [0, 2, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 1]])
    
    """

    def __init__(self, arg1, shape=None, dtype=None, copy=False, dims=None):
        _data_matrix.__init__(self)

        if dims is not None:
            warn("dims is deprecated, use shape instead", DeprecationWarning)
            shape=dims

        if isinstance(arg1, tuple):
            if isshape(arg1):
                M, N = arg1
                self.shape = (M,N)
                self.row  = array([], dtype=intc)
                self.col  = array([], dtype=intc)
                self.data = array([], getdtype(dtype, default=float))
            else:
                try:
                    obj, ij = arg1
                except:
                    raise TypeError, "invalid input format"

                try:
                    if len(ij) != 2:
                        raise TypeError
                except TypeError:
                    raise TypeError, "invalid input format"

                self.row = array(ij[0],copy=copy)
                self.col = array(ij[1],copy=copy)
                self.data = array(obj,copy=copy)

                if shape is None:
                    if len(self.row) == 0 or len(self.col) == 0:
                        raise ValueError, "cannot infer dimensions from zero sized index arrays"
                    M = self.row.max() + 1
                    N = self.col.max() + 1
                    self.shape = (M, N)
                else:
                    # Use 2 steps to ensure shape has length 2.
                    M, N = shape
                    self.shape = (M, N)

        elif arg1 is None:
            # Initialize an empty matrix.
            if not isinstance(shape, tuple) or not isintlike(shape[0]):
                raise TypeError, "dimensions not understood"
            warn('coo_matrix(None, shape=(M,N)) is deprecated, ' \
                    'use coo_matrix( (M,N) ) instead', DeprecationWarning)
            self.shape = shape
            self.data = array([],getdtype(dtype, default=float))
            self.row = array([],dtype=intc)
            self.col = array([],dtype=intc)
        else:
            if isspmatrix(arg1):
                if isspmatrix_coo(arg1) and copy:
                    self.row   = arg1.row.copy()
                    self.col   = arg1.col.copy()
                    self.data  = arg1.data.copy()
                    self.shape = arg1.shape
                else:
                    coo = arg1.tocoo()
                    self.row   = coo.row
                    self.col   = coo.col
                    self.data  = coo.data
                    self.shape = coo.shape
            else:
                #dense argument
                try:
                    M = atleast_2d(asarray(arg1))
                except:
                    raise TypeError, "invalid input format"
    
                if len(M.shape) != 2:
                    raise TypeError, "expected rank <= 2 array or matrix"
                self.shape = M.shape
                self.row,self.col = (M != 0).nonzero()
                self.data  = M[self.row,self.col]

        self._check()


    def _check(self):
        """ Checks for consistency and stores the number of non-zeros as
        self.nnz.
        """
        nnz = len(self.data)
        if (nnz != len(self.row)) or (nnz != len(self.col)):
            raise ValueError, "row, column, and data array must all be "\
                  "the same length"

        # index arrays should have integer data types
        if self.row.dtype.kind != 'i':
            warn("row index array has non-integer dtype (%s)  " \
                    % self.row.dtype.name )
        if self.col.dtype.kind != 'i':
            warn("col index array has non-integer dtype (%s) " \
                    % self.col.dtype.name )
       
        # only support 32-bit ints for now
        self.row  = asarray(self.row,dtype=intc)
        self.col  = asarray(self.col,dtype=intc)
        self.data = to_native(self.data)

        if nnz > 0:
            if(self.row.max() >= self.shape[0]):
                raise ValueError, "row index exceedes matrix dimensions"
            if(self.col.max() >= self.shape[1]):
                raise ValueError, "column index exceedes matrix dimensions"
            if(self.row.min() < 0):
                raise ValueError, "negative row index found"
            if(self.col.min() < 0):
                raise ValueError, "negative column index found"

        # some functions pass floats
        self.shape = tuple([int(x) for x in self.shape])
        self.nnz = nnz

    def rowcol(self, num):
        return (self.row[num], self.col[num])

    def getdata(self, num):
        return self.data[num]
    
    def transpose(self,copy=False):
        M,N = self.shape
        return coo_matrix((self.data,(self.col,self.row)),(N,M),copy=copy)

    def toarray(self):
        A = self.tocsr().tocoo(copy=False) #eliminate (i,j) duplicates
        M = zeros(self.shape, dtype=self.dtype)
        M[A.row, A.col] = A.data
        return M

    def tocsc(self,sum_duplicates=True):
        """Return a copy of this matrix in Compressed Sparse Column format

            By default sum_duplicates=True and any duplicate 
            matrix entries are added together.

        """
        from csc import csc_matrix
        if self.nnz == 0:
            return csc_matrix(self.shape, dtype=self.dtype)
        else:
            indptr  = empty(self.shape[1] + 1,dtype=intc)
            indices = empty(self.nnz, dtype=intc)
            data    = empty(self.nnz, dtype=upcast(self.dtype))

            coo_tocsr(self.shape[1], self.shape[0], self.nnz, \
                      self.col, self.row, self.data, \
                      indptr, indices, data)

            A = csc_matrix((data, indices, indptr), self.shape)
            if sum_duplicates:
                A.sum_duplicates()
            return A

    def tocsr(self,sum_duplicates=True):
        """Return a copy of this matrix in Compressed Sparse Row format

            By default sum_duplicates=True and any duplicate 
            matrix entries are added together.

        """
        from csr import csr_matrix
        if self.nnz == 0:
            return csr_matrix(self.shape, dtype=self.dtype)
        else:
            indptr  = empty(self.shape[0] + 1,dtype=intc)
            indices = empty(self.nnz, dtype=intc)
            data    = empty(self.nnz, dtype=upcast(self.dtype))

            coo_tocsr(self.shape[0], self.shape[1], self.nnz, \
                      self.row, self.col, self.data, \
                      indptr, indices, data)

            A = csr_matrix((data, indices, indptr), self.shape)
            if sum_duplicates:
                A.sum_duplicates()
            return A
    

    def tocoo(self, copy=False):
        if copy:
            return self.copy()
        else:
            return self

    def todia(self):
        from dia import dia_matrix

        ks = self.col - self.row  #the diagonal for each nonzero          
        diags = unique(ks)

        if len(diags) > 100:
            #probably undesired, should we do something?
            #should todia() have a maxdiags parameter?
            pass

        #initialize and fill in data array
        data = zeros( (len(diags), self.col.max()+1), dtype=self.dtype)
        data[ searchsorted(diags,ks), self.col ] = self.data

        return dia_matrix((data,diags),shape=self.shape)

    def todok(self):
        from dok import dok_matrix

        dok = dok_matrix((self.shape),dtype=self.dtype)

        try:
            dok.update( izip(izip(self.row,self.col),self.data) ) 
        except AttributeError:
            # the dict() call is for Python 2.3 compatibility
            # ideally dok_matrix would accept an iterator
            dok.update( dict( izip(izip(self.row,self.col),self.data) ) )

        return dok

#    def tobsc(self,blocksize=None):
#        if blocksize in [None, (1,1)]:
#            return self.tocsc().tobsc(blocksize)
#        else:
#            return self.transpose().tobsr().transpose()

    def tobsr(self,blocksize=None):
        from bsr import bsr_matrix

        if self.nnz == 0:
            return bsr_matrix(self.shape,blocksize=blocksize,dtype=self.dtype)

        if blocksize is None:
            blocksize = estimate_blocksize(self)
        elif blocksize in (1,1):
            return self.tocsr().tobsr(blocksize)

        M,N = self.shape
        X,Y = blocksize
    
        if (M % X) != 0 or (N % Y) != 0:
            raise ValueError, 'shape must be multiple of blocksize'
    
        i_block,i_sub = divmod(self.row, X)
        j_block,j_sub = divmod(self.col, Y)
    
        perm = lexsort( keys=[j_block,i_block] )
    
        i_block = i_block[perm]
        j_block = j_block[perm]
    
        mask = (i_block[1:] != i_block[:-1]) + (j_block[1:] != j_block[:-1])
        mask = concatenate((array([True]),mask))
    
        #map self.data[n] -> data[map[n],i_sub[n],j_sub[n]]
        map = cumsum(mask)
        num_blocks = map[-1]
        map -= 1
        
        iperm = empty_like(perm) #inverse permutation
        iperm[perm] = arange(len(perm))
        
        data = zeros( (num_blocks,X,Y), dtype=self.dtype )
        data[map[iperm],i_sub,j_sub] = self.data
    
        row = i_block[mask]
        col = j_block[mask]
    
        # now row,col,data form BOO format 
    
        temp = cumsum(bincount(row))
        indptr = zeros( M/X + 1, dtype=intc )
        indptr[1:len(temp)+1] = temp
        indptr[len(temp)+1:] = temp[-1]
       
        return bsr_matrix((data,col,indptr),shape=self.shape)

    # needed by _data_matrix
    def _with_data(self,data,copy=True):
        """Returns a matrix with the same sparsity structure as self,
        but with different data.  By default the index arrays
        (i.e. .row and .col) are copied.
        """
        if copy:
            return coo_matrix( (data, (self.row.copy(), self.col.copy()) ), \
                                   shape=self.shape, dtype=data.dtype)
        else:
            return coo_matrix( (data, (self.row, self.col) ), \
                                   shape=self.shape, dtype=data.dtype)


from sputils import _isinstance

def isspmatrix_coo( x ):
    return _isinstance(x, coo_matrix)

