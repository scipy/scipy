""" A sparse matrix in COOrdinate format """

__all__ = ['coo_matrix', 'isspmatrix_coo']

from itertools import izip
from warnings import warn 

from numpy import array, asarray, empty, intc, zeros, bincount

from sparsetools import coo_tocsr, coo_tocsc
from base import spmatrix, isspmatrix
from sputils import upcast, to_native, isshape, getdtype

class coo_matrix(spmatrix):
    """ A sparse matrix in coordinate list format.

    COO matrices are created either as:
        A = coo_matrix( (m, n), [dtype])
    for a zero matrix, or as:
        A = coo_matrix(S)
    where S is a sparse matrix, or as:
        A = coo_matrix(M)
    where M is a dense matrix or rank 2 ndarray, or as:
        A = coo_matrix((obj, ij), [shape])
    where the dimensions are optional.  If supplied, we set (M, N) = shape.
    If not supplied, we infer these from the index arrays:
        ij[0][:] and ij[1][:]

    The arguments 'obj' and 'ij' represent three arrays:
        1. obj[:]    the entries of the matrix, in any order
        2. ij[0][:]  the row indices of the matrix entries
        3. ij[1][:]  the column indices of the matrix entries

    So the following holds:
        A[ij[0][k], ij[1][k] = obj[k]

    Note:
        When converting to CSR or CSC format, duplicate (i,j) entries
        will be summed together.  This facilitates efficient construction
        of finite element matrices and the like.

    """
    def __init__(self, arg1, shape=None, dtype=None, copy=False, dims=None):
        spmatrix.__init__(self)

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
                    M = asarray(arg1)
                except:
                    raise TypeError, "invalid input format"
    
                if len(M.shape) != 2:
                    raise TypeError, "expected rank 2 array or matrix"
                self.shape = M.shape
                self.row,self.col = (M != 0).nonzero()
                self.data  = M[self.row,self.col]

        self._check()

    def _get_dtype(self):
        return self.data.dtype
    def _set_dtype(self,newtype):
        self.data.dtype = newtype
    dtype = property(fget=_get_dtype,fset=_set_dtype)

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
        self.row  = self.row.astype('intc')
        self.col  = self.col.astype('intc')
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

    def tocsc(self):
        from csc import csc_matrix
        if self.nnz == 0:
            return csc_matrix(self.shape, dtype=self.dtype)
        else:
            indptr  = empty(self.shape[1] + 1,dtype=intc)
            indices = empty(self.nnz, dtype=intc)
            data    = empty(self.nnz, dtype=upcast(self.dtype))

            coo_tocsc(self.shape[0], self.shape[1], self.nnz, \
                      self.row, self.col, self.data, \
                      indptr, indices, data)

            return csc_matrix((data, indices, indptr), self.shape)

    def tocsr(self):
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

            return csr_matrix((data, indices, indptr), self.shape)
    

    def tocoo(self, copy=False):
        if copy:
            return self.copy()
        else:
            return self

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



from sputils import _isinstance

def isspmatrix_coo( x ):
    return _isinstance(x, coo_matrix)

