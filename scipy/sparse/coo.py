""" A sparse matrix in COOrdinate or 'triplet' format"""

__docformat__ = "restructuredtext en"

__all__ = ['coo_matrix', 'isspmatrix_coo']

from warnings import warn

import numpy as np

from sparsetools import coo_tocsr, coo_todense, coo_matvec
from base import isspmatrix
from data import _data_matrix
from sputils import upcast, upcast_char, to_native, isshape, getdtype, isintlike

class coo_matrix(_data_matrix):
    """
    A sparse matrix in COOrdinate format.

    Also known as the 'ijv' or 'triplet' format.

    This can be instantiated in several ways:
        coo_matrix(D)
            with a dense matrix D

        coo_matrix(S)
            with another sparse matrix S (equivalent to S.tocoo())

        coo_matrix((M, N), [dtype])
            to construct an empty matrix with shape (M, N)
            dtype is optional, defaulting to dtype='d'.

        coo_matrix((data, (i, j)), [shape=(M, N)])
            to construct from three arrays:
                1. data[:]   the entries of the matrix, in any order
                2. i[:]      the row indices of the matrix entries
                3. j[:]      the column indices of the matrix entries

            Where ``A[i[k], j[k]] = data[k]``.  When shape is not
            specified, it is inferred from the index arrays

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
        COO format data array of the matrix
    row
        COO format row index array of the matrix
    col
        COO format column index array of the matrix

    Notes
    -----

    Sparse matrices can be used in arithmetic operations: they support
    addition, subtraction, multiplication, division, and matrix power.

    Advantages of the COO format
        - facilitates fast conversion among sparse formats
        - permits duplicate entries (see example)
        - very fast conversion to and from CSR/CSC formats

    Disadvantages of the COO format
        - does not directly support:
            + arithmetic operations
            + slicing

    Intended Usage
        - COO is a fast format for constructing sparse matrices
        - Once a matrix has been constructed, convert to CSR or
          CSC format for fast arithmetic and matrix vector operations
        - By default when converting to CSR or CSC format, duplicate (i,j)
          entries will be summed together.  This facilitates efficient
          construction of finite element matrices and the like. (see example)

    Examples
    --------
    >>> from scipy.sparse import coo_matrix
    >>> coo_matrix((3,4), dtype=np.int8).todense()
    matrix([[0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]], dtype=int8)

    >>> row  = np.array([0,3,1,0])
    >>> col  = np.array([0,3,1,2])
    >>> data = np.array([4,5,7,9])
    >>> coo_matrix((data,(row,col)), shape=(4,4)).todense()
    matrix([[4, 0, 9, 0],
            [0, 7, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 5]])

    >>> # example with duplicates
    >>> row  = np.array([0,0,1,3,1,0,0])
    >>> col  = np.array([0,2,1,3,1,0,0])
    >>> data = np.array([1,1,1,1,1,1,1])
    >>> coo_matrix((data, (row,col)), shape=(4,4)).todense()
    matrix([[3, 0, 1, 0],
            [0, 2, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 1]])

    """
    def __init__(self, arg1, shape=None, dtype=None, copy=False):
        _data_matrix.__init__(self)

        if isinstance(arg1, tuple):
            if isshape(arg1):
                M, N = arg1
                self.shape = (M,N)
                self.row  = np.array([], dtype=np.intc)
                self.col  = np.array([], dtype=np.intc)
                self.data = np.array([], getdtype(dtype, default=float))
            else:
                try:
                    obj, ij = arg1
                except:
                    raise TypeError('invalid input format')

                try:
                    if len(ij) != 2:
                        raise TypeError
                except TypeError:
                    raise TypeError('invalid input format')

                self.row  = np.array(ij[0], copy=copy, dtype=np.intc)
                self.col  = np.array(ij[1], copy=copy, dtype=np.intc)
                self.data = np.array(  obj, copy=copy)

                if shape is None:
                    if len(self.row) == 0 or len(self.col) == 0:
                        raise ValueError('cannot infer dimensions from zero sized index arrays')
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
                raise TypeError('dimensions not understood')
            warn('coo_matrix(None, shape=(M,N)) is deprecated, ' \
                    'use coo_matrix( (M,N) ) instead', DeprecationWarning)
            self.shape = shape
            self.data = np.array([], getdtype(dtype, default=float))
            self.row  = np.array([], dtype=np.intc)
            self.col  = np.array([], dtype=np.intc)
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
                    M = np.atleast_2d(np.asarray(arg1))
                except:
                    raise TypeError('invalid input format')

                if np.rank(M) != 2:
                    raise TypeError('expected rank <= 2 array or matrix')

                self.shape = M.shape
                self.row, self.col = M.nonzero()
                self.data  = M[self.row, self.col]

        if dtype is not None:
            self.data = self.data.astype(dtype)

        self._check()

    def getnnz(self):
        nnz = len(self.data)
        if nnz != len(self.row) or nnz != len(self.col):
            raise ValueError('row, column, and data array must all be the same length')

        if np.rank(self.data) != 1 or np.rank(self.row) != 1 or np.rank(self.col) != 1:
            raise ValueError('row, column, and data arrays must have rank 1')

        return nnz
    nnz = property(fget=getnnz)

    def _check(self):
        """ Checks data structure for consistency """
        nnz = self.nnz

        # index arrays should have integer data types
        if self.row.dtype.kind != 'i':
            warn("row index array has non-integer dtype (%s)  " \
                    % self.row.dtype.name )
        if self.col.dtype.kind != 'i':
            warn("col index array has non-integer dtype (%s) " \
                    % self.col.dtype.name )

        # only support 32-bit ints for now
        self.row  = np.asarray(self.row, dtype=np.intc)
        self.col  = np.asarray(self.col, dtype=np.intc)
        self.data = to_native(self.data)

        if nnz > 0:
            if self.row.max() >= self.shape[0]:
                raise ValueError('row index exceedes matrix dimensions')
            if self.col.max() >= self.shape[1]:
                raise ValueError('column index exceedes matrix dimensions')
            if self.row.min() < 0:
                raise ValueError('negative row index found')
            if self.col.min() < 0:
                raise ValueError('negative column index found')


    def transpose(self, copy=False):
        M,N = self.shape
        return coo_matrix((self.data, (self.col, self.row)), shape=(N,M), copy=copy)

    def toarray(self, order=None, out=None):
        """See the docstring for `spmatrix.toarray`."""
        B = self._process_toarray_args(order, out)
        fortran = int(B.flags.f_contiguous)
        if not fortran and not B.flags.c_contiguous:
            raise ValueError("Output array must be C or F contiguous")
        M,N = self.shape
        coo_todense(M, N, self.nnz, self.row, self.col, self.data,
                    B.ravel('A'), fortran)
        return B

    def tocsc(self):
        """Return a copy of this matrix in Compressed Sparse Column format

        Duplicate entries will be summed together.

        Examples
        --------
        >>> from numpy import array
        >>> from scipy.sparse import coo_matrix
        >>> row  = array([0,0,1,3,1,0,0])
        >>> col  = array([0,2,1,3,1,0,0])
        >>> data = array([1,1,1,1,1,1,1])
        >>> A = coo_matrix( (data,(row,col)), shape=(4,4)).tocsc()
        >>> A.todense()
        matrix([[3, 0, 1, 0],
                [0, 2, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 1]])

        """
        from csc import csc_matrix
        if self.nnz == 0:
            return csc_matrix(self.shape, dtype=self.dtype)
        else:
            M,N = self.shape
            indptr  = np.empty(N + 1,    dtype=np.intc)
            indices = np.empty(self.nnz, dtype=np.intc)
            data    = np.empty(self.nnz, dtype=upcast(self.dtype))

            coo_tocsr(N, M, self.nnz, \
                      self.col, self.row, self.data, \
                      indptr, indices, data)

            A = csc_matrix((data, indices, indptr), shape=self.shape)
            A.sum_duplicates()

            return A

    def tocsr(self):
        """Return a copy of this matrix in Compressed Sparse Row format

        Duplicate entries will be summed together.

        Examples
        --------
        >>> from numpy import array
        >>> from scipy.sparse import coo_matrix
        >>> row  = array([0,0,1,3,1,0,0])
        >>> col  = array([0,2,1,3,1,0,0])
        >>> data = array([1,1,1,1,1,1,1])
        >>> A = coo_matrix( (data,(row,col)), shape=(4,4)).tocsr()
        >>> A.todense()
        matrix([[3, 0, 1, 0],
                [0, 2, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 1]])

        """
        from csr import csr_matrix
        if self.nnz == 0:
            return csr_matrix(self.shape, dtype=self.dtype)
        else:
            M,N = self.shape
            indptr  = np.empty(M + 1,    dtype=np.intc)
            indices = np.empty(self.nnz, dtype=np.intc)
            data    = np.empty(self.nnz, dtype=upcast(self.dtype))

            coo_tocsr(M, N, self.nnz, \
                      self.row, self.col, self.data, \
                      indptr, indices, data)

            A = csr_matrix((data, indices, indptr), shape=self.shape)
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
        diags = np.unique(ks)

        if len(diags) > 100:
            #probably undesired, should we do something?
            #should todia() have a maxdiags parameter?
            pass

        #initialize and fill in data array
        data = np.zeros( (len(diags), self.col.max()+1), dtype=self.dtype)
        data[ np.searchsorted(diags,ks), self.col ] = self.data

        return dia_matrix((data,diags), shape=self.shape)

    def todok(self):
        from itertools import izip
        from dok import dok_matrix

        dok = dok_matrix((self.shape), dtype=self.dtype)

        dok.update( izip(izip(self.row,self.col),self.data) )

        return dok


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

    ###########################
    # Multiplication handlers #
    ###########################

    def _mul_vector(self, other):
        #output array
        result = np.zeros( self.shape[0], dtype=upcast_char(self.dtype.char,
                                                            other.dtype.char) )
        coo_matvec(self.nnz, self.row, self.col, self.data, other, result)
        return result

    def _mul_multivector(self, other):
        return np.hstack( [ self._mul_vector(col).reshape(-1,1) for col in other.T ] )


def isspmatrix_coo( x ):
    return isinstance(x, coo_matrix)
