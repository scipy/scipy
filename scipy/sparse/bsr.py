"""Compressed Block Sparse Row matrix format"""

__docformat__ = "restructuredtext en"

__all__ = ['bsr_matrix', 'isspmatrix_bsr']

from warnings import warn

import numpy as np

from data import _data_matrix
from compressed import _cs_matrix
from base import isspmatrix, _formats
from sputils import isshape, getdtype, to_native, upcast
import sparsetools
from sparsetools import bsr_matvec, bsr_matvecs, csr_matmat_pass1, \
                        bsr_matmat_pass2, bsr_transpose, bsr_sort_indices

class bsr_matrix(_cs_matrix):
    """Block Sparse Row matrix

    This can be instantiated in several ways:
        bsr_matrix(D, [blocksize=(R,C)])
            with a dense matrix or rank-2 ndarray D

        bsr_matrix(S, [blocksize=(R,C)])
            with another sparse matrix S (equivalent to S.tobsr())

        bsr_matrix((M, N), [blocksize=(R,C), dtype])
            to construct an empty matrix with shape (M, N)
            dtype is optional, defaulting to dtype='d'.

        bsr_matrix((data, ij), [blocksize=(R,C), shape=(M, N)])
            where ``data`` and ``ij`` satisfy ``a[ij[0, k], ij[1, k]] = data[k]``

        bsr_matrix((data, indices, indptr), [shape=(M, N)])
            is the standard BSR representation where the block column
            indices for row i are stored in ``indices[indptr[i]:indices[i+1]]``
            and their corresponding block values are stored in
            ``data[ indptr[i]: indptr[i+1] ]``.  If the shape parameter is not
            supplied, the matrix dimensions are inferred from the index arrays.

    Notes
    -----

    Summary
        - The Block Compressed Row (BSR) format is very similar to the
          Compressed Sparse Row (CSR) format.  BSR is appropriate for
          sparse matrices with dense sub matrices like the last example
          below.  Block matrices often arise in vector-valued finite
          element discretizations.  In such cases, BSR is considerably
          more efficient than CSR and CSC for many sparse arithmetic
          operations.

    Blocksize
        - The blocksize (R,C) must evenly divide the shape of
          the matrix (M,N).  That is, R and C must satisfy the
          relationship M % R = 0 and N % C = 0.
        - If no blocksize is specified, a simple heuristic is applied
          to determine an appropriate blocksize.



    Examples
    --------

    >>> from scipy.sparse import *
    >>> from scipy import *
    >>> bsr_matrix( (3,4), dtype=int8 ).todense()
    matrix([[0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]], dtype=int8)

    >>> row  = array([0,0,1,2,2,2])
    >>> col  = array([0,2,2,0,1,2])
    >>> data = array([1,2,3,4,5,6])
    >>> bsr_matrix( (data,(row,col)), shape=(3,3) ).todense()
    matrix([[1, 0, 2],
            [0, 0, 3],
            [4, 5, 6]])

    >>> indptr  = array([0,2,3,6])
    >>> indices = array([0,2,2,0,1,2])
    >>> data    = array([1,2,3,4,5,6]).repeat(4).reshape(6,2,2)
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
            if isspmatrix_bsr(arg1) and copy:
                arg1 = arg1.copy()
            else:
                arg1 = arg1.tobsr(blocksize=blocksize)
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
                        raise ValueError('invalid blocksize=%s' % blocksize)
                    blocksize = tuple(blocksize)
                self.data    = np.zeros( (0,) + blocksize, getdtype(dtype, default=float) )
                self.indices = np.zeros( 0, dtype=np.intc )

                R,C = blocksize
                if (M % R) != 0 or (N % C) != 0:
                    raise ValueError, 'shape must be multiple of blocksize'

                self.indptr  = np.zeros(M/R + 1, dtype=np.intc )

            elif len(arg1) == 2:
                # (data,(row,col)) format
                from coo import coo_matrix
                self._set_self( coo_matrix(arg1, dtype=dtype).tobsr(blocksize=blocksize) )

            elif len(arg1) == 3:
                # (data,indices,indptr) format
                (data, indices, indptr) = arg1
                self.indices = np.array(indices, copy=copy)
                self.indptr  = np.array(indptr,  copy=copy)
                self.data    = np.array(data,    copy=copy, dtype=getdtype(dtype, data))
            else:
                raise ValueError('unrecognized bsr_matrix constructor usage')
        else:
            #must be dense
            try:
                arg1 = np.asarray(arg1)
            except:
                raise ValueError("unrecognized form for" \
                        " %s_matrix constructor" % self.format)
            from coo import coo_matrix
            arg1 = coo_matrix(arg1, dtype=dtype).tobsr(blocksize=blocksize)
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
                    raise ValueError('unable to infer matrix dimensions')
                else:
                    R,C = self.blocksize
                    self.shape = (M*R,N*C)

        if self.shape is None:
            if shape is None:
                #TODO infer shape here
                raise ValueError('need to infer shape')
            else:
                self.shape = shape

        if dtype is not None:
            self.data = self.data.astype(dtype)

        self.check_format(full_check=False)

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
        self.indptr  = np.asarray(self.indptr, np.intc)
        self.indices = np.asarray(self.indices, np.intc)
        self.data    = to_native(self.data)

        # check array shapes
        if np.rank(self.indices) != 1 or np.rank(self.indptr) != 1:
            raise ValueError,"indices, and indptr should be rank 1"
        if np.rank(self.data) != 3:
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
                    raise ValueError, "column index values must be < %d" % (N/C)
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


    def diagonal(self):
        """Returns the main diagonal of the matrix
        """
        M,N = self.shape
        R,C = self.blocksize
        y = np.empty(min(M,N), dtype=upcast(self.dtype))
        sparsetools.bsr_diagonal(M/R, N/C, R, C, \
                self.indptr, self.indices, np.ravel(self.data), y)
        return y

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

    def matvec(self, other):
        return self * other

    def matmat(self, other):
        return self * other

    def _mul_vector(self, other):
        M,N = self.shape
        R,C = self.blocksize

        result = np.zeros(self.shape[0], dtype=upcast(self.dtype, other.dtype))

        bsr_matvec(M/R, N/C, R, C, \
            self.indptr, self.indices, self.data.ravel(),
            other, result)

        return result

    def _mul_multivector(self,other):
        R,C = self.blocksize
        M,N = self.shape
        n_vecs = other.shape[1] #number of column vectors

        result = np.zeros((M,n_vecs), dtype=upcast(self.dtype,other.dtype))

        bsr_matvecs(M/R, N/C, n_vecs, R, C, \
                self.indptr, self.indices, self.data.ravel(), \
                other.ravel(), result.ravel())

        return result

    def _mul_sparse_matrix(self, other):
        M, K1 = self.shape
        K2, N = other.shape

        indptr = np.empty_like( self.indptr )

        R,n = self.blocksize

        #convert to this format
        if isspmatrix_bsr(other):
            C = other.blocksize[1]
        else:
            C = 1

        from csr import isspmatrix_csr

        if isspmatrix_csr(other) and n == 1:
            other = other.tobsr(blocksize=(n,C), copy=False) #lightweight conversion
        else:
            other = other.tobsr(blocksize=(n,C))

        csr_matmat_pass1( M/R, N/C, \
                self.indptr,  self.indices, \
                other.indptr, other.indices, \
                indptr)

        bnnz = indptr[-1]
        indices = np.empty(bnnz, dtype=np.intc)
        data    = np.empty(R*C*bnnz, dtype=upcast(self.dtype,other.dtype))

        bsr_matmat_pass2( M/R, N/C, R, C, n, \
                self.indptr,  self.indices,  np.ravel(self.data), \
                other.indptr, other.indices, np.ravel(other.data), \
                indptr,       indices,       data)

        data = data.reshape(-1,R,C)

        #TODO eliminate zeros

        return bsr_matrix((data,indices,indptr),shape=(M,N),blocksize=(R,C))




    ######################
    # Conversion methods #
    ######################

    def tobsr(self,blocksize=None,copy=False):
        if blocksize not in [None, self.blocksize]:
            return self.tocsr().tobsr(blocksize=blocksize)
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
        R,C = self.blocksize

        row  = (R * np.arange(M/R)).repeat(np.diff(self.indptr))
        row  = row.repeat(R*C).reshape(-1,R,C)
        row += np.tile(np.arange(R).reshape(-1,1), (1,C))
        row  = row.reshape(-1)

        col  = (C * self.indices).repeat(R*C).reshape(-1,R,C)
        col += np.tile(np.arange(C), (R,1))
        col  = col.reshape(-1)

        data = self.data.reshape(-1)

        if copy:
            data = data.copy()

        from coo import coo_matrix
        return coo_matrix((data,(row,col)), shape=self.shape)


    def transpose(self):

        R,C = self.blocksize
        M,N = self.shape
        NBLK = self.nnz/(R*C)

        if self.nnz == 0:
            return bsr_matrix((N,M), blocksize=(C,R))

        indptr  = np.empty( N/C + 1,    dtype=self.indptr.dtype)
        indices = np.empty( NBLK,       dtype=self.indices.dtype)
        data    = np.empty( (NBLK,C,R), dtype=self.data.dtype)

        bsr_transpose(M/R, N/C, R, C, \
                      self.indptr, self.indices, self.data.ravel(), \
                      indptr,      indices,      data.ravel())

        return bsr_matrix((data,indices,indptr), shape=(N,M))


    ##############################################################
    # methods that examine or modify the internal data structure #
    ##############################################################

    def eliminate_zeros(self):
        R,C = self.blocksize
        M,N = self.shape

        mask = (self.data != 0).reshape(-1,R*C).sum(axis=1) #nonzero blocks

        nonzero_blocks = mask.nonzero()[0]

        if len(nonzero_blocks) == 0:
            return #nothing to do

        self.data[:len(nonzero_blocks)] = self.data[nonzero_blocks]

        from csr import csr_matrix

        # modifies self.indptr and self.indices *in place*
        proxy = csr_matrix((mask,self.indices,self.indptr),shape=(M/R,N/C))
        proxy.eliminate_zeros()

        self.prune()


    def sum_duplicates(self):
        raise NotImplementedError

    def sort_indices(self):
        """Sort the indices of this matrix *in place*
        """
        if self.has_sorted_indices:
            return

        R,C = self.blocksize
        M,N = self.shape

        bsr_sort_indices(M/R, N/C, R, C, self.indptr, self.indices, self.data.ravel())

        self.has_sorted_indices = True

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

        # ideally we'd take the GCDs of the blocksize dimensions
        # and explode self and other to match
        other = self.__class__(other, blocksize=self.blocksize)

        # e.g. bsr_plus_bsr, etc.
        fn = getattr(sparsetools, self.format + op + self.format)

        R,C = self.blocksize

        max_bnnz = len(self.data) + len(other.data)
        indptr  = np.empty_like(self.indptr)
        indices = np.empty(max_bnnz, dtype=np.intc)
        data    = np.empty(R*C*max_bnnz, dtype=upcast(self.dtype,other.dtype))

        fn(self.shape[0]/R, self.shape[1]/C, R, C,
                self.indptr,  self.indices,  np.ravel(self.data),
                other.indptr, other.indices, np.ravel(other.data),
                indptr,       indices,       data)

        actual_bnnz = indptr[-1]
        indices = indices[:actual_bnnz]
        data    = data[:R*C*actual_bnnz]

        if actual_bnnz < max_bnnz/2:
            indices = indices.copy()
            data    = data.copy()

        data = data.reshape(-1,R,C)

        return self.__class__((data, indices, indptr), shape=self.shape)

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
