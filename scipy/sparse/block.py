"""base class for block sparse formats"""

from warnings import warn

from numpy import zeros, intc, array, asarray, arange, diff, tile, rank, \
        prod, ravel, empty, empty_like

from data import _data_matrix
from base import isspmatrix, _formats
from sputils import isshape, getdtype, to_native, isscalarlike, isdense, \
        upcast
import sparsetools

#TODO refactor w/ compressed.py

class _block_matrix(_data_matrix):
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

                self.indptr  = zeros(self._swap((M/X,N/Y))[0] + 1, dtype=intc )

            elif len(arg1) == 3:
                # data,indices,indptr format
                (data, indices, indptr) = arg1
                self.indices = array(indices, copy=copy)
                self.indptr  = array(indptr,  copy=copy)
                self.data    = array(data,    copy=copy, \
                        dtype=getdtype(dtype, data))
            else:
                raise ValueError, "unrecognized form for" \
                        " %s_matrix constructor" % self.format
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
                    major_dim = len(self.indptr) - 1
                    minor_dim = self.indices.max() + 1
                except:
                    raise ValueError,'unable to infer matrix dimensions'
                else:
                    M,N = self._swap((major_dim,minor_dim))
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

        #use _swap to determine proper bounds
        major_name,minor_name = self._swap(('row','column'))
        major_dim,minor_dim = self._swap(self.shape)
        major_blk,minor_blk = self._swap(self.blocksize)

        # index arrays should have integer data types
        if self.indptr.dtype.kind != 'i':
            warn("indptr array has non-integer dtype (%s)" \
                    % self.indptr.dtype.name )
        if self.indices.dtype.kind != 'i':
            warn("indices array has non-integer dtype (%s)" \
                    % self.indices.dtype.name )

        # only support 32-bit ints for now
        if self.indptr.dtype != intc:
            self.indptr  = self.indptr.astype(intc)
        if self.indices.dtype != intc:
            self.indices = self.indices.astype(intc)
        self.data = to_native(self.data)

        # check array shapes
        if (rank(self.indices) != 1) or (rank(self.indptr) != 1):
            raise ValueError,"indices, and indptr should be rank 1"
        if rank(self.data) != 3:
            raise ValueError,"data should be rank 3"

        # check index pointer
        if (len(self.indptr) != major_dim/major_blk + 1 ):
            raise ValueError, \
                "index pointer size (%d) should be (%d)" % \
                 (len(self.indptr), major_dim/major_blk + 1)
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
                if self.indices.max() >= minor_dim/minor_blk:
                    print "max index",self.indices.max()
                    raise ValueError, "%s index values must be < %d" % \
                            (minor_name,minor_dim)
                if self.indices.min() < 0:
                    raise ValueError, "%s index values must be >= 0" % \
                            minor_name
                if diff(self.indptr).min() < 0:
                    raise ValueError,'index pointer values must form a " \
                                        "non-decreasing sequence'

        if not self.has_sorted_indices():
            warn('Indices were not in sorted order. Sorting indices.')
            self.sort_indices(check_first=False)

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


    def __add__(self,other):
        # First check if argument is a scalar
        if isscalarlike(other):
            # Now we would add this scalar to every element.
            raise NotImplementedError, 'adding a scalar to a CSC or CSR ' \
                  'matrix is not supported'
        elif isspmatrix(other):
            if (other.shape != self.shape):
                raise ValueError, "inconsistent shapes"
            
            return self._binopt(other,'_plus_')
        elif isdense(other):
            # Convert this matrix to a dense matrix and add them
            return self.todense() + other
        else:
            raise NotImplementedError

    def __radd__(self,other):
        return self.__add__(other)

    def __sub__(self,other):
        # First check if argument is a scalar
        if isscalarlike(other):
            # Now we would add this scalar to every element.
            raise NotImplementedError, 'adding a scalar to a sparse ' \
                  'matrix is not supported'
        elif isspmatrix(other):
            if (other.shape != self.shape):
                raise ValueError, "inconsistent shapes"

            return self._binopt(other,'_minus_')
        elif isdense(other):
            # Convert this matrix to a dense matrix and subtract them
            return self.todense() - other
        else:
            raise NotImplementedError

    def __rsub__(self,other):  # other - self
        #note: this can't be replaced by other + (-self) for unsigned types
        if isscalarlike(other):
            # Now we would add this scalar to every element.
            raise NotImplementedError, 'adding a scalar to a sparse ' \
                  'matrix is not supported'
        elif isdense(other):
            # Convert this matrix to a dense matrix and subtract them
            return other - self.todense()
        else:
            raise NotImplementedError


    def __mul__(self, other): # self * other
        """ Scalar, vector, or matrix multiplication
        """
        if isscalarlike(other):
            return self._with_data(self.data * other)
        else:
            return self.dot(other)


    def __rmul__(self, other): # other * self
        if isscalarlike(other):
            return self.__mul__(other)
        else:
            # Don't use asarray unless we have to
            try:
                tr = other.transpose()
            except AttributeError:
                tr = asarray(other).transpose()
            return self.transpose().dot(tr).transpose()


    def __truediv__(self,other):
        if isscalarlike(other):
            return self * (1./other)
        elif isspmatrix(other):
            if (other.shape != self.shape):
                raise ValueError, "inconsistent shapes"
            
            return self._binopt(other,'_eldiv_')
        else:
            raise NotImplementedError

    
    def multiply(self, other):
        """Point-wise multiplication by another matrix
        """
        if (other.shape != self.shape):
            raise ValueError, "inconsistent shapes"

        if isdense(other):
            return multiply(self.todense(),other)
        else:
            other = self.__class__(other)
            return self._binopt(other,'_elmul_')








    def _set_self(self, other, copy=False):
        """take the member variables of other and assign them to self"""

        if copy:
            other = other.copy()

        self.data      = other.data
        self.indices   = other.indices
        self.indptr    = other.indptr
        self.shape     = other.shape

   
    #conversion methods
    def toarray(self):
        A = self.tocoo(copy=False)
        M = zeros(self.shape, dtype=self.dtype)
        M[A.row, A.col] = A.data
        return M

    def todia(self):
        return self.tocoo(copy=False).todia()
    
    def todok(self):
        return self.tocoo(copy=False).todok()

    def tocsr(self):
        return self.tocoo(copy=False).tocsr()
        #TODO make this more efficient

    def tocsc(self):
        return self.tocoo(copy=False).tocsc()
        #TODO make this more efficient




    # methods that modify the internal data structure
    def has_sorted_indices(self):
        """Determine whether the matrix has sorted indices
        """
        fn = sparsetools.csr_has_sorted_indices
        return fn( len(self.indptr) - 1, self.indptr, self.indices)

    def sorted_indices(self):
        """Return a copy of this matrix with sorted indices
        """
        A = self.copy()
        A.sort_indices()
        return A

        # an alternative that has linear complexity is the following
        # typically the previous option is faster
        #return self.toother().toother()

    def sort_indices(self, check_first=True):
        """Sort the indices of this matrix *in place*
        """
        from csr import csr_matrix

        X,Y = self.blocksize
        M,N = self.shape

        if self.nnz == 0:
            return

        #use CSR.sort_indices to determine a permutation for BSR<->BSC
        major,minor = self._swap((M/X,N/Y))

        data = arange(len(self.indices), dtype=self.indices.dtype)
        proxy = csr_matrix((data,self.indices,self.indptr),shape=(major,minor))
        proxy.sort_indices()

        self.data[:] = self.data[proxy.data]
        self.indices[:] = proxy.indices

    def prune(self):
        """ Remove empty space after all non-zero elements.
        """
        major_dim = self._swap(self.shape)[0]
        major_blk = self._swap(self.blocksize)[0]

        if len(self.indptr) != major_dim/major_blk + 1:
            raise ValueError, "index pointer has invalid length"
        if len(self.indices) < self.nnz / prod(self.blocksize): 
            raise ValueError, "indices has too few elements"
        if self.data.size < self.nnz:
            raise ValueError, "data array has too few elements"
        
        nnz = self.indptr[-1]

        self.data    = self.data[:nnz]
        self.indices = self.indices[:nnz]

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






    
# test with:
# A = arange(4*6).reshape(4,6) % 5
# A[0:2,2:4] = 0
# A[0:2,:]   = 0




