"""Base class for sparse matrix formats using compressed storage
"""

__all__ = []

from warnings import warn

import numpy as np

from base import spmatrix, isspmatrix, SparseEfficiencyWarning
from data import _data_matrix
import sparsetools
from sputils import upcast, to_native, isdense, isshape, getdtype, \
        isscalarlike, isintlike


class _cs_matrix(_data_matrix):
    """base matrix class for compressed row and column oriented matrices"""

    def __init__(self, arg1, shape=None, dtype=None, copy=False, dims=None, nzmax=None):
        _data_matrix.__init__(self)

        if dims is not None:
            warn("dims= is deprecated, use shape= instead", DeprecationWarning)
            shape=dims

        if nzmax is not None:
            warn("nzmax= is deprecated", DeprecationWarning)


        if isspmatrix(arg1):
            if arg1.format == self.format and copy:
                arg1 = arg1.copy()
            else:
                arg1 = arg1.asformat(self.format)
            self._set_self( arg1 )

        elif isinstance(arg1, tuple):
            if isshape(arg1):
                # It's a tuple of matrix dimensions (M, N)
                # create empty matrix
                self.shape = arg1   #spmatrix checks for errors here
                M, N = self.shape
                self.data    = np.zeros(0, getdtype(dtype, default=float))
                self.indices = np.zeros(0, np.intc)
                self.indptr  = np.zeros(self._swap((M,N))[0] + 1, dtype=np.intc)
            else:
                if len(arg1) == 2:
                    # (data, ij) format
                    from coo import coo_matrix
                    other = self.__class__( coo_matrix(arg1, shape=shape) )
                    self._set_self( other )
                elif len(arg1) == 3:
                    # (data, indices, indptr) format
                    (data, indices, indptr) = arg1
                    self.indices = np.array(indices, copy=copy)
                    self.indptr  = np.array(indptr, copy=copy)
                    self.data    = np.array(data, copy=copy, dtype=getdtype(dtype, data))
                else:
                    raise ValueError, "unrecognized %s_matrix constructor usage" %\
                            self.format

        else:
            #must be dense
            try:
                arg1 = np.asarray(arg1)
            except:
                raise ValueError, "unrecognized %s_matrix constructor usage" % \
                        self.format
            from coo import coo_matrix
            self._set_self( self.__class__(coo_matrix(arg1, dtype=dtype)) )

        # Read matrix dimensions given, if any
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
                    self.shape = self._swap((major_dim,minor_dim))

        if dtype is not None:
            self.data = self.data.astype(dtype)

        self.check_format(full_check=False)

    def getnnz(self):
        return self.indptr[-1]
    nnz = property(fget=getnnz)


    def _set_self(self, other, copy=False):
        """take the member variables of other and assign them to self"""

        if copy:
            other = other.copy()

        self.data    = other.data
        self.indices = other.indices
        self.indptr  = other.indptr
        self.shape   = other.shape

    def check_format(self, full_check=True):
        """check whether the matrix format is valid

        Parameters
        ==========

            - full_check : {bool}
                - True  - rigorous check, O(N) operations : default
                - False - basic check, O(1) operations

        """
        #use _swap to determine proper bounds
        major_name,minor_name = self._swap(('row','column'))
        major_dim,minor_dim = self._swap(self.shape)

        # index arrays should have integer data types
        if self.indptr.dtype.kind != 'i':
            warn("indptr array has non-integer dtype (%s)" \
                    % self.indptr.dtype.name )
        if self.indices.dtype.kind != 'i':
            warn("indices array has non-integer dtype (%s)" \
                    % self.indices.dtype.name )

        # only support 32-bit ints for now
        self.indptr  = np.asarray(self.indptr,  dtype=np.intc)
        self.indices = np.asarray(self.indices, dtype=np.intc)
        self.data    = to_native(self.data)

        # check array shapes
        if np.rank(self.data) != 1 or np.rank(self.indices) != 1 or np.rank(self.indptr) != 1:
            raise ValueError('data, indices, and indptr should be rank 1')

        # check index pointer
        if (len(self.indptr) != major_dim + 1 ):
            raise ValueError, \
                "index pointer size (%d) should be (%d)" % \
                 (len(self.indptr), major_dim + 1)
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
                if self.indices.max() >= minor_dim:
                    raise ValueError, "%s index values must be < %d" % \
                            (minor_name,minor_dim)
                if self.indices.min() < 0:
                    raise ValueError, "%s index values must be >= 0" % \
                            minor_name
                if np.diff(self.indptr).min() < 0:
                    raise ValueError,'index pointer values must form a " \
                                        "non-decreasing sequence'

        #if not self.has_sorted_indices():
        #    warn('Indices were not in sorted order.  Sorting indices.')
        #    self.sort_indices()
        #    assert(self.has_sorted_indices())
        #TODO check for duplicates?


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


    def __truediv__(self,other):
        if isscalarlike(other):
            return self * (1./other)

        elif isspmatrix(other):
            if other.shape != self.shape:
                raise ValueError('inconsistent shapes')

            return self._binopt(other,'_eldiv_')

        else:
            raise NotImplementedError


    def multiply(self, other):
        """Point-wise multiplication by another matrix
        """
        if other.shape != self.shape:
            raise ValueError('inconsistent shapes')

        if isdense(other):
            return np.multiply(self.todense(),other)
        else:
            other = self.__class__(other)
            return self._binopt(other,'_elmul_')


    ###########################
    # Multiplication handlers #
    ###########################

    def _mul_vector(self, other):
        M,N = self.shape

        #output array
        result = np.zeros( self.shape[0], dtype=upcast(self.dtype,other.dtype) )

        # csr_matvec or csc_matvec
        fn = getattr(sparsetools,self.format + '_matvec')
        fn(M, N, self.indptr, self.indices, self.data, other, result)

        return result


    def _mul_multivector(self, other):
        M,N = self.shape
        n_vecs = other.shape[1] #number of column vectors

        result = np.zeros( (M,n_vecs), dtype=upcast(self.dtype,other.dtype) )

        # csr_matvecs or csc_matvecs
        fn = getattr(sparsetools,self.format + '_matvecs')
        fn(M, N, n_vecs, self.indptr, self.indices, self.data, other.ravel(), result.ravel())

        return result


    def _mul_sparse_matrix(self, other):
        M, K1 = self.shape
        K2, N = other.shape

        major_axis = self._swap((M,N))[0]
        indptr = np.empty(major_axis + 1, dtype=np.intc)

        other = self.__class__(other) #convert to this format
        fn = getattr(sparsetools, self.format + '_matmat_pass1')
        fn( M, N, self.indptr, self.indices, \
                  other.indptr, other.indices, \
                  indptr)

        nnz = indptr[-1]
        indices = np.empty(nnz, dtype=np.intc)
        data    = np.empty(nnz, dtype=upcast(self.dtype,other.dtype))

        fn = getattr(sparsetools, self.format + '_matmat_pass2')
        fn( M, N, self.indptr, self.indices, self.data, \
                  other.indptr, other.indices, other.data, \
                  indptr, indices, data)

        return self.__class__((data,indices,indptr),shape=(M,N))


    @np.deprecate
    def getdata(self, ind):
        return self.data[ind]

    def diagonal(self):
        """Returns the main diagonal of the matrix
        """
        #TODO support k-th diagonal
        fn = getattr(sparsetools, self.format + "_diagonal")
        y = np.empty( min(self.shape), dtype=upcast(self.dtype) )
        fn(self.shape[0], self.shape[1], self.indptr, self.indices, self.data, y)
        return y

    def sum(self, axis=None):
        """Sum the matrix over the given axis.  If the axis is None, sum
        over both rows and columns, returning a scalar.
        """
        # The spmatrix base class already does axis=0 and axis=1 efficiently
        # so we only do the case axis=None here
        if axis is None:
            return self.data.sum()
        else:
            return spmatrix.sum(self,axis)
            raise ValueError, "axis out of bounds"

    #######################
    # Getting and Setting #
    #######################

    def __getitem__(self, key):
        if isinstance(key, tuple):
            row = key[0]
            col = key[1]

            #TODO implement CSR[ [1,2,3], X ] with sparse matmat
            #TODO make use of sorted indices

            if isintlike(row) and isintlike(col):
                return self._get_single_element(row,col)
            else:
                major,minor = self._swap((row,col))
                if isintlike(major) and isinstance(minor,slice):
                    minor_shape = self._swap(self.shape)[1]
                    start, stop, stride = minor.indices(minor_shape)
                    out_shape   = self._swap( (1, stop-start) )
                    return self._get_slice( major, start, stop, stride, out_shape)
                elif isinstance( row, slice) or isinstance(col, slice):
                    return self._get_submatrix( row, col )
                else:
                    raise NotImplementedError

        elif isintlike(key):
            return self[key, :]
        else:
            raise IndexError, "invalid index"


    def _get_single_element(self,row,col):
        M, N = self.shape
        if (row < 0):
            row += M
        if (col < 0):
            col += N
        if not (0<=row<M) or not (0<=col<N):
            raise IndexError("index out of bounds")

        major_index, minor_index = self._swap((row,col))

        start = self.indptr[major_index]
        end   = self.indptr[major_index+1]
        indxs = np.where(minor_index == self.indices[start:end])[0]

        num_matches = len(indxs)

        if num_matches == 0:
            # entry does not appear in the matrix
            return 0
        elif num_matches == 1:
            return self.data[start:end][indxs[0]]
        else:
            raise ValueError('nonzero entry (%d,%d) occurs more than once' % (row,col))

    def _get_slice(self, i, start, stop, stride, shape):
        """Returns a copy of the elements
            [i, start:stop:string] for row-oriented matrices
            [start:stop:string, i] for column-oriented matrices
        """
        if stride != 1:
            raise ValueError, "slicing with step != 1 not supported"
        if stop <= start:
            raise ValueError, "slice width must be >= 1"

        #TODO make [i,:] faster
        #TODO implement [i,x:y:z]

        indices = []

        for ind in xrange(self.indptr[i], self.indptr[i+1]):
            if self.indices[ind] >= start and self.indices[ind] < stop:
                indices.append(ind)

        index  = self.indices[indices] - start
        data   = self.data[indices]
        indptr = np.array([0, len(indices)])
        return self.__class__((data, index, indptr), shape=shape, \
                              dtype=self.dtype)

    def _get_submatrix( self, slice0, slice1 ):
        """Return a submatrix of this matrix (new matrix is created)."""

        slice0, slice1 = self._swap((slice0,slice1))
        shape0, shape1 = self._swap(self.shape)
        def _process_slice( sl, num ):
            if isinstance( sl, slice ):
                i0, i1 = sl.start, sl.stop
                if i0 is None:
                    i0 = 0
                elif i0 < 0:
                    i0 = num + i0

                if i1 is None:
                    i1 = num
                elif i1 < 0:
                    i1 = num + i1

                return i0, i1

            elif np.isscalar( sl ):
                if sl < 0:
                    sl += num

                return sl, sl + 1

            else:
                return sl[0], sl[1]

        def _in_bounds( i0, i1, num ):
            if not (0<=i0<num) or not (0<i1<=num) or not (i0<i1):
                raise IndexError,\
                      "index out of bounds: 0<=%d<%d, 0<=%d<%d, %d<%d" %\
                      (i0, num, i1, num, i0, i1)

        i0, i1 = _process_slice( slice0, shape0 )
        j0, j1 = _process_slice( slice1, shape1 )
        _in_bounds( i0, i1, shape0 )
        _in_bounds( j0, j1, shape1 )

        aux = sparsetools.get_csr_submatrix( shape0, shape1,
                                             self.indptr, self.indices,
                                             self.data,
                                             i0, i1, j0, j1 )

        data, indices, indptr = aux[2], aux[1], aux[0]
        shape = self._swap( (i1 - i0, j1 - j0) )

        return self.__class__( (data,indices,indptr), shape=shape )


    def __setitem__(self, key, val):
        if isinstance(key, tuple):
            row,col = key
            if not (isscalarlike(row) and isscalarlike(col)):
                raise NotImplementedError("Fancy indexing in assignment not "
                                          "supported for csr matrices.")
            M, N = self.shape
            if (row < 0):
                row += M
            if (col < 0):
                col += N
            if not (0<=row<M) or not (0<=col<N):
                raise IndexError, "index out of bounds"

            major_index, minor_index = self._swap((row,col))

            start = self.indptr[major_index]
            end   = self.indptr[major_index+1]
            indxs = np.where(minor_index == self.indices[start:end])[0]

            num_matches = len(indxs)
    
            
            if not np.isscalar(val):
                raise ValueError('setting an array element with a sequence')

            val = self.dtype.type(val) 

            if num_matches == 0:
                #entry not already present
                warn('changing the sparsity structure of a %s_matrix is expensive. ' \
                        'lil_matrix is more efficient.' % self.format, \
                        SparseEfficiencyWarning)

                if self.has_sorted_indices:
                    # preserve sorted order
                    newindx = start + self.indices[start:end].searchsorted(minor_index)
                else:
                    newindx = start

                val         = np.array([val],         dtype=self.data.dtype)
                minor_index = np.array([minor_index], dtype=self.indices.dtype)

                self.data    = np.concatenate((self.data[:newindx],    val,         self.data[newindx:]))
                self.indices = np.concatenate((self.indices[:newindx], minor_index, self.indices[newindx:]))
                self.indptr  = self.indptr.copy()

                self.indptr[major_index+1:] += 1

            elif num_matches == 1:
                #entry appears exactly once
                self.data[start:end][indxs[0]] = val
            else:
                #entry appears more than once
                raise ValueError,'nonzero entry (%d,%d) occurs more than once' % (row,col)

            self.check_format(full_check=True)
        else:
            # We should allow slices here!
            raise IndexError, "invalid index"

    ######################
    # Conversion methods #
    ######################

    def todia(self):
        return self.tocoo(copy=False).todia()

    def todok(self):
        return self.tocoo(copy=False).todok()

    def tocoo(self,copy=True):
        """Return a COOrdinate representation of this matrix

        When copy=False the index and data arrays are not copied.
        """
        major_dim,minor_dim = self._swap(self.shape)

        data = self.data
        minor_indices = self.indices

        if copy:
            data = data.copy()
            minor_indices = minor_indices.copy()

        major_indices = np.empty(len(minor_indices), dtype=np.intc)

        sparsetools.expandptr(major_dim,self.indptr,major_indices)

        row,col = self._swap((major_indices,minor_indices))

        from coo import coo_matrix
        return coo_matrix((data,(row,col)), self.shape)

    def toarray(self):
        return self.tocoo(copy=False).toarray()

    ##############################################################
    # methods that examine or modify the internal data structure #
    ##############################################################

    def eliminate_zeros(self):
        """Remove zero entries from the matrix

        The is an *in place* operation
        """
        fn = sparsetools.csr_eliminate_zeros
        M,N = self._swap(self.shape)
        fn( M, N, self.indptr, self.indices, self.data)

        self.prune() #nnz may have changed

    def sum_duplicates(self):
        """Eliminate duplicate matrix entries by adding them together

        The is an *in place* operation
        """
        self.sort_indices()

        fn = sparsetools.csr_sum_duplicates
        M,N = self._swap(self.shape)
        fn( M, N, self.indptr, self.indices, self.data)

        self.prune() #nnz may have changed


    def __get_sorted(self):
        """Determine whether the matrix has sorted indices

        Returns
            - True: if the indices of the matrix are in sorted order
            - False: otherwise

        """

        #first check to see if result was cached
        if not hasattr(self,'__has_sorted_indices'):
            fn = sparsetools.csr_has_sorted_indices
            self.__has_sorted_indices = \
                    fn( len(self.indptr) - 1, self.indptr, self.indices)
        return self.__has_sorted_indices

    def __set_sorted(self, val):
        self.__has_sorted_indices = bool(val)

    has_sorted_indices = property(fget=__get_sorted, fset=__set_sorted)

    def sorted_indices(self):
        """Return a copy of this matrix with sorted indices
        """
        A = self.copy()
        A.sort_indices()
        return A

        # an alternative that has linear complexity is the following
        # although the previous option is typically faster
        #return self.toother().toother()

    def sort_indices(self):
        """Sort the indices of this matrix *in place*
        """

        if not self.has_sorted_indices:
            fn = sparsetools.csr_sort_indices
            fn( len(self.indptr) - 1, self.indptr, self.indices, self.data)
            self.has_sorted_indices = True

    #TODO remove after 0.7
    def ensure_sorted_indices(self, inplace=False):
        """Return a copy of this matrix where the column indices are sorted
        """
        warn('ensure_sorted_indices is deprecated, ' \
                'use sorted_indices() or sort_indices() instead', \
                DeprecationWarning)

        if inplace:
            self.sort_indices()
        else:
            return self.sorted_indices()

    def prune(self):
        """Remove empty space after all non-zero elements.
        """
        major_dim = self._swap(self.shape)[0]

        if len(self.indptr) != major_dim + 1:
            raise ValueError('index pointer has invalid length')
        if len(self.indices) < self.nnz:
            raise ValueError('indices array has fewer than nnz elements')
        if len(self.data) < self.nnz:
            raise ValueError('data array has fewer than nnz elements')

        self.data    = self.data[:self.nnz]
        self.indices = self.indices[:self.nnz]


    ###################
    # utility methods #
    ###################

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

    def _binopt(self, other, op):
        """apply the binary operation fn to two sparse matrices"""
        other = self.__class__(other)

        # e.g. csr_plus_csr, csr_minus_csr, etc.
        fn = getattr(sparsetools, self.format + op + self.format)

        maxnnz  = self.nnz + other.nnz
        indptr  = np.empty_like(self.indptr)
        indices = np.empty(maxnnz, dtype=np.intc)
        data    = np.empty(maxnnz, dtype=upcast(self.dtype,other.dtype))

        fn(self.shape[0], self.shape[1], \
                self.indptr,  self.indices,  self.data,
                other.indptr, other.indices, other.data,
                indptr, indices, data)

        actual_nnz = indptr[-1]
        indices = indices[:actual_nnz]
        data    = data[:actual_nnz]
        if actual_nnz < maxnnz / 2:
            #too much waste, trim arrays
            indices = indices.copy()
            data    = data.copy()

        A = self.__class__((data, indices, indptr), shape=self.shape)

        return A
