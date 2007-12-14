"""Sparse matrix classes using compressed storage
"""


__all__ = ['csr_matrix', 'csc_matrix', 'isspmatrix_csr', 'isspmatrix_csc' ]

import numpy
from numpy import array, matrix, asarray, asmatrix, zeros, rank, intc, \
        empty, hstack, isscalar, ndarray, shape, searchsorted

from base import spmatrix,isspmatrix
import sparsetools
from sparsetools import csrtodense, cootocsr, cootocsc, csctocsr, csrtocsc
from sputils import upcast, to_native, isdense, isshape, getdtype, \
        isscalarlike



def resize1d(arr, newlen):
    old = len(arr)
    new = zeros((newlen,), arr.dtype)
    new[:old] = arr
    return new




class _cs_matrix(spmatrix):
    """base matrix class for compressed row and column oriented matrices"""
    
    def __init__(self, arg1, dims=None, dtype=None, copy=False):
        spmatrix.__init__(self)

        if isdense(arg1):
            # Convert the dense array or matrix arg1 to sparse format
            if rank(arg1) == 1:
                # Convert to a row vector
                arg1 = arg1.reshape(1, arg1.shape[0])
            if rank(arg1) == 2:
                from coo import coo_matrix
                self.shape = arg1.shape
                self._set_self( self._tothis(coo_matrix(arg1)) )
            else:
                raise ValueError, "dense array must have rank 1 or 2"

        elif isspmatrix(arg1):
            if copy:
                arg1 = arg1.copy()
            self._set_self( self._tothis(arg1) )

        elif isinstance(arg1, tuple):
            if isshape(arg1):
                # It's a tuple of matrix dimensions (M, N)
                # create empty matrix
                self.shape = arg1   #spmatrix checks for errors here
                M, N = self.shape
                self.data    = zeros(0, getdtype(dtype, default=float))
                self.indices = zeros(0, intc)
                self.indptr  = zeros(self._swap(self.shape)[0] + 1, dtype='intc')
            else:
                try:
                    # Try interpreting it as (data, ij)
                    (data, ij) = arg1
                    assert isinstance(ij, ndarray) and (rank(ij) == 2) \
                           and (shape(ij) == (2, len(data)))
                except (AssertionError, TypeError, ValueError, AttributeError):
                    try:
                        # Try interpreting it as (data, indices, indptr)
                        (data, indices, indptr) = arg1
                    except:
                        raise ValueError, "unrecognized form for csr_matrix constructor"
                    else:
                        self.data    = array(data, copy=copy, dtype=getdtype(dtype, data))
                        self.indices = array(indices, copy=copy)
                        self.indptr  = array(indptr, copy=copy)
                else:
                    # (data, ij) format
                    from coo import coo_matrix
                    other = coo_matrix((data, ij), dims=dims )
                    other = self._tothis(other)
                    self._set_self( other )

        else:
            raise ValueError, "unrecognized form for csr_matrix constructor"

        # Read matrix dimensions given, if any
        if dims is not None:
            self.shape = dims   # spmatrix will check for errors
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

        self.check_format(full_check=False)

    def getnnz(self):
        return self.indptr[-1]
    nnz = property(fget=getnnz)
    
    def _get_dtype(self):
        return self.data.dtype
    def _set_dtype(self,newtype):
        self.data.dtype = newtype
    dtype = property(fget=_get_dtype,fset=_set_dtype)
    
    def _set_self(self, other, copy=False):
        """take the member variables of other and assign them to self"""

        if copy:
            other = other.copy()

        self.data    = other.data
        self.indices = other.indices
        self.indptr  = other.indptr
        self.shape   = other.shape
    
    def _check_format(self, full_check):
        self.shape = tuple([int(x) for x in self.shape])  # for floats etc.

        #use _swap to determine proper bounds
        major_name,minor_name = self._swap(('row','column'))
        major_dim,minor_dim = self._swap(self.shape)

        # index arrays should have integer data types
        if self.indptr.dtype.kind != 'i':
            warnings.warn("indptr array has non-integer dtype (%s)"  \
                                            % self.indptr.dtype.name )
        if self.indices.dtype.kind != 'i':
            warnings.warn("indices array has non-integer dtype (%s)" \
                                            % self.indices.dtype.name )

        # only support 32-bit ints for now
        self.indptr  = self.indptr.astype('intc')
        self.indices = self.indices.astype('intc')
        self.data    = to_native(self.data)

        # check array shapes
        if (rank(self.data) != 1) or (rank(self.indices) != 1) or \
           (rank(self.indptr) != 1):
            raise ValueError,"data, indices, and indptr should be rank 1"

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
                if amax(self.indices) >= minor_dim:
                    raise ValueError, "%s index values must be < %d" % \
                            (minor_name,minor_dim)
                if amin(self.indices) < 0:
                    raise ValueError, "%s index values must be >= 0" % \
                            minor_name
                if numpy.diff(self.indptr).min() < 0:
                    raise ValueError,'index pointer values must form a " \
                                        "non-decreasing sequence'

    def astype(self, t):
        return self._with_data(self.data.astype(t))

    def _with_data(self,data,copy=True):
        """Returns a matrix with the same sparsity structure as self,
        but with different data.  By default the structure arrays
        (i.e. .indptr and .indices) are copied.
        """
        if copy:
            return self.__class__((data,self.indices.copy(),self.indptr.copy()), \
                                   dims=self.shape,dtype=data.dtype)
        else:
            return self.__class__((data,self.indices,self.indptr), \
                                   dims=self.shape,dtype=data.dtype)

    def __abs__(self):
        return self._with_data(abs(self.data))

    def _real(self):
        return self._with_data(numpy.real(self.data))

    def _imag(self):
        return self._with_data(numpy.imag(self.data))

    def _binopt(self, other, op, in_shape=None, out_shape=None):
        """apply the binary operation fn to two sparse matrices"""
        other = self._tothis(other)

        if in_shape is None:
            in_shape = self.shape
        if out_shape is None:
            out_shape = self.shape

        # e.g. csr_plus_csr, cscmucsc, etc.
        fn = getattr(sparsetools, self.format + op + self.format)

        indptr, ind, data = fn(in_shape[0], in_shape[1], \
                               self.indptr, self.indices, self.data,
                               other.indptr, other.indices, other.data)
        return self.__class__((data, ind, indptr), dims=out_shape)


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
            raise NotImplementedError, 'adding a scalar to a CSC or CSR ' \
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
            raise NotImplementedError, 'adding a scalar to a CSC or CSR ' \
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

    def __imul__(self, other): #self *= other
        if isscalarlike(other):
            self.data *= other
            return self
        else:
            raise NotImplementedError

    def __neg__(self):
        return self._with_data(-self.data)

    def __truediv__(self,other):
        if isscalarlike(other):
            return self * (1./other)
        elif isspmatrix(other):
            if (other.shape != self.shape):
                raise ValueError, "inconsistent shapes"
            
            return self._binopt(other,'_eldiv_')
        else:
            raise NotImplementedError

    def __itruediv__(self, other): #self *= other
        if isscalarlike(other):
            recip = 1.0 / other
            self.data *= recip
            return self
        else:
            raise NotImplementedError

    def __pow__(self, other):
        """ Element-by-element power (unless other is a scalar, in which
        case return the matrix power.)
        """
        if isscalarlike(other):
            return self._with_data(self.data**other)
        elif isspmatrix(other):
            if (other.shape != self.shape):
                raise ValueError, "inconsistent shapes"
            
            return self._binopt(other,'_elmul_')
        else:
            raise NotImplementedError


    def matmat(self, other):
        if isspmatrix(other):
            M, K1 = self.shape
            K2, N = other.shape
            if (K1 != K2):
                raise ValueError, "shape mismatch error"
            other = self._tothis(other)

            return self._binopt(other,'mu',in_shape=(M,N),out_shape=(M,N))
        elif isdense(other):
            # TODO make sparse * dense matrix multiplication more efficient

            # matvec each column of other
            return hstack( [ self * col.reshape(-1,1) for col in other.T ] )
        else:
            raise TypeError, "need a dense or sparse matrix"

    def matvec(self, other):
        if isdense(other):
            if other.size != self.shape[1] or \
                    (other.ndim == 2 and self.shape[1] != other.shape[0]):
                raise ValueError, "dimension mismatch"

            # csrmux, cscmux
            fn = getattr(sparsetools,self.format + 'mux')
    
            #output array
            y = empty( self.shape[0], dtype=upcast(self.dtype,other.dtype) )

            fn(self.shape[0], self.shape[1], \
                self.indptr, self.indices, self.data, numpy.ravel(other), y)

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

    def rmatvec(self, other, conjugate=True):
        if conjugate:
            return self.transpose().conj() * other
        else:
            return self.transpose() * other

    def getdata(self, ind):
        return self.data[ind]

#    def _other_format(self):
#        if self.format == 'csr':
#            return 'csc'
#        elif self.format == 'csc':
#            return 'csr'
#        else:
#            raise TypeError,'unrecognized type'
#    
#    def _other__class__(self):
#        if self.format == 'csr':
#            return csc_matrix
#        elif self.format == 'csc':
#            return csr_matrix
#        else:
#            raise TypeError,'unrecognized type'



    def sum(self, axis=None):
        """Sum the matrix over the given axis.  If the axis is None, sum
        over both rows and columns, returning a scalar.
        """
        # The spmatrix base class already does axis=0 and axis=1 efficiently
        # so we only do the case axis=None here
        if axis is None:
            return self.data[:self.indptr[-1]].sum()
        else:
            return spmatrix.sum(self,axis)
            raise ValueError, "axis out of bounds"


    def copy(self):
        return self._with_data(self.data.copy(),copy=True)

    def _get_slice(self, i, start, stop, stride, dims):
        """Returns a view of the elements [i, myslice.start:myslice.stop].
        """
        if stride != 1:
            raise ValueError, "slicing with step != 1 not supported"
        if stop <= start:
            raise ValueError, "slice width must be >= 1"

        indices = []

        for ind in xrange(self.indptr[i], self.indptr[i+1]):
            if self.indices[ind] >= start and self.indices[ind] < stop:
                indices.append(ind)

        index = self.indices[indices] - start
        data   = self.data[indices]
        indptr = numpy.array([0, len(indices)])
        return self.__class__((data, index, indptr), dims=dims, \
                              dtype=self.dtype)


    def _transpose(self, cls, copy=False):
        M, N = self.shape
        return cls((self.data,self.indices,self.indptr),(N,M),copy=copy)
    
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

        major_indices = empty(len(minor_indices),dtype=intc)

        sparsetools.expandptr(major_dim,self.indptr,major_indices)

        row,col = self._swap((major_indices,minor_indices))

        from coo import coo_matrix
        return coo_matrix((data,(row,col)), self.shape)

    def conj(self, copy=False):
        return self._with_data(self.data.conj(),copy=copy)

    def sorted_indices(self):
        """Return a copy of this matrix with sorted indices
        """
        A = self.copy()
        A.sort_indices()
        return A

        # an alternative that has linear complexity is the following
        # typically the previous option is faster
        #return self._toother()._toother()

    def sort_indices(self):
        """Sort the indices of this matrix *in place*
        """
        fn = getattr(sparsetools,'sort_' + self.format + '_indices')

        M,N = self.shape
        fn( M, N, self.indptr, self.indices, self.data)

    def ensure_sorted_indices(self, inplace=False):
        """Return a copy of this matrix where the column indices are sorted
        """
        warnings.warn('ensure_sorted_indices is deprecated, ' \
                      'use sorted_indices() or sort_indices() instead', \
                      DeprecationWarning)
        
        if inplace:
            self.sort_indices()
        else:
            return self.sorted_indices()
    
    def prune(self):
        """ Remove empty space after all non-zero elements.
        """
        major_dim = self._swap(self.shape)[0]

        if len(self.indptr) != major_dim + 1:
            raise ValueError, "index pointer has invalid length"
        if len(self.indices) < self.nnz: 
            raise ValueError, "indices array has fewer than nnz elements"
        if len(self.data) < self.nnz:
            raise ValueError, "data array has fewer than nnz elements"
        
        self.data    = self.data[:self.nnz]
        self.indices = self.indices[:self.nnz]

    def _get_submatrix( self, shape0, shape1, slice0, slice1 ):
        """Return a submatrix of this matrix (new matrix is created)."""
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

            elif isscalar( sl ):
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
        return data, indices, indptr, i1 - i0, j1 - j0


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
    
    def check_format(self,full_check=True):
        """check whether matrix is in valid CSC format

            *Parameters*:
                full_check:
                    True  - rigorous check, O(N) operations : default
                    False - basic check, O(1) operations

        """
        _cs_matrix._check_format(self,full_check)

    def __getattr__(self, attr):
        if attr == 'rowind':
            warnings.warn("rowind attribute no longer in use. Use .indices instead",
                          DeprecationWarning)
            return self.indices
        else:
            return _cs_matrix.__getattr__(self, attr)

    def __iter__(self):
        csr = self.tocsr()
        for r in xrange(self.shape[0]):
            yield csr[r,:]

    def transpose(self, copy=False):
        return _cs_matrix._transpose(self, csr_matrix, copy)

    def __getitem__(self, key):
        if isinstance(key, tuple):
            #TODO use _swap() to unify this in _cs_matrix
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
            #this was implemented in fortran before - is there a noticable performance advantage?
            indxs = numpy.where(row == self.indices[self.indptr[col]:self.indptr[col+1]])
            if len(indxs[0]) == 0:
                return 0
            else:
                return self.data[self.indptr[col]:self.indptr[col+1]][indxs[0]]
        elif isintlike(key):
            # Was: return self.data[key]
            # If this is allowed, it should return the relevant row, as for
            # dense matrices (and __len__ should be supported again).
            raise IndexError, "integer index not supported for csc_matrix"
        else:
            raise IndexError, "invalid index"


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

    def _getslice(self, i, myslice):
        return self._getcolslice(i, myslice)

    def _getcolslice(self, myslice, j):
        """Returns a view of the elements [myslice.start:myslice.stop, j].
        """
        start, stop, stride = myslice.indices(self.shape[0])
        return _cs_matrix._get_slice(self, j, start, stop, stride, (stop - start, 1))

    def rowcol(self, ind):
        row = self.indices[ind]
        col = searchsorted(self.indptr, ind+1)-1
        return (row, col)


    def tocsc(self, copy=False):
        return self.toself(copy)
    
    def tocsr(self):
        indptr  = empty(self.shape[0] + 1, dtype=intc)
        indices = empty(self.nnz, dtype=intc)
        data    = empty(self.nnz, dtype=upcast(self.dtype))

        csctocsr(self.shape[0], self.shape[1], \
                self.indptr, self.indices, self.data, \
                indptr, indices, data)

        return csr_matrix((data, indices, indptr), self.shape)

    def toarray(self):
        return self.tocsr().toarray()

    def get_submatrix( self, slice0, slice1 ):
        """Return a submatrix of this matrix (new matrix is created).
        Rows and columns can be selected using slice instances, tuples,
        or scalars."""
        aux = _cs_matrix._get_submatrix( self, self.shape[1], self.shape[0],
                                         slice1, slice0 )
        nr, nc = aux[3:]
        return self.__class__( aux[:3], dims = (nc, nr) )
    
    # these functions are used by the parent class (_cs_matrix)
    # to remove redudancy between csc_matrix and csr_matrix
    def _swap(self,x):
        """swap the members of x if this is a column-oriented matrix
        """
        return (x[1],x[0])

    def _toother(self):
        return self.tocsr()

    def _tothis(self, other):
        return other.tocsc()

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

          - csr_matrix((data, ij), [dims=(M, N)])
            where data, ij satisfy:
                a[ij[0, k], ij[1, k]] = data[k]

          - csr_matrix((data, col, ptr), [dims=(M, N)])
            standard CSR representation
    """

    def check_format(self,full_check=True):
        """check whether matrix is in valid CSR format

            *Parameters*:
                full_check:
                    True  - perform rigorous checking - default
                    False - perform basic format check

        """
        _cs_matrix._check_format(self,full_check)

    def __getattr__(self, attr):
        if attr == 'colind':
            warnings.warn("colind attribute no longer in use. Use .indices instead",
                          DeprecationWarning)
            return self.indices
        else:
            return _cs_matrix.__getattr__(self, attr)

    def transpose(self, copy=False):
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
        return self.toself(copy)

    def tocsc(self):
        indptr  = empty(self.shape[1] + 1, dtype=intc)
        indices = empty(self.nnz, dtype=intc)
        data    = empty(self.nnz, dtype=upcast(self.dtype))

        csrtocsc(self.shape[0], self.shape[1], \
                 self.indptr, self.indices, self.data, \
                 indptr, indices, data)

        return csc_matrix((data, indices, indptr), self.shape)
    
    def toarray(self):
        data = numpy.zeros(self.shape, dtype=upcast(self.data.dtype))
        csrtodense(self.shape[0], self.shape[1], self.indptr, self.indices,
                   self.data, data)
        return data
    
    def get_submatrix( self, slice0, slice1 ):
        """Return a submatrix of this matrix (new matrix is created)..
        Rows and columns can be selected using slice instances, tuples,
        or scalars."""
        aux = _cs_matrix._get_submatrix( self, self.shape[0], self.shape[1],
                                         slice0, slice1 )
        nr, nc = aux[3:]
        return self.__class__( aux[:3], dims = (nr, nc) )

    # these functions are used by the parent class (_cs_matrix)
    # to remove redudancy between csc_matrix and csr_matrix
    def _swap(self,x):
        """swap the members of x if this is a column-oriented matrix
        """
        return (x[0],x[1])

    def _toother(self):
        return self.tocsc()

    def _tothis(self, other):
        return other.tocsr()


from sputils import _isinstance

def isspmatrix_csr(x):
    return _isinstance(x, csr_matrix)

def isspmatrix_csc(x):
    return _isinstance(x, csc_matrix)


