""" Scipy 2D sparse matrix module.

Original code by Travis Oliphant.
Modified and extended by Ed Schofield, Robert Cimrman, and Nathan Bell
"""


__all__ = ['spmatrix','csc_matrix','csr_matrix','coo_matrix',
            'lil_matrix','dok_matrix',
            'spdiags','speye','spidentity','spkron','extract_diagonal',
            'isspmatrix','issparse','isspmatrix_csc','isspmatrix_csr',
            'isspmatrix_lil','isspmatrix_dok', 'isspmatrix_coo',
            'lil_eye', 'lil_diags' ]

import warnings

from numpy import zeros, isscalar, real, imag, asarray, asmatrix, matrix, \
                  ndarray, amax, amin, rank, conj, searchsorted, ndarray,   \
                  less, where, greater, array, transpose, empty, ones, \
                  arange, shape, intc, clip, prod, unravel_index, hstack
import numpy
from scipy.sparse.sparsetools import cscmux, csrmux, \
     cootocsr, csrtocoo, cootocsc, csctocoo, csctocsr, csrtocsc, \
     densetocsr, csrtodense, \
     csrmucsr, cscmucsc, \
     csr_plus_csr, csc_plus_csc, csr_minus_csr, csc_minus_csc, \
     csr_elmul_csr, csc_elmul_csc, csr_eldiv_csr, csc_eldiv_csc

import sparsetools
import itertools, operator, copy
from bisect import bisect_left

def resize1d(arr, newlen):
    old = len(arr)
    new = zeros((newlen,), arr.dtype)
    new[:old] = arr
    return new

MAXPRINT = 50
ALLOCSIZE = 1000
NZMAX = 100

#TODO handle this in SWIG
def to_native(A):
    if not A.dtype.isnative:
        return A.astype(A.dtype.newbyteorder('native'))
    else:
        return A

# The formats that we might potentially understand.
_formats = {'csc':[0, "Compressed Sparse Column"],
            'csr':[1, "Compressed Sparse Row"],
            'dok':[2, "Dictionary Of Keys"],
            'lil':[3, "LInked List"],
            'dod':[4, "Dictionary of Dictionaries"],
            'sss':[5, "Symmetric Sparse Skyline"],
            'coo':[6, "COOrdinate"],
            'lba':[7, "Linpack BAnded"],
            'egd':[8, "Ellpack-itpack Generalized Diagonal"],
            'dia':[9, "DIAgonal"],
            'bsr':[10, "Block Sparse Row"],
            'msr':[11, "Modified compressed Sparse Row"],
            'bsc':[12, "Block Sparse Column"],
            'msc':[13, "Modified compressed Sparse Column"],
            'ssk':[14, "Symmetric SKyline"],
            'nsk':[15, "Nonsymmetric SKyline"],
            'jad':[16, "JAgged Diagonal"],
            'uss':[17, "Unsymmetric Sparse Skyline"],
            'vbr':[18, "Variable Block Row"],
            'und':[19, "Undefined"]
            }





class spmatrix(object):
    """ This class provides a base class for all sparse matrices.  It
    cannot be instantiated.  Most of the work is provided by subclasses.
    """

    __array_priority__ = 10.1
    ndim = 2
    def __init__(self, maxprint=MAXPRINT, allocsize=ALLOCSIZE):
        self.format = self.__class__.__name__[:3]
        self._shape = None
        if self.format == 'spm':
            raise ValueError, "This class is not intended" \
                  " to be instantiated directly."
        self.maxprint = maxprint
        self.allocsize = allocsize

    def set_shape(self,shape):
        shape = tuple(shape)

        if len(shape) != 2:
            raise ValueError("Only two-dimensional sparse arrays "
                                     "are supported.")
        try:
            shape = int(shape[0]),int(shape[1]) #floats, other weirdness
        except:
            raise TypeError,'invalid shape'

        if not (shape[0] >= 1 and shape[1] >= 1):
            raise TypeError,'invalid shape'

        if (self._shape != shape) and (self._shape is not None):
            try:
                self = self.reshape(shape)
            except NotImplementedError:
                raise NotImplementedError("Reshaping not implemented for %s." %
                                          self.__class__.__name__)
        self._shape = shape

    def get_shape(self):
        return self._shape

    shape = property(fget=get_shape, fset=set_shape)

    def reshape(self,shape):
        raise NotImplementedError

    def astype(self, t):
        csc = self.tocsc()
        return csc.astype(t)

    def asfptype(self):
        """Upcast matrix to a floating point format (if necessary)"""

        fp_types = ['f','d','F','D']

        if self.dtype.char in fp_types:
            return self
        else:
            for fp_type in fp_types:
                if self.dtype <= numpy.dtype(fp_type):
                    return self.astype(fp_type)

            raise TypeError,'cannot upcast [%s] to a floating \
                             point format' % self.dtype.name

    def __iter__(self):
        for r in xrange(self.shape[0]):
            yield self[r,:]

    def getmaxprint(self):
        try:
            maxprint = self.maxprint
        except AttributeError:
            maxprint = MAXPRINT
        return maxprint

    #def typecode(self):
    #    try:
    #        typ = self.dtype.char
    #    except AttributeError:
    #        typ = None
    #    return typ

    def getnnz(self):
        try:
            return self.nnz
        except AttributeError:
            raise AttributeError, "nnz not defined"

    def getnzmax(self):
        try:
            nzmax = self.nzmax
        except AttributeError:
            try:
                nzmax = self.nnz
            except AttributeError:
                nzmax = 0
        return nzmax

    def getformat(self):
        try:
            format = self.format
        except AttributeError:
            format = 'und'
        return format

    def rowcol(self, num):
        return (None, None)

    def getdata(self, num):
        return None

    def listprint(self, start, stop):
        """Provides a way to print over a single index.
        """
        return '\n'.join(['  %s\t%s' % (self.rowcol(ind), self.getdata(ind))
                         for ind in xrange(start,stop)]) + '\n'

    def __repr__(self):
        nnz = self.getnnz()
        format = self.getformat()
        return "<%dx%d sparse matrix of type '%s'\n" \
               "\twith %d stored elements in %s format>" % \
               (self.shape + (self.dtype.type, nnz, _formats[format][1]))

    def __str__(self):
        nnz = self.getnnz()
        maxprint = self.getmaxprint()
        val = ''
        if nnz > maxprint:
            val = val + self.listprint(0, maxprint/2)
            val = val + "  :\t:\n"
            val = val + self.listprint(nnz-maxprint//2, nnz)
        else:
            val = val + self.listprint(0, nnz)
        return val[:-1]

    def __nonzero__(self):  # Simple -- other ideas?
        return self.getnnz() > 0

    # What should len(sparse) return? For consistency with dense matrices,
    # perhaps it should be the number of rows?  But for some uses the number of
    # non-zeros is more important.  For now, raise an exception!
    def __len__(self):
        # return self.getnnz()
        raise TypeError, "sparse matrix length is ambiguous; use getnnz()" \
                         " or shape[0]"

    def asformat(self, format):
        # default converter goes through the CSC format
        csc = self.tocsc()
        return eval('%s_matrix' % format)(csc)

    # default operations use the CSC format as a base
    #   and operations return in csc format
    #  thus, a new sparse matrix format just needs to define
    #  a tocsc method

    def __abs__(self):
        csc = self.tocsc()
        return abs(csc)

    def __add__(self, other):   # self + other
        csc = self.tocsc()
        return csc.__add__(other)

    def __radd__(self, other):  # other + self
        csc = self.tocsc()
        return csc.__radd__(other)

    def __sub__(self, other):   # self - other
        #note: this can't be replaced by self + (-other) for unsigned types
        csc = self.tocsc()
        return csc.__sub__(other)

    def __rsub__(self, other):  # other - self
        csc = self.tocsc()
        return csc.__rsub__(other)

    def __mul__(self, other):
        csc = self.tocsc()
        return csc.__mul__(other)

    def __rmul__(self, other):
        csc = self.tocsc()
        return csc.__rmul__(other)

    def __truediv__(self, other):
        if isscalarlike(other):
            return self * (1./other)
        else:
            csc = self.tocsc()
            return csc.__truediv__(other)

    def __div__(self, other):
        # Always do true division
        return self.__truediv__(other)

    def __pow__(self, other):
        csc = self.tocsc()
        return csc ** other

    def __neg__(self):
        csc = self.tocsc()
        return -csc

    def __iadd__(self, other):
        raise NotImplementedError

    def __isub__(self, other):
        raise NotImplementedError

    def __imul__(self, other):
        raise NotImplementedError

    def __idiv__(self, other):
        return self.__itruediv__(other)

    def __itruediv__(self, other):
        raise NotImplementedError

    def __getattr__(self, attr):
        if attr == 'A':
            return self.toarray()
        elif attr == 'T':
            return self.transpose()
        elif attr == 'H':
            return self.getH()
        elif attr == 'real':
            return self._real()
        elif attr == 'imag':
            return self._imag()
        elif attr == 'size':
            return self.getnnz()
        else:
            raise AttributeError, attr + " not found"

    def transpose(self):
        csc = self.tocsc()
        return csc.transpose()

    def conj(self):
        csc = self.tocsc()
        return csc.conj()

    def conjugate(self):
        csc = self.tocsc()
        return csc.conj()

    # Renamed conjtranspose() -> getH() for compatibility with dense matrices
    def getH(self):
        return self.transpose().conj()

    def _real(self):
        csc = self.tocsc()
        return csc._real()

    def _imag(self):
        csc = self.tocsc()
        return csc._imag()

    def getcol(self, j):
        """Returns a copy of column j of the matrix, as an (m x 1) sparse
        matrix (column vector).
        """
        # Spmatrix subclasses should override this method for efficiency.
        # Post-multiply by a (n x 1) column vector 'a' containing all zeros
        # except for a_j = 1
        n = self.shape[1]
        a = csc_matrix((n, 1), dtype=self.dtype)
        a[j, 0] = 1
        return self * a

    def getrow(self, i):
        """Returns a copy of row i of the matrix, as a (1 x n) sparse
        matrix (row vector).
        """
        # Spmatrix subclasses should override this method for efficiency.
        # Pre-multiply by a (1 x m) row vector 'a' containing all zeros
        # except for a_i = 1
        m = self.shape[0]
        a = csr_matrix((1, m), dtype=self.dtype)
        a[0, i] = 1
        return a * self

    def dot(self, other):
        """ A generic interface for matrix-matrix or matrix-vector
        multiplication.
        """

        try:
            other.shape
        except AttributeError:
            # If it's a list or whatever, treat it like a matrix
            other = asmatrix(other)

        if isdense(other) and asarray(other).squeeze().ndim <= 1:
            # it's a dense row or column vector
            return self.matvec(other)
        elif len(other.shape) == 2:
            # it's a 2d dense array, dense matrix, or sparse matrix
            return self.matmat(other)
        else:
            raise ValueError, "could not interpret dimensions"


    def matmat(self, other):
        csc = self.tocsc()
        return csc.matmat(other)

    def matvec(self, other):
        """Multiplies the sparse matrix by the vector 'other', returning a
        dense vector as a result.
        """
        csc = self.tocsc()
        return csc.matvec(other)

    def rmatvec(self, other, conjugate=True):
        """Multiplies the vector 'other' by the sparse matrix, returning a
        dense vector as a result.

        If 'conjugate' is True:
            returns A.transpose().conj() * other
        Otherwise:
            returns A.transpose() * other.
        """
        csc = self.tocsc()
        return csc.rmatvec(other, conjugate=conjugate)

    #def rmatmat(self, other, conjugate=True):
    #    """ If 'conjugate' is True:
    #        returns other * A.transpose().conj(),
    #    where 'other' is a matrix.  Otherwise:
    #        returns other * A.transpose().
    #    """
    #    other = csc_matrix(other)
    #    if conjugate:
    #        return other.matmat(self.transpose()).conj()
    #    else:
    #        return other.matmat(self.transpose())

    def todense(self):
        return asmatrix(self.toarray())

    def toarray(self):
        csc = self.tocsc()
        return csc.toarray()

    def tocoo(self):
        csc = self.tocsc()
        return csc.tocoo()

    def tolil(self):
        return lil_matrix(self.tocsr())

    def toself(self, copy=False):
        if copy:
            new = self.copy()
        else:
            new = self
        return new

    def copy(self):
        csc = self.tocsc()
        return csc.copy()

    def sum(self, axis=None):
        """Sum the matrix over the given axis.  If the axis is None, sum
        over both rows and columns, returning a scalar.
        """
        # We use multiplication by an array of ones to achieve this.
        # For some sparse matrix formats more efficient methods are
        # possible -- these should override this function.
        m, n = self.shape
        if axis == 0:
            # sum over columns
            return asmatrix(ones((1, m), dtype=self.dtype)) * self
        elif axis == 1:
            # sum over rows
            return self * asmatrix(ones((n, 1), dtype=self.dtype))
        elif axis is None:
            # sum over rows and columns
            return ( self * asmatrix(ones((n, 1), dtype=self.dtype)) ).sum()
        else:
            raise ValueError, "axis out of bounds"

    def mean(self, axis=None):
        """Average the matrix over the given axis.  If the axis is None,
        average over both rows and columns, returning a scalar.
        """
        if axis == 0:
            mean = self.sum(0)
            mean *= 1.0 / self.shape[0]
            return mean
        elif axis == 1:
            mean = self.sum(1)
            mean *= 1.0 / self.shape[1]
            return mean
        elif axis is None:
            return self.sum(None) * 1.0 / (self.shape[0]*self.shape[1])
        else:
            raise ValueError, "axis out of bounds"

    def setdiag(self, values, k=0):
        """Fills the diagonal elements {a_ii} with the values from the
        given sequence.  If k != 0, fills the off-diagonal elements
        {a_{i,i+k}} instead.

        values may have any length.  If the diagonal is longer than values,
        then the remaining diagonal entries will not be set.  If values if
        longer than the diagonal, then the remaining values are ignored.
        """
        M, N = self.shape
        if (k > 0 and k >= N) or (k < 0 and -k >= M):
            raise ValueError, "k exceedes matrix dimensions"
        if k < 0:
            max_index = min(M+k, N, len(values))
            for i,v in enumerate(values[:max_index]):
                self[i - k, i] = v
        else:
            max_index = min(M, N-k, len(values))
            for i,v in enumerate(values[:max_index]):
                self[i, i + k] = v


    def save(self, file_name, format = '%d %d %f\n'):
        try:
            fd = open(file_name, 'w')
        except Exception, e:
            raise e, file_name

        fd.write('%d %d\n' % self.shape)
        fd.write('%d\n' % self.size)
        for ii in xrange(self.size):
            ir, ic = self.rowcol(ii)
            data = self.getdata(ii)
            fd.write(format % (ir, ic, data))
        fd.close()

class _cs_matrix(spmatrix):
    """base matrix class for compressed row and column oriented matrices"""

    def _set_self(self, other, copy=False):
        """take the member variables of other and assign them to self"""

        if copy:
            other = other.copy()

        self.data    = other.data
        self.indices = other.indices
        self.indptr  = other.indptr
        self.shape   = other.shape
        self.dtype   = other.data.dtype


    def _check_format(self, orientation, full_check):
        # some functions pass floats
        self.shape = tuple([int(x) for x in self.shape])

        assert(orientation in ['row','column'])
        if orientation == 'row':
            primary,secondary = 'row','column'
            indptr_size = self.shape[0] + 1
            indices_bound = self.shape[1]
        else:
            primary,secondary = 'column','row'
            indptr_size = self.shape[1] + 1
            indices_bound = self.shape[0]


        # index arrays should have integer data types
        if self.indptr.dtype.kind != 'i':
            warnings.warn("indptr array has non-integer dtype.  " \
                          "Casting from %s to int" % self.indptr.dtype.name )
            self.indptr = self.indptr.astype('i')
        if self.indices.dtype.kind != 'i':
            warnings.warn("indices array has non-integer dtype.  " \
                          "Casting from %s to int" % self.indices.dtype.name )
            self.indices = self.indices.astype('i')

        #TODO handle non-native byteorder better
        self.indptr  = to_native(self.indptr)
        self.indices = to_native(self.indices)
        self.data    = to_native(self.data)

        # set the data type
        self.dtype = self.data.dtype


        # check array shapes
        if (rank(self.data) != 1) or (rank(self.indices) != 1) or \
           (rank(self.indptr) != 1):
            raise ValueError,"data, indices, and indptr should be rank 1"

        # check index pointer
        if (len(self.indptr) != indptr_size ):
            raise ValueError, \
                "index pointer size (%d) should be (%d)" % \
                 (len(self.indptr), indptr_size)
        if (self.indptr[0] != 0):
            raise ValueError,"index pointer should start with 0"

        # check index and data arrays
        if (len(self.indices) != len(self.data)):
            raise ValueError,"indices and data should have the same size"
        if (self.indptr[-1] > len(self.indices)):
            raise ValueError, \
                  "Last value of index pointer should be less than "\
                  "the size of index and data arrays"

        self.nnz   =  self.indptr[-1]
        self.nzmax = len(self.indices)

        if full_check:
            #check format validity (more expensive)
            if self.nnz > 0:
                if amax(self.indices[:self.nnz]) >= indices_bound:
                    raise ValueError, "%s index values must be < %d" % \
                            (secondary,indices_bound)
                if amin(self.indices[:self.nnz]) < 0:
                    raise ValueError, "%s index values must be >= 0" % \
                            secondary
            if numpy.diff(self.indptr).min() < 0:
                raise ValueError,'index pointer values must form a " \
                                    "non-decreasing sequence'



    def astype(self, t):
        return self._with_data(self.data.astype(t))

    def __repr__(self):
        format = self.getformat()
        return "<%dx%d sparse matrix of type '%s'\n\twith %d stored "\
               "elements (space for %d)\n\tin %s format>" % \
               (self.shape + (self.dtype.type, self.getnnz(), self.nzmax, \
                   _formats[format][1]))

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

    def _binopt(self, other, fn, in_shape=None, out_shape=None):
        """apply the binary operation fn to two sparse matrices"""
        other = self._tothis(other)

        if in_shape is None:
            in_shape = self.shape
        if out_shape is None:
            out_shape = self.shape

        indptr, ind, data = fn(in_shape[0], in_shape[1], \
                               self.indptr, self.indices, self.data,
                               other.indptr, other.indices, other.data)
        return self.__class__((data, ind, indptr), dims=out_shape)


    def __add__(self,other,fn):
        # First check if argument is a scalar
        if isscalarlike(other):
            # Now we would add this scalar to every element.
            raise NotImplementedError, 'adding a scalar to a CSC or CSR ' \
                  'matrix is not supported'
        elif isspmatrix(other):
            if (other.shape != self.shape):
                raise ValueError, "inconsistent shapes"
            return self._binopt(other,fn)
        elif isdense(other):
            # Convert this matrix to a dense matrix and add them
            return self.todense() + other
        else:
            raise NotImplementedError

    def __radd__(self,other):
        return self.__add__(other)

    def __sub__(self,other,fn):
        # First check if argument is a scalar
        if isscalarlike(other):
            # Now we would add this scalar to every element.
            raise NotImplementedError, 'adding a scalar to a CSC or CSR ' \
                  'matrix is not supported'
        elif isspmatrix(other):
            if (other.shape != self.shape):
                raise ValueError, "inconsistent shapes"
            return self._binopt(other,fn)
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

    def __truediv__(self,other,fn):
        if isscalarlike(other):
            return self * (1./other)
        elif isspmatrix(other):
            other = self._tothis(other)
            if (other.shape != self.shape):
                raise ValueError, "inconsistent shapes"
            return self._binopt(other,fn)
        else:
            raise NotImplementedError

    def __itruediv__(self, other): #self *= other
        if isscalarlike(other):
            recip = 1.0 / other
            self.data *= recip
            return self
        else:
            raise NotImplementedError

    def __pow__(self, other, fn):
        """ Element-by-element power (unless other is a scalar, in which
        case return the matrix power.)
        """
        if isscalarlike(other):
            return self._with_data(self.data**other)
        elif isspmatrix(other):
            return self._binopt(other,fn)
        else:
            raise NotImplementedError


    def _matmat(self, other, fn):
        if isspmatrix(other):
            M, K1 = self.shape
            K2, N = other.shape
            if (K1 != K2):
                raise ValueError, "shape mismatch error"
            other = self._tothis(other)
            return self._binopt(other,fn,in_shape=(M,N),out_shape=(M,N))
        elif isdense(other):
            # TODO make sparse * dense matrix multiplication more efficient

            # matvec each column of other
            return hstack( [ self * col.reshape(-1,1) for col in other.T ] )
        else:
            raise TypeError, "need a dense or sparse matrix"

    def _matvec(self, other, fn):
        if isdense(other):
            if other.size != self.shape[1] or \
                    (other.ndim == 2 and self.shape[1] != other.shape[0]):
                raise ValueError, "dimension mismatch"

            y = fn(self.shape[0], self.shape[1], \
                   self.indptr, self.indices, self.data, numpy.ravel(other))

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

    def _tocoo(self, fn):
        rows, cols, data = fn(self.shape[0], self.shape[1], \
                              self.indptr, self.indices, self.data)
        return coo_matrix((data, (rows, cols)), self.shape)


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


    def conj(self, copy=False):
        return self._with_data(self.data.conj(),copy=copy)

    def _ensure_sorted_indices(self, shape0, shape1, inplace=False):
        """Return a copy of this matrix where the row indices are sorted
        """
        if inplace:
            sparsetools.sort_csr_indices(shape0, shape1,
                                         self.indptr, self.indices,
                                         self.data )
        else:
            return self._toother()._toother()

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

          - csc_matrix((M, N), [nzmax, dtype])
            to construct a container, where (M, N) are dimensions and
            nzmax, dtype are optional, defaulting to nzmax=sparse.NZMAX
            and dtype='d'.

          - csc_matrix((data, ij), [(M, N), nzmax])
            where data, ij satisfy:
                a[ij[0, k], ij[1, k]] = data[k]

          - csc_matrix((data, row, ptr), [(M, N)])
            standard CSC representation
    """
    def __init__(self, arg1, dims=None, nzmax=NZMAX, dtype=None, copy=False):
        _cs_matrix.__init__(self)

        if isdense(arg1):
            self.dtype = getdtype(dtype, arg1)
            # Convert the dense array or matrix arg1 to sparse format
            if rank(arg1) == 1:
                # Convert to a row vector
                arg1 = arg1.reshape(1, arg1.shape[0])
            if rank(arg1) == 2:
                self.shape = arg1.shape
                self.indptr, self.indices, self.data = \
                        densetocsr(arg1.shape[1], arg1.shape[0], arg1.T)
            else:
                raise ValueError, "dense array must have rank 1 or 2"

        elif isspmatrix(arg1):
            if copy:
                arg1 = arg1.copy()
            self._set_self( self._tothis(arg1) )

        elif isinstance(arg1, tuple):
            if isshape(arg1):
                # It's a tuple of matrix dimensions (M, N)
                self.shape = arg1   #spmatrix checks for errors here
                M, N = self.shape
                self.dtype   = getdtype(dtype, default=float)
                self.data    = zeros((nzmax,), self.dtype)
                self.indices = zeros((nzmax,), intc)
                self.indptr  = zeros((N+1,), intc)
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
                        raise ValueError, "unrecognized form for csc_matrix constructor"
                    else:
                        self.dtype = getdtype(dtype, data)
                        self.dtype   = getdtype(dtype, data)
                        self.data    = array(data, copy=copy, dtype=self.dtype)
                        self.indices = array(indices, copy=copy)
                        self.indptr  = array(indptr, copy=copy)
                else:
                    # (data, ij) format
                    other = coo_matrix((data, ij), dims=dims )
                    other = self._tothis(other)
                    self._set_self( other )

        else:
            raise ValueError, "unrecognized form for csc_matrix constructor"


        # Read matrix dimensions given, if any
        if dims is not None:
            self.shape = dims   # spmatrix will check for errors
        else:
            if self.shape is None:
                # shape not already set, try to infer dimensions
                try:
                    M = self.indices.max() + 1
                    N = len(self.indptr) - 1
                    self.shape = (M,N)
                except:
                    raise ValueError,'unable to infer matrix dimensions'

        self.check_format(full_check=False)

    def check_format(self,full_check=True):
        """check whether matrix is in valid CSC format

            *Parameters*:
                full_check:
                    True  - rigorous check, O(N) operations : default
                    False - basic check, O(1) operations

        """

        _cs_matrix._check_format(self,'column',full_check)

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

    def __add__(self, other):
        return _cs_matrix.__add__(self, other, csc_plus_csc)

    def __sub__(self, other):
        return _cs_matrix.__sub__(self, other, csc_minus_csc)

    def __truediv__(self,other):
        return _cs_matrix.__truediv__(self,other, csc_eldiv_csc)

    def __pow__(self, other):
        return _cs_matrix.__pow__(self, other, csc_elmul_csc)

    def transpose(self, copy=False):
        return _cs_matrix._transpose(self, csr_matrix, copy)

    def conj(self, copy=False):
        return _cs_matrix.conj(self, copy)

    def matvec(self, other):
        return _cs_matrix._matvec(self, other, cscmux)

    def matmat(self, other):
        return _cs_matrix._matmat(self, other, cscmucsc)

    def __getitem__(self, key):
        if isinstance(key, tuple):
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
                nzmax = self.nzmax
                if (nzmax < self.nnz+1):  # need more room
                    alloc = max(1, self.allocsize)
                    self.data = resize1d(self.data, nzmax + alloc)
                    self.indices = resize1d(self.indices, nzmax + alloc)

                newindex = self.indptr[col]
                self.data[newindex+1:]   = self.data[newindex:-1]
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

    def tocoo(self):
        return _cs_matrix._tocoo(self, csctocoo)

    def tocsr(self):
        indptr, colind, data = csctocsr(self.shape[0], self.shape[1], \
                                        self.indptr, self.indices, self.data)
        return csr_matrix((data, colind, indptr), self.shape)

    def _toother(self):
        return self.tocsr()

    def _tothis(self, other):
        return other.tocsc()

    def toarray(self):
        return self.tocsr().toarray()

    def prune(self):
        """ Remove empty space after all non-zero elements.
        """
        nnz = self.indptr[-1]
        if self.nzmax <= nnz:
            if self.nzmax < nnz:
                raise RuntimeError, "should never have nnz > nzmax"
            return
        self.nnz = nnz
        self.data = self.data[:nnz]
        self.indices = self.indices[:nnz]
        self.nzmax = nnz

    def ensure_sorted_indices(self, inplace=False):
        """Return a copy of this matrix where the row indices are sorted
        """
        return _cs_matrix._ensure_sorted_indices(self, self.shape[1], self.shape[0], inplace)

    def get_submatrix( self, slice0, slice1 ):
        """Return a submatrix of this matrix (new matrix is created).
        Rows and columns can be selected using slice instances, tuples,
        or scalars."""
        aux = _cs_matrix._get_submatrix( self, self.shape[1], self.shape[0],
                                         slice1, slice0 )
        nr, nc = aux[3:]
        return self.__class__( aux[:3], dims = (nc, nr) )

class csr_matrix(_cs_matrix):
    """ Compressed sparse row matrix
        This can be instantiated in several ways:
          - csr_matrix(d)
            with a dense matrix d

          - csr_matrix(s)
            with another sparse matrix s (sugar for .tocsr())

          - csr_matrix((M, N), [nzmax, dtype])
            to construct a container, where (M, N) are dimensions and
            nzmax, dtype are optional, defaulting to nzmax=sparse.NZMAX
            and dtype='d'.

          - csr_matrix((data, ij), [dims=(M, N), nzmax=nzmax])
            where data, ij satisfy:
                a[ij[0, k], ij[1, k]] = data[k]

          - csr_matrix((data, col, ptr), [dims=(M, N)])
            standard CSR representation
    """
    def __init__(self, arg1, dims=None, nzmax=NZMAX, dtype=None, copy=False):
        _cs_matrix.__init__(self)

        if isdense(arg1):
            self.dtype = getdtype(dtype, arg1)
            # Convert the dense array or matrix arg1 to sparse format
            if rank(arg1) == 1:
                # Convert to a row vector
                arg1 = arg1.reshape(1, arg1.shape[0])
            if rank(arg1) == 2:
                self.shape = arg1.shape
                self.indptr, self.indices, self.data = \
                        densetocsr(arg1.shape[0], arg1.shape[1], arg1)
            else:
                raise ValueError, "dense array must have rank 1 or 2"

        elif isspmatrix(arg1):
            if copy:
                arg1 = arg1.copy()
            self._set_self( self._tothis(arg1) )

        elif isinstance(arg1, tuple):
            if isshape(arg1):
                # It's a tuple of matrix dimensions (M, N)
                self.shape = arg1   #spmatrix checks for errors here
                M, N = self.shape
                self.dtype   = getdtype(dtype, default=float)
                self.data    = zeros((nzmax,), self.dtype)
                self.indices = zeros((nzmax,), intc)
                self.indptr  = zeros((M+1,), intc)
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
                        self.dtype   = getdtype(dtype, data)
                        self.data    = array(data, copy=copy, dtype=self.dtype)
                        self.indices = array(indices, copy=copy)
                        self.indptr  = array(indptr, copy=copy)
                else:
                    # (data, ij) format
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
                    M = len(self.indptr) - 1
                    N = self.indices.max() + 1
                    self.shape = (M,N)
                except:
                    raise ValueError,'unable to infer matrix dimensions'

        self.check_format(full_check=False)

    def check_format(self,full_check=True):
        """check whether matrix is in valid CSR format

            *Parameters*:
                full_check:
                    True  - perform rigorous checking - default
                    False - perform basic format check

        """

        _cs_matrix._check_format(self,'row',full_check)

    def __getattr__(self, attr):
        if attr == 'colind':
            warnings.warn("colind attribute no longer in use. Use .indices instead",
                          DeprecationWarning)
            return self.indices
        else:
            return _cs_matrix.__getattr__(self, attr)

    def __add__(self, other):
        return _cs_matrix.__add__(self, other, csr_plus_csr)

    def __sub__(self, other):
        return _cs_matrix.__sub__(self, other, csr_minus_csr)

    def __truediv__(self,other):
        return _cs_matrix.__truediv__(self,other, csr_eldiv_csr)

    def __pow__(self, other):
        return _cs_matrix.__pow__(self, other, csr_elmul_csr)

    def transpose(self, copy=False):
        return _cs_matrix._transpose(self, csc_matrix, copy)

    def conj(self, copy=False):
        return _cs_matrix.conj(self, copy)

    def matvec(self, other):
        return _cs_matrix._matvec(self, other, csrmux)

    def matmat(self, other):
        return _cs_matrix._matmat(self, other, csrmucsr)

    def __getitem__(self, key):
        if isinstance(key, tuple):
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
                nzmax = self.nzmax
                if (nzmax < self.nnz+1):  # need more room
                    alloc = max(1, self.allocsize)
                    self.data = resize1d(self.data, nzmax + alloc)
                    self.indices = resize1d(self.indices, nzmax + alloc)

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

    def tocsr(self, copy=False):
        return self.toself(copy)

    def tocoo(self):
        return _cs_matrix._tocoo(self, csrtocoo)

    def tocsc(self):
        indptr, rowind, data = csrtocsc(self.shape[0], self.shape[1], \
                                        self.indptr, self.indices, self.data)
        return csc_matrix((data, rowind, indptr), self.shape)

    def _toother(self):
        return self.tocsc()

    def _tothis(self, other):
        return other.tocsr()

    def toarray(self):
        data = numpy.zeros(self.shape, self.data.dtype)
        csrtodense(self.shape[0], self.shape[1], self.indptr, self.indices,
                   self.data, data)
        return data

    def prune(self):
        """ Eliminate non-zero entries, leaving space for at least
        newnzmax elements.
        """
        nnz = self.indptr[-1]
        if self.nzmax <= nnz:
            if self.nzmax < nnz:
                raise RuntimeError, "should never have nnz > nzmax"
            return
        self.data = self.data[:nnz]
        self.indices = self.indices[:nnz]
        self.nzmax = nnz

    def ensure_sorted_indices(self, inplace=False):
        """Return a copy of this matrix where the column indices are sorted
        """
        return _cs_matrix._ensure_sorted_indices(self, self.shape[0], self.shape[1], inplace)

    def get_submatrix( self, slice0, slice1 ):
        """Return a submatrix of this matrix (new matrix is created)..
        Rows and columns can be selected using slice instances, tuples,
        or scalars."""
        aux = _cs_matrix._get_submatrix( self, self.shape[0], self.shape[1],
                                         slice0, slice1 )
        nr, nc = aux[3:]
        return self.__class__( aux[:3], dims = (nr, nc) )

# This function was for sorting dictionary keys by the second tuple element.
# (We now use the Schwartzian transform instead for efficiency.)
# def csc_cmp(x, y):
#     if (x == y): return 0
#     elif (x[1] == y[1]):
#         if (x[0] > y[0]): return 1
#         elif (x[0] == y[0]): return 0
#         else: return -1
#     elif (x[1] > y[1]): return 1
#     else: return -1

# dictionary of keys based matrix
class dok_matrix(spmatrix, dict):
    """ A dictionary of keys based matrix.  This is an efficient
    structure for constructing sparse matrices for conversion to other
    sparse matrix types.
    """
    def __init__(self, A=None, shape=None, dtype=None):
        """ Create a new dictionary-of-keys sparse matrix.  An optional
        argument A is accepted, which initializes the dok_matrix with it.
        This can be a tuple of dimensions (M, N) or a (dense) array
        to copy.
        """
        dict.__init__(self)
        spmatrix.__init__(self,shape)
        self.dtype = getdtype(dtype, A, default=float)
        if A is not None:
            if isinstance(A, tuple):
                # Interpret as dimensions
                if not isshape(A):
                    raise TypeError, "dimensions must be a 2-tuple of positive"\
                            " integers"
                self.shape = A
            elif isspmatrix(A):
                # For sparse matrices, this is too inefficient; we need
                # something else.
                raise NotImplementedError, "initializing a dok_matrix with " \
                        "a sparse matrix is not yet supported"
            elif isdense(A):
                # Convert to a (1 x n) row vector
                if rank(A) == 1:
                    A = A.reshape(1, len(A))
                if rank(A) == 2:
                    M, N = A.shape
                    self.shape = (M, N)
                    for i in xrange(M):
                        for j in xrange(N):
                            if A[i, j] != 0:
                                self[i, j] = A[i, j]
                else:
                    raise ValueError, "dense array must have rank 1 or 2"
            else:
                raise TypeError, "argument should be a tuple of dimensions " \
                        "or a sparse or dense matrix"

    def getnnz(self):
        return dict.__len__(self)

    def __len__(self):
        return dict.__len__(self)

    def __str__(self):
        val = ''
        nnz = self.getnnz()
        keys = self.keys()
        keys.sort()
        if nnz > self.maxprint:
            for k in xrange(self.maxprint / 2):
                key = keys[k]
                val += "  %s\t%s\n" % (str(key), str(self[key]))
            val = val + "   :    \t  :\n"
            for k in xrange(nnz-self.maxprint/2, nnz):
                key = keys[k]
                val += "  %s\t%s\n" % (str(key), str(self[key]))
        else:
            for k in xrange(nnz):
                key = keys[k]
                val += "  %s\t%s\n" % (str(key), str(self[key]))
        return val[:-1]

    def get(self, key, default=0.):
        """This overrides the dict.get method, providing type checking
        but otherwise equivalent functionality.
        """
        try:
            i, j = key
            assert isintlike(i) and isintlike(j)
        except (AssertionError, TypeError, ValueError):
            raise IndexError, "index must be a pair of integers"
        try:
            assert not (i < 0 or i >= self.shape[0] or j < 0 or
                     j >= self.shape[1])
        except AssertionError:
            raise IndexError, "index out of bounds"
        return dict.get(self, key, default)

    def  __getitem__(self, key):
        """If key=(i,j) is a pair of integers, return the corresponding
        element.  If either i or j is a slice or sequence, return a new sparse
        matrix with just these elements.
        """
        try:
            i, j = key
        except (ValueError, TypeError):
            raise TypeError, "index must be a pair of integers or slices"


        # Bounds checking
        if isintlike(i):
            if i < 0:
                i += self.shape[0]
            if i < 0 or i >= self.shape[0]:
                raise IndexError, "index out of bounds"
        if isintlike(j):
            if j < 0:
                j += self.shape[1]
            if j < 0 or j >= self.shape[1]:
                raise IndexError, "index out of bounds"

        # First deal with the case where both i and j are integers
        if isintlike(i) and isintlike(j):
            return dict.get(self, key, 0.)
        else:
            # Either i or j is a slice, sequence, or invalid.  If i is a slice
            # or sequence, unfold it first and call __getitem__ recursively.

            if isinstance(i, slice):
                # Is there an easier way to do this?
                seq = xrange(i.start or 0, i.stop or self.shape[0], i.step or 1)
            elif operator.isSequenceType(i):
                seq = i
            else:
                # Make sure i is an integer. (But allow it to be a subclass of int).
                if not isintlike(i):
                    raise TypeError, "index must be a pair of integers or slices"
                seq = None
            if seq is not None:
                # i is a seq
                if isintlike(j):
                    # Create a new matrix of the correct dimensions
                    first = seq[0]
                    last = seq[-1]
                    if first < 0 or first >= self.shape[0] or last < 0 \
                                 or last >= self.shape[0]:
                        raise IndexError, "index out of bounds"
                    newshape = (last-first+1, 1)
                    new = dok_matrix(newshape)
                    # ** This uses linear time in the size m of dimension 0:
                    # new[0:seq[-1]-seq[0]+1, 0] = \
                    #         [self.get((element, j), 0) for element in seq]
                    # ** Instead just add the non-zero elements.  This uses
                    # ** linear time in the number of non-zeros:
                    for (ii, jj) in self.keys():
                        if jj == j and ii >= first and ii <= last:
                            dict.__setitem__(new, (ii-first, 0), \
                                             dict.__getitem__(self, (ii,jj)))
                else:
                    ###################################
                    # We should reshape the new matrix here!
                    ###################################
                    raise NotImplementedError, "fancy indexing supported over" \
                            " one axis only"
                return new

            # Below here, j is a sequence, but i is an integer
            if isinstance(j, slice):
                # Is there an easier way to do this?
                seq = xrange(j.start or 0, j.stop or self.shape[1], j.step or 1)
            elif operator.isSequenceType(j):
                seq = j
            else:
                # j is not an integer
                raise TypeError, "index must be a pair of integers or slices"

            # Create a new matrix of the correct dimensions
            first = seq[0]
            last = seq[-1]
            if first < 0 or first >= self.shape[1] or last < 0 \
                         or last >= self.shape[1]:
                raise IndexError, "index out of bounds"
            newshape = (1, last-first+1)
            new = dok_matrix(newshape)
            # ** This uses linear time in the size n of dimension 1:
            # new[0, 0:seq[-1]-seq[0]+1] = \
            #         [self.get((i, element), 0) for element in seq]
            # ** Instead loop over the non-zero elements.  This is slower
            # ** if there are many non-zeros
            for (ii, jj) in self.keys():
                if ii == i and jj >= first and jj <= last:
                    dict.__setitem__(new, (0, jj-first), \
                                     dict.__getitem__(self, (ii,jj)))
            return new


    def __setitem__(self, key, value):
        try:
            assert len(key) == 2
        except (AssertionError, TypeError):
            raise TypeError, "index must be a pair of integers, slices, or" \
                    " sequences"
        i, j = key


        # First deal with the case where both i and j are integers
        if isintlike(i) and isintlike(j):
            if i < 0:
                i += self.shape[0]
            if j < 0:
                j += self.shape[1]

            if i < 0 or i >= self.shape[0] or j < 0 or j >= self.shape[1]:
                raise IndexError, "index out of bounds"
            if isintlike(value) and value == 0:
                if key in self.keys():  # get rid of it something already there
                    del self[key]
            else:
                # Ensure value is a single element, not a sequence
                if isinstance(value, float) or isintlike(value) or \
                        isinstance(value, complex):
                    dict.__setitem__(self, key, self.dtype.type(value))
                    newrows = max(self.shape[0], int(key[0])+1)
                    newcols = max(self.shape[1], int(key[1])+1)
                    self.shape = (newrows, newcols)
                else:
                    raise TypeError, "cannot set matrix element to non-scalar"
            return                 # done
        else:
            # Either i or j is a slice, sequence, or invalid.  If i is a slice
            # or sequence, unfold it first and call __setitem__ recursively.
            if isinstance(i, slice):
                # Is there an easier way to do this?
                seq = xrange(i.start or 0, i.stop or self.shape[0], i.step or 1)
            elif operator.isSequenceType(i):
                seq = i
            else:
                # Make sure i is an integer. (But allow it to be a subclass of int).
                if not isintlike(i):
                    raise TypeError, "index must be a pair of integers or slices"
                seq = None
            if seq is not None:
                # First see if 'value' is another dok_matrix of the appropriate
                # dimensions
                if isinstance(value, dok_matrix):
                    if value.shape[1] == 1:
                        for element in seq:
                            self[element, j] = value[element, 0]
                    else:
                        raise NotImplementedError, "setting a 2-d slice of" \
                                " a dok_matrix is not yet supported"
                elif isscalar(value):
                    for element in seq:
                        self[element, j] = value
                else:
                    # See if value is a sequence
                    try:
                        if len(seq) != len(value):
                            raise ValueError, "index and value ranges must" \
                                              " have the same length"
                    except TypeError:
                        # Not a sequence
                        raise TypeError, "unsupported type for" \
                                         " dok_matrix.__setitem__"

                    # Value is a sequence
                    for element, val in itertools.izip(seq, value):
                        self[element, j] = val   # don't use dict.__setitem__
                            # here, since we still want to be able to delete
                            # 0-valued keys, do type checking on 'val' (e.g. if
                            # it's a rank-1 dense array), etc.
            else:
                # Process j
                if isinstance(j, slice):
                    seq = xrange(j.start or 0, j.stop or self.shape[1], j.step or 1)
                elif operator.isSequenceType(j):
                    seq = j
                else:
                    # j is not an integer
                    raise TypeError, "index must be a pair of integers or slices"

                # First see if 'value' is another dok_matrix of the appropriate
                # dimensions
                if isinstance(value, dok_matrix):
                    if value.shape[0] == 1:
                        for element in seq:
                            self[i, element] = value[0, element]
                    else:
                        raise NotImplementedError, "setting a 2-d slice of" \
                                " a dok_matrix is not yet supported"
                elif isscalar(value):
                    for element in seq:
                        self[i, element] = value
                else:
                    # See if value is a sequence
                    try:
                        if len(seq) != len(value):
                            raise ValueError, "index and value ranges must have" \
                                              " the same length"
                    except TypeError:
                        # Not a sequence
                        raise TypeError, "unsupported type for dok_matrix.__setitem__"
                    else:
                        for element, val in itertools.izip(seq, value):
                            self[i, element] = val


    def __add__(self, other):
        # First check if argument is a scalar
        if isscalarlike(other):
            new = dok_matrix(self.shape, dtype=self.dtype)
            # Add this scalar to every element.
            M, N = self.shape
            for i in xrange(M):
                for j in xrange(N):
                    aij = self.get((i, j), 0) + other
                    if aij != 0:
                        new[i, j] = aij
            #new.dtype.char = self.dtype.char
        elif isinstance(other, dok_matrix):
            if other.shape != self.shape:
                raise ValueError, "matrix dimensions are not equal"
            # We could alternatively set the dimensions to the the largest of
            # the two matrices to be summed.  Would this be a good idea?
            new = dok_matrix(self.shape, dtype=self.dtype)
            new.update(self)
            for key in other.keys():
                new[key] += other[key]
        elif isspmatrix(other):
            csc = self.tocsc()
            new = csc + other
        elif isdense(other):
            new = self.todense() + other
        else:
            raise TypeError, "data type not understood"
        return new

    def __radd__(self, other):
        # First check if argument is a scalar
        if isscalarlike(other):
            new = dok_matrix(self.shape, dtype=self.dtype)
            # Add this scalar to every element.
            M, N = self.shape
            for i in xrange(M):
                for j in xrange(N):
                    aij = self.get((i, j), 0) + other
                    if aij != 0:
                        new[i, j] = aij
        elif isinstance(other, dok_matrix):
            if other.shape != self.shape:
                raise ValueError, "matrix dimensions are not equal"
            new = dok_matrix(self.shape, dtype=self.dtype)
            new.update(self)
            for key in other:
                new[key] += other[key]
        elif isspmatrix(other):
            csc = self.tocsc()
            new = csc + other
        elif isdense(other):
            new = other + self.todense()
        else:
            raise TypeError, "data type not understood"
        return new

    def __neg__(self):
        new = dok_matrix(self.shape, dtype=self.dtype)
        for key in self.keys():
            new[key] = -self[key]
        return new

    def __mul__(self, other):           # self * other
        if isscalarlike(other):
            new = dok_matrix(self.shape, dtype=self.dtype)
            # Multiply this scalar by every element.
            for (key, val) in self.iteritems():
                new[key] = val * other
            #new.dtype.char = self.dtype.char
            return new
        else:
            return self.dot(other)

    def __imul__(self, other):           # self * other
        if isscalarlike(other):
            # Multiply this scalar by every element.
            for (key, val) in self.iteritems():
                self[key] = val * other
            #new.dtype.char = self.dtype.char
            return self
        else:
            return NotImplementedError

    def __rmul__(self, other):          # other * self
        if isscalarlike(other):
            new = dok_matrix(self.shape, dtype=self.dtype)
            # Multiply this scalar by every element.
            for (key, val) in self.iteritems():
                new[key] = other * val
            #new.dtype.char = self.dtype.char
            return new
        else:
            # Don't use asarray unless we have to
            try:
                tr = other.transpose()
            except AttributeError:
                tr = asarray(other).transpose()
            return self.transpose().dot(tr).transpose()

    def __truediv__(self, other):           # self * other
        if isscalarlike(other):
            new = dok_matrix(self.shape, dtype=self.dtype)
            # Multiply this scalar by every element.
            for (key, val) in self.iteritems():
                new[key] = val / other
            #new.dtype.char = self.dtype.char
            return new
        else:
            return self.tocsr() / other


    def __itruediv__(self, other):           # self * other
        if isscalarlike(other):
            # Multiply this scalar by every element.
            for (key, val) in self.iteritems():
                self[key] = val / other
            return self
        else:
            return NotImplementedError

    # What should len(sparse) return? For consistency with dense matrices,
    # perhaps it should be the number of rows?  For now it returns the number
    # of non-zeros.

    def transpose(self):
        """ Return the transpose
        """
        M, N = self.shape
        new = dok_matrix((N, M), dtype=self.dtype)
        for key, value in self.iteritems():
            new[key[1], key[0]] = value
        return new

    def conjtransp(self):
        """ Return the conjugate transpose
        """
        M, N = self.shape
        new = dok_matrix((N, M), dtype=self.dtype)
        for key, value in self.iteritems():
            new[key[1], key[0]] = conj(value)
        return new

    def copy(self):
        new = dok_matrix(self.shape, dtype=self.dtype)
        new.update(self)
        new.shape = self.shape
        return new

    def take(self, cols_or_rows, columns=1):
        # Extract columns or rows as indictated from matrix
        # assume cols_or_rows is sorted
        new = dok_matrix(dtype=self.dtype)    # what should the dimensions be ?!
        indx = int((columns == 1))
        N = len(cols_or_rows)
        if indx: # columns
            for key in self.keys():
                num = searchsorted(cols_or_rows, key[1])
                if num < N:
                    newkey = (key[0], num)
                    new[newkey] = self[key]
        else:
            for key in self.keys():
                num = searchsorted(cols_or_rows, key[0])
                if num < N:
                    newkey = (num, key[1])
                    new[newkey] = self[key]
        return new

    def split(self, cols_or_rows, columns=1):
        # Similar to take but returns two arrays, the extracted columns plus
        # the resulting array.  Assumes cols_or_rows is sorted
        base = dok_matrix()
        ext = dok_matrix()
        indx = int((columns == 1))
        if indx:
            for key in self.keys():
                num = searchsorted(cols_or_rows, key[1])
                if cols_or_rows[num] == key[1]:
                    newkey = (key[0], num)
                    ext[newkey] = self[key]
                else:
                    newkey = (key[0], key[1]-num)
                    base[newkey] = self[key]
        else:
            for key in self.keys():
                num = searchsorted(cols_or_rows, key[0])
                if cols_or_rows[num] == key[0]:
                    newkey = (num, key[1])
                    ext[newkey] = self[key]
                else:
                    newkey = (key[0]-num, key[1])
                    base[newkey] = self[key]
        return base, ext


    def matvec(self, other):
        if isdense(other):
            if other.shape[0] != self.shape[1]:
                raise ValueError, "dimensions do not match"
            new = [0] * self.shape[0]
            for key in self.keys():
                new[int(key[0])] += self[key] * other[int(key[1]), ...]
            new = array(new)
            if isinstance(other, matrix):
                new = asmatrix(new)
                # Do we need to return the transpose?
                if other.shape[1] == 1:
                    new = new.T
            return new
        else:
            raise TypeError, "need a dense vector"

    def rmatvec(self, other, conjugate=True):
        if isdense(other):
            if other.shape[-1] != self.shape[0]:
                raise ValueError, "dimensions do not match"
            new = [0] * self.shape[1]
            for key in self.keys():
                new[int(key[1])] += other[..., int(key[0])] * conj(self[key])
            new = array(new)
            if isinstance(other, matrix):
                new = asmatrix(new)
                # Do we need to return the transpose?
                if other.shape[1] == 1:
                    new = new.T
            return new
        else:
            raise TypeError, "need a dense vector"

    def tocsr(self, nzmax=None):
        """ Return Compressed Sparse Row format arrays for this matrix
        """
        keys = self.keys()
        keys.sort()

        nnz = len(keys)
        nzmax = max(nnz, nzmax)
        data = zeros(nzmax, dtype=self.dtype)
        colind = zeros(nzmax, dtype=intc)
        # Empty rows will leave row_ptr dangling.  We assign row_ptr[i]
        # for each empty row i to point off the end.  Is this sufficient??
        row_ptr = empty(self.shape[0]+1, dtype=intc)
        row_ptr[:] = nnz
        current_row = -1
        k = 0
        for key in keys:
            ikey0 = int(key[0])
            ikey1 = int(key[1])
            if ikey0 != current_row:
                row_ptr[current_row+1:ikey0+1] = k
                current_row = ikey0
            data[k] = dict.__getitem__(self, key)
            colind[k] = ikey1
            k += 1
        data = array(data)
        colind = array(colind)
        row_ptr = array(row_ptr)
        return csr_matrix((data, colind, row_ptr), dims=self.shape, nzmax=nzmax)

    def tocsc(self, nzmax=None):
        """ Return Compressed Sparse Column format arrays for this matrix
        """
        # Fast sort on columns using the Schwartzian transform
        keys = [(k[1], k[0]) for k in self.keys()]
        keys.sort()
        keys = [(k[1], k[0]) for k in keys]

        nnz = len(keys)
        nzmax = max(nnz, nzmax)
        data = zeros(nzmax, dtype=self.dtype)
        rowind = zeros(nzmax, dtype=intc)
        # Empty columns will leave col_ptr dangling.  We assign col_ptr[j]
        # for each empty column j to point off the end.  Is this sufficient??
        col_ptr = empty(self.shape[1]+1, dtype=intc)
        col_ptr[:] = nnz
        current_col = -1
        k = 0
        for key in keys:
            ikey0 = int(key[0])
            ikey1 = int(key[1])
            if ikey1 != current_col:
                col_ptr[current_col+1:ikey1+1] = k
                current_col = ikey1
            data[k] = self[key]
            rowind[k] = ikey0
            k += 1
        return csc_matrix((data, rowind, col_ptr), dims=self.shape, nzmax=nzmax)

    def toarray(self):
        new = zeros(self.shape, dtype=self.dtype)
        for key in self.keys():
            ikey0 = int(key[0])
            ikey1 = int(key[1])
            new[ikey0, ikey1] = self[key]
        if abs(new.imag).max() == 0:
            new = new.real
        return new

    def resize(self, shape):
        """ Resize the matrix to dimensions given by 'shape', removing any
        non-zero elements that lie outside.
        """
        if not isshape(shape):
            raise TypeError, "dimensions must be a 2-tuple of positive"\
                             " integers"
        newM, newN = shape
        M, N = self.shape
        if newM < M or newN < N:
            # Remove all elements outside new dimensions
            for (i, j) in self.keys():
                if i >= newM or j >= newN:
                    del self[i, j]
        self.shape = shape



class coo_matrix(spmatrix):
    """ A sparse matrix in coordinate list format.

    COO matrices are created either as:
        A = coo_matrix(None, dims=(m, n), [dtype])
    for a zero matrix, or as:
        A = coo_matrix(M)
    where M is a dense matrix or rank 2 ndarray, or as:
        A = coo_matrix((obj, ij), [dims])
    where the dimensions are optional.  If supplied, we set (M, N) = dims.
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
    def __init__(self, arg1, dims=None, dtype=None):
        spmatrix.__init__(self)
        if isinstance(arg1, tuple):
            try:
                obj, ij = arg1
            except:
                raise TypeError, "invalid input format"

            try:
                if len(ij) != 2:
                    raise TypeError
            except TypeError:
                raise TypeError, "invalid input format"

            self.row = asarray(ij[0])
            self.col = asarray(ij[1])
            self.data = asarray(obj)
            self.dtype = self.data.dtype

            if dims is None:
                if len(self.row) == 0 or len(self.col) == 0:
                    raise ValueError, "cannot infer dimensions from zero sized index arrays"
                M = self.row.max() + 1
                N = self.col.max() + 1
                self.shape = (M, N)
            else:
                # Use 2 steps to ensure dims has length 2.
                M, N = dims
                self.shape = (M, N)

        elif arg1 is None:      # clumsy!  We should make ALL arguments
                                # keyword arguments instead!
            # Initialize an empty matrix.
            if not isinstance(dims, tuple) or not isintlike(dims[0]):
                raise TypeError, "dimensions not understood"
            self.shape = dims
            self.dtype = getdtype(dtype, default=float)
            self.data = array([])
            self.row = array([],dtype=intc)
            self.col = array([],dtype=intc)
        else:
            #dense argument
            try:
                M = asarray(arg1)
            except:
                raise TypeError, "invalid input format"

            if len(M.shape) != 2:
                raise TypeError, "expected rank 2 array or matrix"
            self.shape = M.shape
            self.dtype = M.dtype
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
            warnings.warn("row index array has non-integer dtype.  " \
                          "Casting from %s to int" % self.row.dtype.name )
            self.row = self.row.astype('i')
        if self.col.dtype.kind != 'i':
            warnings.warn("column index array has non-integer dtype.  " \
                          "Casting from %s to int" % self.col.dtype.name )
            self.col = self.col.astype('i')

        #TODO do this in SWIG
        self.row  = to_native(self.row)
        self.col  = to_native(self.col)
        self.data = to_native(self.data)

        if nnz > 0:
            if(amax(self.row) >= self.shape[0]):
                raise ValueError, "row index exceedes matrix dimensions"
            if(amax(self.col) >= self.shape[1]):
                raise ValueError, "column index exceedes matrix dimensions"
            if(amin(self.row) < 0):
                raise ValueError, "negative row index found"
            if(amin(self.col) < 0):
                raise ValueError, "negative column index found"

        # some functions pass floats
        self.shape = tuple([int(x) for x in self.shape])
        self.nnz = nnz

##    _normalize shouldn't be necessary anymore

##    def _normalize(self, rowfirst=False):
##        if rowfirst:
##            #sort by increasing rows first, columns second
##            if getattr(self, '_is_normalized', None):
##                #columns already sorted, use stable sort for rows
##                P = numpy.argsort(self.row, kind='mergesort')
##                return self.data[P], self.row[P], self.col[P]
##            else:
##                #nothing already sorted
##                P  = numpy.lexsort(keys=(self.col, self.row))
##                return self.data[P], self.row[P], self.col[P]
##        if getattr(self, '_is_normalized', None):
##            return self.data, self.row, self.col
##        #sort by increasing rows first, columns second
##        P  = numpy.lexsort(keys=(self.row, self.col))
##        self.data, self.row, self.col = self.data[P], self.row[P], self.col[P]
##        setattr(self, '_is_normalized', 1)
##        return self.data, self.row, self.col

    def rowcol(self, num):
        return (self.row[num], self.col[num])

    def getdata(self, num):
        return self.data[num]

    def tocsc(self):
        if self.nnz == 0:
            return csc_matrix(self.shape, dtype=self.dtype)
        else:
            indptr, rowind, data = cootocsc(self.shape[0], self.shape[1], \
                                            self.nnz, self.row, self.col, \
                                            self.data)
            return csc_matrix((data, rowind, indptr), self.shape)


    def tocsr(self):
        if self.nnz == 0:
            return csr_matrix(self.shape, dtype=self.dtype)
        else:
            indptr, colind, data = cootocsr(self.shape[0], self.shape[1], \
                                            self.nnz, self.row, self.col, \
                                            self.data)
            return csr_matrix((data, colind, indptr), self.shape)

    def tocoo(self, copy=False):
        return self.toself(copy)


class lil_matrix(spmatrix):
    """Row-based linked list matrix, by Ed Schofield.

    This contains a list (self.rows) of rows, each of which is a sorted
    list of column indices of non-zero elements. It also contains a list
    (self.data) of lists of these elements.
    """

    def __init__(self, A=None, shape=None, dtype=None):
        """ Create a new list-of-lists sparse matrix.  An optional
        argument A is accepted, which initializes the lil_matrix with it.
        This can be a tuple of dimensions (M, N) or a dense array /
        matrix to copy, or a sparse matrix of the following types:
          - csr_matrix
          - lil_matrix
        """
        spmatrix.__init__(self)
        self.dtype = getdtype(dtype, A, default=float)

        # First get the shape
        if A is None:
            if not isshape(shape):
                raise TypeError, "need a valid shape"
            M, N = shape
        else:
            if isshape(A):
                M, N = A
                A = None
            else:
                if not isdense(A):
                    # A is not dense.  If it's a spmatrix, ensure it's a
                    # csr_matrix or lil_matrix, which are the only types that
                    # support row slices.
                    if isinstance(A, spmatrix):
                        if not isinstance(A, lil_matrix) and \
                                not isinstance(A, csr_matrix):
                            raise TypeError, "unsupported matrix type"

                    # Otherwise, try converting to a matrix.  So if it's
                    # a list (rank 1), it will become a row vector
                    else:
                        try:
                            A = asmatrix(A)
                        except TypeError:
                            raise TypeError, "unsupported matrix type"
                elif rank(A) == 1:
                    # Construct a row vector
                    A = asmatrix(A)
                if rank(A) != 2:
                    raise ValueError, "can only initialize with a rank 1 or" \
                            " 2 array"
                if shape is None:
                    shape = (None, None)   # simplifies max() operation
                M = max(shape[0], A.shape[0])
                N = max(shape[1], A.shape[1])
        self.shape = (M, N)

        # Pluck out all non-zeros from the dense array/matrix A
        self.rows = numpy.empty((M,), dtype=object)
        self.data = numpy.empty((M,), dtype=object)
        for i in range(M):
            self.rows[i] = []
            self.data[i] = []

        if A is not None:
            for i in xrange(A.shape[0]):
                self[i, :] = A[i, :]

    def __iadd__(self,other):
        self[:,:] = self + other
        return self

    def __isub__(self,other):
        self[:,:] = self - other
        return self

    def __imul__(self,other):
        if isscalarlike(other):
            self[:,:] = self * other
            return self
        else:
            raise NotImplementedError

    def __itruediv__(self,other):
        if isscalarlike(other):
            self[:,:] = self / other
            return self
        else:
            raise NotImplementedError

    # Whenever the dimensions change, empty lists should be created for each
    # row

    def getnnz(self):
        return sum([len(rowvals) for rowvals in self.data])

    def __str__(self):
        val = ''
        for i, row in enumerate(self.rows):
            for pos, j in enumerate(row):
                val += "  %s\t%s\n" % (str((i, j)), str(self.data[i][pos]))
        return val[:-1]

    #def __repr__(self):
    #    format = self.getformat()
    #    return "<%dx%d sparse matrix with %d stored "\
    #           "elements in %s format>" % \
    #           (self.shape + (self.getnnz(), _formats[format][1]))

    def getrowview(self, i):
        """Returns a view of the 'i'th row (without copying).
        """
        new = lil_matrix((1, self.shape[1]), dtype=self.dtype)
        new.rows[0] = self.rows[i]
        new.data[0] = self.data[i]
        return new

    def getrow(self, i):
        """Returns a copy of the 'i'th row.
        """
        new = lil_matrix((1, self.shape[1]), dtype=self.dtype)
        new.rows[0] = self.rows[i][:]
        new.data[0] = self.data[i][:]
        return new

    def _get1(self, i, j):
        row = self.rows[i]
        data = self.data[i]

        if j < 0:
            j += self.shape[1]

        if j < 0 or j > self.shape[1]:
            raise IndexError,'column index out of bounds'

        pos = bisect_left(row, j)
        if pos != len(data) and row[pos] == j:
            return data[pos]
        else:
            return 0

    def _slicetoseq(self, j, shape):
        if j.start is not None and j.start < 0:
            start =  shape + j.start
        elif j.start is None:
            start = 0
        else:
            start = j.start
        if j.stop is not None and j.stop < 0:
            stop = shape + j.stop
        elif j.stop is None:
            stop = shape
        else:
            stop = j.stop
        j = range(start, stop, j.step or 1)
        return j


    def __getitem__(self, index):
        """Return the element(s) index=(i, j), where j may be a slice.
        This always returns a copy for consistency, since slices into
        Python lists return copies.
        """
        try:
            i, j = index
        except (AssertionError, TypeError):
            raise IndexError, "invalid index"

        if isscalar(i):
            if isscalar(j):
                return self._get1(i, j)
            if isinstance(j, slice):
                j = self._slicetoseq(j, self.shape[1])
            if issequence(j):
                return self.__class__([[self._get1(i, jj) for jj in j]])
        elif issequence(i) and issequence(j):
            return self.__class__([[self._get1(ii, jj) for (ii, jj) in zip(i, j)]])
        elif issequence(i) or isinstance(i, slice):
            if isinstance(i, slice):
                i = self._slicetoseq(i, self.shape[0])
            if isscalar(j):
                return self.__class__([[self._get1(ii, j)] for ii in i])
            if isinstance(j, slice):
                j = self._slicetoseq(j, self.shape[1])
            if issequence(j):
                return self.__class__([[self._get1(ii, jj) for jj in j] for ii in i])
        else:
            raise IndexError


    def _insertat(self, i, j, x):
        """ helper for __setitem__: insert a value at (i,j) where i, j and x
        are all scalars """
        row = self.rows[i]
        data = self.data[i]
        self._insertat2(row, data, j, x)

    def _insertat2(self, row, data, j, x):
        """ helper for __setitem__: insert a value in the given row/data at
        column j. """

        if j < 0: #handle negative column indices
            j += self.shape[1]

        if j < 0 or j >= self.shape[1]:
            raise IndexError,'column index out of bounds'

        pos = bisect_left(row, j)
        if x != 0:
            if pos == len(row):
                row.append(j)
                data.append(x)
            elif row[pos] != j:
                row.insert(pos, j)
                data.insert(pos, x)
            else:
                data[pos] = x
        else:
            if pos < len(row) and row[pos] == j:
                del row[pos]
                del data[pos]


    def _insertat3(self, row, data, j, x):
        """ helper for __setitem__ """
        if isinstance(j, slice):
            j = self._slicetoseq(j, self.shape[1])
        if issequence(j):
            if isinstance(x, spmatrix):
                x = x.todense()
            x = numpy.asarray(x).squeeze()
            if isscalar(x) or x.size == 1:
                for jj in j:
                    self._insertat2(row, data, jj, x)
            else:
                # x must be one D. maybe check these things out
                for jj, xx in zip(j, x):
                    self._insertat2(row, data, jj, xx)
        elif isscalar(j):
            self._insertat2(row, data, j, x)
        else:
            raise ValueError, "invalid column value: %s" % str(j)


    def __setitem__(self, index, x):
        if isscalar(x):
            x = self.dtype.type(x)
        elif not isinstance(x, spmatrix):
            x = numpy.asarray(x, dtype=self.dtype)

        try:
            i, j = index
        except (ValueError, TypeError):
            raise IndexError, "invalid index"

        if isscalar(i):
            row = self.rows[i]
            data = self.data[i]
            self._insertat3(row, data, j, x)
        elif issequence(i) and issequence(j):
            if isscalar(x):
                for ii, jj in zip(i, j):
                    self._insertat(ii, jj, x)
            else:
                for ii, jj, xx in zip(i, j, x):
                    self._insertat(ii, jj, xx)
        elif isinstance(i, slice) or issequence(i):
            rows = self.rows[i]
            datas = self.data[i]
            if isscalar(x):
                for row, data in zip(rows, datas):
                    self._insertat3(row, data, j, x)
            else:
                for row, data, xx in zip(rows, datas, x):
                    self._insertat3(row, data, j, xx)
        else:
            raise ValueError, "invalid index value: %s" % str((i, j))

    def __mul__(self, other):           # self * other
        if isscalarlike(other):
            if other == 0:
                # Multiply by zero: return the zero matrix
                new = lil_matrix(shape=self.shape, dtype=self.dtype)
            else:
                new = self.copy()
                # Multiply this scalar by every element.
                new.data = numpy.array([[val*other for val in rowvals] for
                                        rowvals in new.data], dtype=object)
            return new
        else:
            return self.dot(other)

    def __truediv__(self, other):           # self / other
        if isscalarlike(other):
            new = self.copy()
            # Divide every element by this scalar
            new.data = numpy.array([[val/other for val in rowvals] for
                                    rowvals in new.data], dtype=object)
            return new
        else:
            return self.tocsr() / other

    def multiply(self, other):
        """Point-wise multiplication by another lil_matrix.

        """
        if isscalar(other):
            return self.__mul__(other)

        if isspmatrix_lil(other):
            reference,target = self,other

            if reference.shape != target.shape:
                raise ValueError("Dimensions do not match.")

            if len(reference.data) > len(target.data):
                reference,target = target,reference

            new = lil_matrix(reference.shape)
            for r,row in enumerate(reference.rows):
                tr = target.rows[r]
                td = target.data[r]
                rd = reference.data[r]
                L = len(tr)
                for c,column in enumerate(row):
                    ix = bisect_left(tr,column)
                    if ix < L and tr[ix] == column:
                        new.rows[r].append(column)
                        new.data[r].append(rd[c] * td[ix])
            return new
        else:
            raise ValueError("Point-wise multiplication only allowed "
                             "with another lil_matrix.")

    def copy(self):
        new = lil_matrix(self.shape, dtype=self.dtype)
        new.data = copy.deepcopy(self.data)
        new.rows = copy.deepcopy(self.rows)
        return new

    def reshape(self,shape):
        new = lil_matrix(shape,dtype=self.dtype)
        j_max = self.shape[1]
        for i,row in enumerate(self.rows):
            for col,j in enumerate(row):
                new_r,new_c = unravel_index(i*j_max + j,shape)
                new[new_r,new_c] = self[i,j]
        return new

    def __add__(self, other):
        if isscalar(other) and other != 0:
            raise ValueError("Refusing to destroy sparsity. "
                             "Use x.todense() + c instead.")
        else:
            return spmatrix.__add__(self, other)

    def __rmul__(self, other):          # other * self
        if isscalarlike(other):
            # Multiplication by a scalar is symmetric
            return self.__mul__(other)
        else:
            return spmatrix.__rmul__(self, other)


    def toarray(self):
        d = zeros(self.shape, dtype=self.dtype)
        for i, row in enumerate(self.rows):
            for pos, j in enumerate(row):
                d[i, j] = self.data[i][pos]
        return d

    def transpose(self):
        """ Return the transpose as a csc_matrix.
        """
        # Overriding the spmatrix.transpose method here prevents an unnecessary
        # csr -> csc conversion
        return self.tocsr().transpose()

    def tocsr(self, nzmax=None):
        """ Return Compressed Sparse Row format arrays for this matrix.
        """
        nnz = self.getnnz()
        nzmax = max(nnz, nzmax)
        data = zeros(nzmax, dtype=self.dtype)
        colind = zeros(nzmax, dtype=intc)
        row_ptr = empty(self.shape[0]+1, dtype=intc)
        row_ptr[:] = nnz
        k = 0
        for i, row in enumerate(self.rows):
            data[k : k+len(row)] = self.data[i]
            colind[k : k+len(row)] = self.rows[i]
            row_ptr[i] = k
            k += len(row)

        row_ptr[-1] = nnz           # last row number + 1
        return csr_matrix((data, colind, row_ptr), dims=self.shape, nzmax=nzmax)

    def tocsc(self, nzmax=None):
        """ Return Compressed Sparse Column format arrays for this matrix.
        """
        return self.tocsr(nzmax).tocsc()



# symmetric sparse skyline
# diagonal (banded) matrix
# ellpack-itpack generalized diagonal
# block sparse row
# modified compressed sparse row
# block sparse column
# modified compressed sparse column
# symmetric skyline
# nonsymmetric skyline
# jagged diagonal
# unsymmetric sparse skyline
# variable block row

def _isinstance(x, _class):
    ##
    # This makes scipy.sparse.sparse.csc_matrix == __main__.csc_matrix.
    c1 = ('%s' % x.__class__).split( '.' )
    c2 = ('%s' % _class).split( '.' )
    aux = c1[-1] == c2[-1]
    return isinstance(x, _class) or aux

def isspmatrix(x):
    return _isinstance(x, spmatrix)

issparse = isspmatrix

def isspmatrix_csr(x):
    return _isinstance(x, csr_matrix)

def isspmatrix_csc(x):
    return _isinstance(x, csc_matrix)

def isspmatrix_dok(x):
    return _isinstance(x, dok_matrix)

def isspmatrix_lil( x ):
    return _isinstance(x, lil_matrix)

def isspmatrix_coo( x ):
    return _isinstance(x, coo_matrix)

def isdense(x):
    return _isinstance(x, ndarray)

def isscalarlike(x):
    """Is x either a scalar, an array scalar, or a 0-dim array?"""
    return isscalar(x) or (isdense(x) and x.ndim == 0)

def isintlike(x):
    """Is x appropriate as an index into a sparse matrix? Returns True
    if it can be cast safely to a machine int.
    """
    try:
        if int(x) == x:
            return True
        else:
            return False
    except TypeError:
        return False

def isshape(x):
    """Is x a valid 2-tuple of dimensions?
    """
    try:
        # Assume it's a tuple of matrix dimensions (M, N)
        (M, N) = x
        assert isintlike(M) and isintlike(N)   # raises TypeError unless integers
        assert M > 0 and N > 0
    except (ValueError, TypeError, AssertionError):
        return False
    else:
        return True

def getdtype(dtype, a=None, default=None):
    """Function used to simplify argument processing.  If 'dtype' is not
    specified (is None), returns a.dtype; otherwise returns a numpy.dtype
    object created from the specified dtype argument.  If 'dtype' and 'a'
    are both None, construct a data type out of the 'default' parameter.
    Furthermore, 'dtype' must be in 'allowed' set.
    """
    canCast = True
    if dtype is None:
        try:
            newdtype = a.dtype
        except AttributeError:
            if default is not None:
                newdtype = numpy.dtype(default)
                canCast = False
            else:
                raise TypeError, "could not interpret data type"
    else:
        newdtype = numpy.dtype(dtype)

    return newdtype


def _spdiags_tosub(diag_num, a, b):
    part1 = where(less(diag_num, a), abs(diag_num-a), 0)
    part2 = where(greater(diag_num, b), abs(diag_num-b), 0)
    return part1+part2

# Note: sparsetools only offers diagonal -> CSC matrix conversion functions,
# not to CSR
def spdiags(diags, offsets, M, N):
    """Return a sparse matrix in CSC format given its diagonals.

    B = spdiags(diags, offsets, M, N)

    Inputs:
        diags  --  rows contain diagonal values
        offsets -- diagonals to set (0 is main)
        M, N    -- sparse matrix returned is M X N
    """
    #    diags = array(transpose(diags), copy=True)
    diags = numpy.array(diags, copy = True)
    if diags.dtype.char not in 'fdFD':
        diags = diags.astype('d')
    if not hasattr(offsets, '__len__' ):
        offsets = (offsets,)
    offsets = array(offsets, copy=False, dtype=numpy.intc)
    assert(len(offsets) == diags.shape[0])
    indptr, rowind, data = sparsetools.spdiags(M, N, len(offsets), offsets, diags)
    return csc_matrix((data, rowind, indptr), (M, N))

def extract_diagonal(A):
    """
    extract_diagonal(A) returns the main diagonal of A.
    """
    if isspmatrix_csr(A):
        return sparsetools.extract_csr_diagonal(A.shape[0],A.shape[1],
                                                A.indptr,A.indices,A.data)
    elif isspmatrix_csc(A):
        return sparsetools.extract_csc_diagonal(A.shape[0],A.shape[1],
                                                A.indptr,A.indices,A.data)
    elif isspmatrix(A):
        return extract_diagonal(A.tocsr())
    else:
        raise ValueError,'expected sparse matrix'




def spidentity(n, dtype='d'):
    """
    spidentity( n ) returns the identity matrix of shape (n, n) stored
    in CSC sparse matrix format.
    """
    return csc_matrix((ones(n,dtype=dtype),arange(n),arange(n+1)),(n,n))


def speye(n, m, k = 0, dtype = 'd'):
    """
    speye(n, m) returns a (n x m) matrix stored
    in CSC sparse matrix format, where the  k-th diagonal is all ones,
    and everything else is zeros.
    """
    diags = ones((1, n), dtype = dtype)
    return spdiags(diags, k, n, m)

def spkron(a,b):
    """kronecker product of sparse matrices a and b

    *Parameters*:
        a,b : sparse matrices
            E.g. csr_matrix, csc_matrix, coo_matrix, etc.

    *Returns*:
        coo_matrix
            kronecker product in COOrdinate format

    *Example*:
    -------

    >>> a = csr_matrix(array([[0,2],[5,0]]))
    >>> b = csr_matrix(array([[1,2],[3,4]]))
    >>> spkron(a,b).todense()
    matrix([[  0.,   0.,   2.,   4.],
            [  0.,   0.,   6.,   8.],
            [  5.,  10.,   0.,   0.],
            [ 15.,  20.,   0.,   0.]])

    """
    if not isspmatrix(a) and isspmatrix(b):
        raise ValueError,'expected sparse matrix'

    a,b = a.tocoo(),b.tocoo()
    output_shape = (a.shape[0]*b.shape[0],a.shape[1]*b.shape[1])

    if a.nnz == 0 or b.nnz == 0:
        # kronecker product is the zero matrix
        return coo_matrix(None, dims=output_shape)


    # expand entries of a into blocks
    row  = a.row.repeat(b.nnz)
    col  = a.col.repeat(b.nnz)
    data = a.data.repeat(b.nnz)

    row *= b.shape[0]
    col *= b.shape[1]

    # increment block indices
    row,col = row.reshape(-1,b.nnz),col.reshape(-1,b.nnz)
    row += b.row
    col += b.col
    row,col = row.reshape(-1),col.reshape(-1)

    # compute block entries
    data = data.reshape(-1,b.nnz)
    data *= b.data
    data = data.reshape(-1)

    return coo_matrix((data,(row,col)), dims=output_shape)




def lil_eye((r,c), k=0, dtype=float):
    """Generate a lil_matrix of dimensions (r,c) with the k-th
    diagonal set to 1.

    :Parameters:
        r,c : int
            Row and column-dimensions of the output.
        k : int
            Diagonal offset.  In the output matrix,
            out[m,m+k] == 1 for all m.
        dtype : dtype
            Data-type of the output array.

    """
    out = lil_matrix((r,c),dtype=dtype)
    for c in xrange(clip(k,0,c),clip(r+k,0,c)):
        out.rows[c-k].append(c)
        out.data[c-k].append(1)
    return out

def lil_diags(diags,offsets,(m,n),dtype=float):
    """Generate a lil_matrix with the given diagonals.

    :Parameters:
        diags : list of list of values e.g. [[1,2,3],[4,5]]
            Values to be placed on each indicated diagonal.
        offsets : list of ints
            Diagonal offsets.  This indicates the diagonal on which
            the given values should be placed.
        (r,c) : tuple of ints
            Row and column dimensions of the output.
        dtype : dtype
           Output data-type.

    Example:
    -------

    >>> lil_diags([[1,2,3],[4,5],[6]],[0,1,2],(3,3)).todense()
    matrix([[ 1.,  4.,  6.],
            [ 0.,  2.,  5.],
            [ 0.,  0.,  3.]])

    """
    offsets_unsorted = list(offsets)
    diags_unsorted = list(diags)
    if len(diags) != len(offsets):
        raise ValueError("Number of diagonals provided should "
                         "agree with offsets.")

    sort_indices = numpy.argsort(offsets_unsorted)
    diags = [diags_unsorted[k] for k in sort_indices]
    offsets = [offsets_unsorted[k] for k in sort_indices]

    for i,k in enumerate(offsets):
        if len(diags[i]) < m-abs(k):
            raise ValueError("Not enough values specified to fill "
                             "diagonal %s." % k)

    out = lil_matrix((m,n),dtype=dtype)
    for k,diag in itertools.izip(offsets,diags):
        for ix,c in enumerate(xrange(clip(k,0,n),clip(m+k,0,n))):
            out.rows[c-k].append(c)
            out.data[c-k].append(diag[ix])
    return out

def issequence(t):
    return isinstance(t, (list, tuple))\
           or (isinstance(t, ndarray) and (t.ndim == 1))
