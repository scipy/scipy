"""Base class for sparse matrices"""

__all__ = ['spmatrix','isspmatrix','issparse']

from warnings import warn

from numpy import asarray, asmatrix, asanyarray, ones

from sputils import isdense, isscalarlike 


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


MAXPRINT = 50

class spmatrix(object):
    """ This class provides a base class for all sparse matrices.  It
    cannot be instantiated.  Most of the work is provided by subclasses.
    """

    __array_priority__ = 10.1
    ndim = 2
    def __init__(self, maxprint=MAXPRINT):
        self.format = self.__class__.__name__[:3]
        self._shape = None
        if self.format == 'spm':
            raise ValueError, "This class is not intended" \
                  " to be instantiated directly."
        self.maxprint = maxprint

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
        return self.tocsr().astype(t)

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
        """Return this matrix in a given sparse format

        *Parameters*:
            format : desired sparse matrix format
                If format is None then no conversion is performed
                Other possible values include:
                    "csr" for csr_matrix format
                    "csc" for csc_matrix format
                    "dok" for dok_matrix format and so on
        """

        if format is None or format == self.format:
            return self
        else:
            return getattr(self,'to' + format)()

    # default operations use the CSR format as a base
    #   and operations return in csr format
    #  thus, a new sparse matrix format just needs to define
    #  a tocsr method

    def __abs__(self):
        return abs(self.tocsr())

    def __add__(self, other):   # self + other
        return self.tocsr().__add__(other)

    def __radd__(self, other):  # other + self
        return self.tocsr().__radd__(other)

    def __sub__(self, other):   # self - other
        #note: this can't be replaced by self + (-other) for unsigned types
        return self.tocsr().__sub__(other)

    def __rsub__(self, other):  # other - self
        return self.tocsr().__rsub__(other)
    
    def multiply(self, other):
        """Point-wise multiplication by another matrix
        """
        return self.tocsr().multiply(other)

    def __mul__(self, other):
        return self.tocsr().__mul__(other)

    def __rmul__(self, other):
        return self.tocsr().__rmul__(other)

    def __truediv__(self, other):
        if isscalarlike(other):
            return self * (1./other)
        else:
            return self.tocsr().__truediv__(other)

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
        return self.tocsr().transpose()

    def conj(self):
        return self.tocsr().conj()

    def conjugate(self):
        return self.tocsr().conj()

    # Renamed conjtranspose() -> getH() for compatibility with dense matrices
    def getH(self):
        return self.transpose().conj()

    def _real(self):
        return self.tocsr()._real()

    def _imag(self):
        return self.tocsr()._imag()

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
            other = asanyarray(other)

        if isdense(other) and asarray(other).squeeze().ndim <= 1:
            # it's a dense row or column vector
            return self.matvec(other)
        elif len(other.shape) == 2:
            # it's a 2d dense array, dense matrix, or sparse matrix
            return self.matmat(other)
        else:
            raise ValueError, "could not interpret dimensions"


    def matmat(self, other):
        return self.tocsr().matmat(other)

    def matvec(self, other):
        """Multiplies the sparse matrix by the vector 'other', returning a
        dense vector as a result.
        """
        return self.tocsr().matvec(other)

    def rmatvec(self, other, conjugate=True):
        """Multiplies the vector 'other' by the sparse matrix, returning a
        dense vector as a result.

        If 'conjugate' is True:
            returns A.transpose().conj() * other
        Otherwise:
            returns A.transpose() * other.
        """
        return self.tocsr().rmatvec(other, conjugate=conjugate)

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
        return self.tocsr().toarray()

    def todok(self):
        return self.tocoo().todok()
    
    def tocoo(self):
        return self.tocsr().tocoo()

    def tolil(self):
        return self.tocsr().tolil()

    def todia(self):
        return self.tocoo().todia()

    def copy(self):
        return self.__class__(self,copy=True)

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
        #deprecated on Dec 14 2007
        warn('save() is deprecated, consider using mmwrite() or savemat()' \
                ' provided by scipy.io instead',
                DeprecationWarning)
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


from sputils import _isinstance

def isspmatrix(x):
    return _isinstance(x, spmatrix)

issparse = isspmatrix

