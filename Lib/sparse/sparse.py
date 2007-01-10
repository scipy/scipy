""" Scipy 2D sparse matrix module.

Original code by Travis Oliphant.
Modified and extended by Ed Schofield and Robert Cimrman.
Revision of sparsetools by Nathan Bell
"""

from numpy import zeros, isscalar, real, imag, asarray, asmatrix, matrix, \
                  ndarray, amax, rank, conj, searchsorted, ndarray,   \
                  less, where, greater, array, transpose, empty, ones, \
                  arange, shape, intc
import numpy
from scipy.sparse.sparsetools import densetocsr, csrtocsc, csrtodense, cscplcsc, \
     cscelmulcsc, cscmux, csrmux, csrmucsr, csrtocoo, cootocsc, cootocsr, \
     cscmucsc, csctocoo, csctocsr, csrplcsr, csrelmulcsr
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



_coerce_rules = {('f', 'f'):'f', ('f', 'd'):'d', ('f', 'F'):'F',
                 ('f', 'D'):'D', ('d', 'f'):'d', ('d', 'd'):'d',
                 ('d', 'F'):'D', ('d', 'D'):'D', ('F', 'f'):'F',
                 ('F', 'd'):'D', ('F', 'F'):'F', ('F', 'D'):'D',
                 ('D', 'f'):'D', ('D', 'd'):'d', ('D', 'F'):'D',
                 ('D', 'D'):'D'}
_transtabl = {'f':'s', 'd':'d', 'F':'c', 'D':'z'}
_itranstabl = {'s':'f', 'd':'d', 'c':'F', 'z':'D'}
def _convert_data(data1, data2, newtype):
    if data1.dtype.char != newtype:
        data1 = data1.astype(newtype)
    if data2.dtype.char != newtype:
        data2 = data2.astype(newtype)
    return data1, data2




class spmatrix(object):
    """ This class provides a base class for all sparse matrices.  It
    cannot be instantiated.  Most of the work is provided by subclasses.
    """

    __array_priority__ = 10.1
    ndim = 2
    def __init__(self, maxprint=MAXPRINT, allocsize=ALLOCSIZE):
        self.format = self.__class__.__name__[:3]
        if self.format == 'spm':
            raise ValueError, "This class is not intended" \
                  " to be instantiated directly."
        self.maxprint = maxprint
        self.allocsize = allocsize

    def astype(self, t):
        csc = self.tocsc()
        return csc.astype(t)

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
        val = ''
        for ind in xrange(start, stop):
            val = val + '  %s\t%s\n' % (self.rowcol(ind), self.getdata(ind))
        return val

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

    def __cmp__(self, other):
        raise NotImplementedError, "comparison of sparse matrices not implemented"

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

    def __add__(self, other):
        csc = self.tocsc()
        return csc + other

    def __radd__(self, other):  # other + self
        csc = self.tocsc()
        return csc.__radd__(other)

    def __sub__(self, other):   # self - other
        neg_other = -other
        return self + neg_other

    def __rsub__(self, other):  # other - self
        neg_self = -self
        return other + neg_self

    def __mul__(self, other):
        csc = self.tocsc()
        return csc * other

    def __truediv__(self, other):
        if isscalarlike(other):
            return self * (1./other)
        else:
            raise NotImplementedError, "sparse matrix division not yet supported"

    def __div__(self, other):
        # Always do true division
        if isscalarlike(other):
            return self * (1./other)
        else:
            raise NotImplementedError, "sparse matrix division not yet supported"

    def __pow__(self, other):
        csc = self.tocsc()
        return csc ** other

    def __rmul__(self, other):
        csc = self.tocsc()
        return csc.__rmul__(other)

    def __neg__(self):
        csc = self.tocsc()
        return -csc

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
        csc.data = real(csc.data)
        csc.dtype = csc.data.dtype
        csc.ftype = _transtabl[csc.dtype.char]
        return csc

    def _imag(self):
        csc = self.tocsc()
        csc.data = imag(csc.data)
        csc.dtype = csc.data.dtype
        csc.ftype = _transtabl[csc.dtype.char]
        return csc

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
        multiplication.  Returns A.transpose().conj() * other or
        A.transpose() * other.
        """

        try:
            other.shape
        except AttributeError:
            # If it's a list or whatever, treat it like a matrix
            other = asmatrix(other)

        if len(other.shape) == 1:
            result = self.matvec(other)
        elif isdense(other) and asarray(other).squeeze().ndim == 1:
            # If it's a row or column vector, return a DENSE result
            result = self.matvec(other)
        elif len(other.shape) == 2:
            # Return a sparse result
            result = self.matmat(other)
        else:
            raise ValueError, "could not interpret dimensions"

        if isinstance(other, matrix) and isdense(result):
            return asmatrix(result)
        else:
            # if the result is sparse or 'other' is an array:
            return result

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
            # Does the following multiplication work in NumPy now?
            o = asmatrix(ones((1, m), dtype=self.dtype))
            return o * self
            # o = ones(m, dtype=self.dtype)
            # return asmatrix(self.rmatvec(o))
        elif axis == 1:
            # sum over rows
            o = asmatrix(ones((n, 1), dtype=self.dtype))
            return self * o
        elif axis == None:
            # sum over rows and columns
            m, n = self.shape
            o0 = asmatrix(ones((1, m), dtype=self.dtype))
            o1 = asmatrix(ones((n, 1), dtype=self.dtype))
            return (o0 * (self * o1)).A.squeeze()
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

    def astype(self, t):
        out = self.copy()
        out.data = out.data.astype(t)
        out.dtype = out.data.dtype
        out.ftype = _transtabl[out.dtype.char]
        return out

    def __repr__(self):
        format = self.getformat()
        return "<%dx%d sparse matrix of type '%s'\n\twith %d stored "\
               "elements (space for %d)\n\tin %s format>" % \
               (self.shape + (self.dtype.type, self.getnnz(), self.nzmax, \
                   _formats[format][1]))



    def __add__(self, other, self_ind, other_ind, fn, cls):
        # First check if argument is a scalar
        if isscalarlike(other):
            # Now we would add this scalar to every element.
            raise NotImplementedError, 'adding a scalar to a CSC or CSR ' \
                  'matrix is not yet supported'
        elif isspmatrix(other):
            other = other.tocsc()
            if (other.shape != self.shape):
                raise ValueError, "inconsistent shapes"
            indptr, ind, data = fn(self.shape[0], self.shape[1], \
                                         self.indptr, self_ind, \
                                         self.data, other.indptr, \
                                         other_ind, other.data)
            return cls((data, ind, indptr), self.shape)
        elif isdense(other):
            # Convert this matrix to a dense matrix and add them
            return other + self.todense()
        else:
            raise TypeError, "unsupported type for sparse matrix addition"


    def __mul__(self, other):
        """ Scalar, vector, or matrix multiplication
        """
        if isscalarlike(other):
            new = self.copy()
            new.data = other * new.data         # allows type conversion
            new.dtype = new.data.dtype
            new.ftype = _transtabl[new.dtype.char]
            return new
        else:
            return self.dot(other)


    def __rmul__(self, other):  # other * self
        if isscalarlike(other):
            new = self.copy()
            new.data = other * new.data         # allows type conversion
            new.dtype = new.data.dtype
            new.ftype = _transtabl[new.dtype.char]
            return new
        else:
            # Don't use asarray unless we have to
            try:
                tr = other.transpose()
            except AttributeError:
                tr = asarray(other).transpose()
            return self.transpose().dot(tr).transpose()


    def __neg__(self):
        new = self.copy()
        new.data *= -1
        return new

    def __pow__(self, other, self_ind, other_ind, fn, cls):
        """ Element-by-element power (unless other is a scalar, in which
        case return the matrix power.)
        """
        if isscalarlike(other):
            new = self.copy()
            new.data = new.data ** other
            new.dtype = new.data.dtype
            new.ftype = _transtabl[new.dtype.char]
            return new
        elif isspmatrix(other):
            other = other.tocsr()
            if (other.shape != self.shape):
                raise ValueError, "inconsistent shapes"
            indptr, ind, data = fn(self.shape[0], self.shape[1], \
                                               self.indptr, self_ind, \
                                               self.data, other.indptr, \
                                               other_ind, other.data)
            return cls((data, ind, indptr), (self.shape[0], other.shape[1]))
        else:
            raise TypeError, "unsupported type for sparse matrix power"


    def _matmat(self, other, self_ind, other_ind, fn, cls):
        if isspmatrix(other):
            M, K1 = self.shape
            K2, N = other.shape
            if (K1 != K2):
                raise ValueError, "shape mismatch error"
            other = other.tocsc()
            indptr, ind, data = fn(M, N, self.indptr, self_ind, \
                                   self.data, other.indptr, \
                                   other_ind, other.data)
            return cls((data, ind, indptr), (M, N))      
        elif isdense(other):
            # This is SLOW!  We need a more efficient implementation
            # of sparse * dense matrix multiplication!
            return self.matmat(csc_matrix(other))
        else:
            raise TypeError, "need a dense or sparse matrix"

    def _matvec(self, other, self_ind, fn):
        if isdense(other):
            # This check is too harsh -- it prevents a column vector from
            # being created on-the-fly like dense matrix objects can.
            #if len(other) != self.shape[1]:
            #    raise ValueError, "dimension mismatch"
            oth = numpy.ravel(other)            
            y = fn(self.shape[0], self.shape[1], \
                   self.indptr, self_ind, self.data, oth)
            if isinstance(other, matrix):
                y = asmatrix(y)
                # If 'other' was an (nx1) column vector, transpose the result
                # to obtain an (mx1) column vector.
                if other.ndim == 2 and other.shape[1] == 1:
                    y = y.T
            return y
        elif isspmatrix(other):
            raise TypeError, "use matmat() for sparse * sparse"
        else:
            raise TypeError, "need a dense vector"

    def _rmatvec(self, other, shape0, shape1, fn, conjugate=True):
        if isdense(other):
            # This check is too harsh -- it prevents a column vector from
            # being created on-the-fly like dense matrix objects can.
            # if len(other) != self.shape[0]:
            #    raise ValueError, "dimension mismatch"
            if conjugate:
                cd = conj(self.data)
            else:
                cd = self.data
            oth = numpy.ravel(other)            
            y = fn(shape0, shape1, self.indptr, self.rowind, cd, oth)
            if isinstance(other, matrix):
                y = asmatrix(y)
                # In the (unlikely) event that this matrix is 1x1 and 'other'
                # was an (mx1) column vector, transpose the result.
                if other.ndim == 2 and other.shape[1] == 1:
                    y = y.T
            return y
        elif isspmatrix(other):
            raise TypeError, "use matmat() for sparse * sparse"
        else:
            raise TypeError, "need a dense vector"

    def getdata(self, ind):
        return self.data[ind]

    def _tocoo(self, fn, self_ind):
        rows, cols, data = fn(self.shape[0], self.shape[1], \
                              self.indptr, self_ind, self.data)
        return coo_matrix((data, (rows, cols)), self.shape)



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
                a[ij[k, 0], ij[k, 1]] = data[k]

          - csc_matrix((data, row, ptr), [(M, N)])
            standard CSC representation
    """
    def __init__(self, arg1, dims=None, nzmax=NZMAX, dtype=None, copy=False):
        _cs_matrix.__init__(self)
        if isdense(arg1):
            self.dtype = getdtype(dtype, arg1)
            # Convert the dense array or matrix arg1 to CSC format
            if rank(arg1) == 1:
                # Convert to a row vector
                arg1 = arg1.reshape(1, arg1.shape[0])
            if rank(arg1) == 2:
                #s = asarray(arg1)
                s = arg1
                if s.dtype.char not in 'fdFD':
                    # Use a double array as the source (but leave it alone)
                    s = s*1.0
                if (rank(s) == 2):
                    self.shape = s.shape
                    self.indptr, self.rowind, self.data = densetocsr(s.shape[1], \
                                                                     s.shape[0], \
                                                                     s.T)
            else:
                raise ValueError, "dense array must have rank 1 or 2"
        elif isspmatrix(arg1):
            s = arg1
            self.dtype = getdtype(dtype, s)
            if isinstance(s, csc_matrix):
                # do nothing but copy information
                self.shape = s.shape
                if copy:
                    self.data = s.data.copy()
                    self.rowind = s.rowind.copy()
                    self.indptr = s.indptr.copy()
                else:
                    self.data = s.data
                    self.rowind = s.rowind
                    self.indptr = s.indptr
            elif isinstance(s, csr_matrix):
                self.shape = s.shape
                self.indptr, self.rowind, self.data = csrtocsc(s.shape[0],
                                                               s.shape[1],
                                                               s.indptr,
                                                               s.colind,
                                                               s.data)
            else:
                temp = s.tocsc()
                self.data = temp.data
                self.rowind = temp.rowind
                self.indptr = temp.indptr
                self.shape = temp.shape
        elif type(arg1) == tuple:
            if isshape(arg1):
                self.dtype = getdtype(dtype, default=float)
                # It's a tuple of matrix dimensions (M, N)
                M, N = arg1
                self.data = zeros((nzmax,), self.dtype)
                self.rowind = zeros((nzmax,), intc)
                self.indptr = zeros((N+1,), intc)
                self.shape = (M, N)
            else:
                try:
                    # Try interpreting it as (data, ij)
                    (s, ij) = arg1
                    assert isinstance(ij, ndarray) and (rank(ij) == 2) \
                            and (shape(ij) == (2, len(s)))
                except (AssertionError, TypeError, ValueError):
                    try:
                        # Try interpreting it as (data, rowind, indptr)
                        (s, rowind, indptr) = arg1
                        self.dtype = getdtype(dtype, s)
                        if copy:
                            self.data = array(s)
                            self.rowind = array(rowind)
                            self.indptr = array(indptr, dtype=intc)
                        else:
                            self.data = asarray(s)
                            self.rowind = asarray(rowind)
                            self.indptr = asarray(indptr, dtype=intc)
                    except:
                        raise ValueError, "unrecognized form for csc_matrix constructor"
                else:
                    # (data, ij) format
                    self.dtype = getdtype(dtype, s)
                    ijnew = array(ij, copy=copy)
                    temp = coo_matrix((s, ijnew), dims=dims, \
                                      dtype=self.dtype).tocsc()
                    self.shape = temp.shape
                    self.data = temp.data
                    self.rowind = temp.rowind
                    self.indptr = temp.indptr
        else:
            raise ValueError, "unrecognized form for csc_matrix constructor"

        # Read existing matrix dimensions
        try:
            (oldM, oldN) = self.shape
        except:
            oldM = oldN = None
        # Read matrix dimensions given, if any
        if dims is not None:
            try:
                (M, N) = dims
            except (TypeError, ValueError), e:
                raise TypeError, "dimensions not understood"
        else:
            M = N = None
        if len(self.rowind) > 0:
            M = max(oldM, M, int(amax(self.rowind)) + 1)
        else:
            # Matrix is completely empty
            M = max(oldM, M)
        N = max(0, oldN, N, len(self.indptr) - 1)
        self.shape = (M, N)
        self._check()

    def _check(self):
        # some functions pass floats
        self.shape = tuple([int(x) for x in self.shape])

        M, N = self.shape
        nnz = self.indptr[-1]
        nzmax = len(self.rowind)
        if (rank(self.data) != 1) or (rank(self.rowind) != 1) or \
           (rank(self.indptr) != 1):
            raise ValueError, "data, rowind, and indptr arrays "\
                  "should be rank 1"
        if (len(self.data) != nzmax):
            raise ValueError, "data and row list should have same length"
        if (len(self.indptr) != N+1):
            raise ValueError, "index pointer should be of of size N+1"
        if (nzmax < nnz):
            raise ValueError, "nzmax must not be less than nnz"
        if (nnz>0) and (amax(self.rowind[:nnz]) >= M):
            raise ValueError, "row values must be < M"
        if (self.indptr[-1] > len(self.rowind)):
            raise ValueError, \
                  "Last value of index list should be less than "\
                  "the size of data list"
        if (self.rowind.dtype != numpy.intc):
            self.rowind = self.rowind.astype(numpy.intc)
            self.indptr = self.indptr.astype(numpy.intc)
        self.nnz = nnz
        self.nzmax = nzmax
        self.dtype = self.data.dtype
        if self.dtype.char not in 'fdFD':
            self.data = 1.0 * self.data
            self.dtype = self.data.dtype

        self.ftype = _transtabl[self.dtype.char]

    
    def __radd__(self, other):
        """ Function supporting the operation: self + other.
        """
        if isscalarlike(other):
            raise NotImplementedError, 'adding a scalar to a CSC matrix is ' \
                    'not yet supported'
        elif isspmatrix(other):
            ocs = other.tocsc()
            if (ocs.shape != self.shape):
                raise ValueError, "inconsistent shapes"
            indptr, rowind, data = cscplcsc(self.shape[0], self.shape[1], \
                                            self.indptr, self.rowind, \
                                            self.data, ocs.indptr, \
                                            ocs.rowind, ocs.data)
            return csc_matrix((data, rowind, indptr), self.shape)
        elif isdense(other):
            # Convert this matrix to a dense matrix and add them.
            return self.todense() + other
        else:
            raise TypeError, "unsupported type for sparse matrix addition"

    def __add__(self, other):
        _cs_matrix.__add__(self, other, self.rowind, other.rowind, cscplcsc, csc_matrix)

    def __pow__(self, other):
        _cs_matrix.__pow__(self, other, self.rowind, other.rowind, cscelmulcsc, csc_matrix)

    def transpose(self, copy=False):
        M, N = self.shape
        new = csr_matrix((N, M), nzmax=self.nzmax, dtype=self.dtype)
        if copy:
            new.data = self.data.copy()
            new.colind = self.rowind.copy()
            new.indptr = self.indptr.copy()
        else:
            new.data = self.data
            new.colind = self.rowind
            new.indptr = self.indptr
        new._check()
        return new

    def conj(self, copy=False):
        new = csc_matrix(self.shape, nzmax=self.nzmax, dtype=self.dtype)
        if copy:
            new.data = self.data.conj().copy()
            new.rowind = self.rowind.conj().copy()
            new.indptr = self.indptr.conj().copy()
        else:
            new.data = self.data.conj()
            new.rowind = self.rowind.conj()
            new.indptr = self.indptr.conj()
        new._check()
        return new

    def sum(self, axis=None):
        # Override the base class sum method for efficiency in the cases
        # axis=0 and axis=None.
        m, n = self.shape
        data = self.data
        if axis in (0, None):
            out = empty(n, dtype=self.dtype)
            # The first element in column j has index indptr[j], the last
            # indptr[j+1]
            indptr = self.indptr
            for j in xrange(n):
                out[j] = data[indptr[j] : indptr[j+1]].sum()
            if axis == 0:
                # Output is a (1 x n) dense matrix
                return asmatrix(out)
            else:
                return out.sum()
        else:
            rowind = self.rowind
            out = zeros(m, dtype=self.dtype)
            # Loop over non-zeros
            for k in xrange(self.nnz):
                out[rowind[k]] += data[k]
            # Output is a (m x 1) dense matrix
            return asmatrix(out).T
                    
    def matvec(self, other):
        _cs_matrix._matvec(self, other, self.rowind, cscmux)

    def rmatvec(self, other, conjugate=True):
        _cs_matrix._rmatvec(self, other, shape[1], shape[0], cscmux, conjugate=conjugate)

    def matmat(self, other):
        _cs_matrix._matmat(self, other, self.rowind, other.rowind, cscmucsc, csc_matrix)


    def __getitem__(self, key):
        if isinstance(key, tuple):
            row = key[0]
            col = key[1]
            if isinstance(col, slice):
                raise IndexError, "csc_matrix supports slices only of a single"\
                                  " column"
            elif isinstance(row, slice):
                return self._getcolslice(row, col)
            M, N = self.shape
            if (row < 0):
                row = M + row
            if (col < 0):
                col = N + col
            if not (0<=row<M) or not (0<=col<N):
                raise IndexError, "index out of bounds"
            #this was implemented in fortran before - is there a noticable performance advangate?
            indxs = numpy.where(row == self.rowind[self.indptr[col]:self.indptr[col+1]])
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

            indxs = numpy.where(row == self.rowind[self.indptr[col]:self.indptr[col+1]])
            if len(indxs[0]) == 0:
                #value not present
                nzmax = self.nzmax
                if (nzmax < self.nnz+1):  # need more room
                    alloc = max(1, self.allocsize)
                    self.data = resize1d(self.data, nzmax + alloc)
                    self.rowind = resize1d(self.rowind, nzmax + alloc)
                    
                newindex = self.indptr[col]
                self.data[newindex+1:]   = self.data[newindex:-1]
                self.rowind[newindex+1:] = self.rowind[newindex:-1]
                
                self.data[newindex]   = val
                self.rowind[newindex] = row
                self.indptr[col+1:] += 1
                
            elif len(indxs[0]) == 1:
                #value already present
                self.data[self.indptr[col]:self.indptr[col+1]][indxs[0]] = val
            else:
                raise IndexError, "row index occurs more than once"
            
            self._check()
        else:
            # We should allow slices here!
            raise IndexError, "invalid index"

    def _getcolslice(self, myslice, j):
        """Returns a view of the elements [myslice.start:myslice.stop, j].
        """
        start, stop, stride = myslice.indices(self.shape[0])
        if stride != 1:
            raise ValueError, "slicing with step != 1 not supported"
        if stop <= start:
            raise ValueError, "slice width must be >= 1"

        indices = []
        
        for ind in xrange(self.indptr[j], self.indptr[j+1]):
            if self.rowind[ind] >= start and self.rowind[ind] < stop:
                indices.append(ind)

        rowind = self.rowind[indices] - start
        data   = self.data[indices]
        indptr = numpy.array([0, len(indices)])
        return csc_matrix((data, rowind, indptr), dims=(stop-start, 1), \
                          dtype=self.dtype)
    
    def rowcol(self, ind):
        row = self.rowind[ind]
        col = searchsorted(self.indptr, ind+1)-1
        return (row, col)

    def tocsc(self, copy=False):
        return self.toself(copy)

    def tocoo(self):
        _cs_matrix._tocoo(self, csctocoo, self.rowind)

    def tocsr(self):
        indptr, colind, data = csctocsr(self.shape[0], self.shape[1], \
                                        self.indptr, self.rowind, self.data)
        return csr_matrix((data, colind, indptr), self.shape)
    
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
        self.rowind = self.rowind[:nnz]
        self.nzmax = nnz
        self._check()

    def ensure_sorted_indices(self, inplace=False):
        """Return a copy of this matrix where the row indices are sorted
        """
        if inplace:
            sparsetools.ensure_sorted_indices(self.shape[1], self.shape[0],
                                              self.indptr, self.rowind,
                                              self.data )
        else:
            return self.tocsr().tocsc()


    def copy(self):
        new = csc_matrix(self.shape, nzmax=self.nzmax, dtype=self.dtype)
        new.data = self.data.copy()
        new.rowind = self.rowind.copy()
        new.indptr = self.indptr.copy()
        new._check()
        return new


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
            # Convert the dense array or matrix arg1 to CSR format
            if rank(arg1) == 1:
                # Convert to a row vector
                arg1 = arg1.reshape(1, arg1.shape[0])
            if rank(arg1) == 2:
                s = arg1
                ocsc = csc_matrix(transpose(s))
                self.colind = ocsc.rowind
                self.indptr = ocsc.indptr
                self.data = ocsc.data
                self.shape = (ocsc.shape[1], ocsc.shape[0])
            else:
                raise ValueError, "dense array must have rank 1 or 2"
        elif isspmatrix(arg1):
            s = arg1
            self.dtype = getdtype(dtype, s)
            if isinstance(s, csr_matrix):
                # do nothing but copy information
                self.shape = s.shape
                if copy:
                    self.data = s.data.copy()
                    self.colind = s.colind.copy()
                    self.indptr = s.indptr.copy()
                else:
                    self.data = s.data
                    self.colind = s.colind
                    self.indptr = s.indptr
            else:
                try:
                    temp = s.tocsr()
                except AttributeError:
                    temp = csr_matrix(s.tocsc())
                self.data = temp.data
                self.colind = temp.colind
                self.indptr = temp.indptr
                self.shape = temp.shape
        elif type(arg1) == tuple:
            if isshape(arg1):
                # It's a tuple of matrix dimensions (M, N)
                M, N = arg1
                self.dtype = getdtype(dtype, default=float)
                self.data = zeros((nzmax,), self.dtype)
                self.colind = zeros((nzmax,), intc)
                self.indptr = zeros((M+1,), intc)
                self.shape = (M, N)
            else:
                try:
                    # Try interpreting it as (data, ij)
                    (s, ij) = arg1
                    assert isinstance(ij, ndarray) and (rank(ij) == 2) \
                           and (shape(ij) == (2, len(s)))
                except (AssertionError, TypeError, ValueError, AttributeError):
                    try:
                        # Try interpreting it as (data, colind, indptr)
                        (s, colind, indptr) = arg1
                    except (TypeError, ValueError):
                        raise ValueError, "unrecognized form for csr_matrix constructor"
                    else:
                        self.dtype = getdtype(dtype, s)
                        if copy:
                            self.data = array(s, dtype=self.dtype)
                            self.colind = array(colind)
                            self.indptr = array(indptr, dtype=intc)
                        else:
                            self.data = asarray(s, dtype=self.dtype)
                            self.colind = asarray(colind)
                            self.indptr = asarray(indptr, dtype=intc)
                else:
                    # (data, ij) format
                    self.dtype = getdtype(dtype, s)
                    ijnew = array([ij[1], ij[0]], copy=copy)
                    temp = coo_matrix((s, ijnew), dims=dims, \
                                      dtype=self.dtype).tocsr()
                    self.shape = temp.shape
                    self.data = temp.data
                    self.colind = temp.colind
                    self.indptr = temp.indptr
        else:
            raise ValueError, "unrecognized form for csr_matrix constructor"

        # Read existing matrix dimensions
        try:
            (oldM, oldN) = self.shape
        except:
            oldM = oldN = None
        # Read matrix dimensions given, if any
        if dims is not None:
            try:
                (M, N) = dims
            except (TypeError, ValueError), e:
                raise TypeError, "dimensions not understood"
        else:
            M = N = None
        M = max(0, oldM, M, len(self.indptr) - 1)
        if len(self.colind) > 0:
            N = max(oldN, N, int(amax(self.colind)) + 1)
        else:
            # Matrix is completely empty
            N = max(oldN, N)
        self.shape = (M, N)
        self._check()

    def _check(self):
        # some functions pass floats
        self.shape = tuple([int(x) for x in self.shape])

        M, N = self.shape
        nnz = self.indptr[-1]
        nzmax = len(self.colind)
        if (rank(self.data) != 1) or (rank(self.colind) != 1) or \
           (rank(self.indptr) != 1):
            raise ValueError, "data, colind, and indptr arrays "\
                  "should be rank 1"
        if (len(self.data) != nzmax):
            raise ValueError, "data and row list should have same length"
        if (len(self.indptr) != M+1):
            raise ValueError, "index pointer should be of length #rows + 1"
        if (nnz>0) and (amax(self.colind[:nnz]) >= N):
            raise ValueError, "column-values must be < N"
        if (nnz > nzmax):
            raise ValueError, \
                  "last value of index list should be less than "\
                  "the size of data list"
        if (self.colind.dtype != numpy.intc):
            self.colind = self.colind.astype(numpy.intc)
            self.indptr = self.indptr.astype(numpy.intc)
        self.nnz = nnz
        self.nzmax = nzmax
        self.dtype = self.data.dtype
        if self.dtype.char not in 'fdFD':
            self.data = self.data + 0.0
            self.dtype = self.data.dtype

        self.ftype = _transtabl[self.dtype.char]

    
    def __add__(self, other):
        _cs_matrix.__add__(self, other, self.colind, other.colind, csrplcsr, csr_matrix)


    def __pow__(self, other):
        _cs_matrix.__pow__(self, other, self.colind, other.colind, csrelmulcsr, csr_matrix)

    def transpose(self, copy=False):
        M, N = self.shape
        new = csc_matrix((N, M), nzmax=self.nzmax, dtype=self.dtype)
        if copy:
            new.data = self.data.copy()
            new.rowind = self.colind.copy()
            new.indptr = self.indptr.copy()
        else:
            new.data = self.data
            new.rowind = self.colind
            new.indptr = self.indptr
        new._check()
        return new

    def sum(self, axis=None):
        # Override the base class sum method for efficiency in the cases
        # axis=1 and axis=None.
        m, n = self.shape
        data = self.data
        if axis in (1, None):
            out = empty(m, dtype=self.dtype)
            # The first element in row i has index indptr[i], the last
            # indptr[i+1]
            indptr = self.indptr
            for i in xrange(m):
                out[i] = data[indptr[i] : indptr[i+1]].sum()
            if axis == 1:
                # Output is a (m x 1) dense matrix
                return asmatrix(out).T
            else:
                return out.sum()
        else:
            colind = self.colind
            out = zeros(n, dtype=self.dtype)
            # Loop over non-zeros
            for k in xrange(self.nnz):
                out[colind[k]] += data[k]
            # Output is a (1 x n) dense matrix
            return asmatrix(out)

    def matvec(self, other):
        _cs_matrix._matvec(self, other, self.colind, csrmux)
        
    def rmatvec(self, other, conjugate=True):
        _cs_matrix._rmatvec(self, other, shape[0], shape[1], csrmux, conjugate=conjugate)

    def matmat(self, other):
        _cs_matrix._matmat(self, other, self.colind, other.colind, csrmucsr, csr_matrix)


    def __getitem__(self, key):
        if isinstance(key, tuple):
            row = key[0]
            col = key[1]
            if isinstance(row, slice):
                raise IndexError, "csr_matrix supports slices only of a single"\
                                  " row"
            elif isinstance(col, slice):
                return self._getrowslice(row, col)
            M, N = self.shape
            if (row < 0):
                row = M + row
            if (col < 0):
                col = N + col
            if not (0<=row<M) or not (0<=col<N):
                raise IndexError, "index out of bounds"

            indxs = numpy.where(col == self.colind[self.indptr[row]:self.indptr[row+1]])
            if len(indxs[0]) == 0:
                return 0
            else:
                return self.data[self.indptr[row]:self.indptr[row+1]][indxs[0]]
        elif isintlike(key):
            return self[key, :]
        else:
            raise IndexError, "invalid index"
    
    
    def _getrowslice(self, i, myslice):
        """Returns a view of the elements [i, myslice.start:myslice.stop].
        """
        start, stop, stride = myslice.indices(self.shape[1])
        if stride != 1:
            raise ValueError, "slicing with step != 1 not supported"
        if stop <= start:
            raise ValueError, "slice width must be >= 1"

        indices = []
        
        for ind in xrange(self.indptr[i], self.indptr[i+1]):
            if self.colind[ind] >= start and self.colind[ind] < stop:
                indices.append(ind)

        colind = self.colind[indices] - start
        data   = self.data[indices]
        indptr = numpy.array([0, len(indices)])
        return csr_matrix((data, colind, indptr), dims=(1, stop-start), \
                          dtype=self.dtype)
    
    def __setitem__(self, key, val):
        if isinstance(key, tuple):
            row = key[0]
            col = key[1]
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

            indxs = numpy.where(col == self.colind[self.indptr[row]:self.indptr[row+1]])
            if len(indxs[0]) == 0:
                #value not present
                nzmax = self.nzmax
                if (nzmax < self.nnz+1):  # need more room
                    alloc = max(1, self.allocsize)
                    self.data = resize1d(self.data, nzmax + alloc)
                    self.colind = resize1d(self.colind, nzmax + alloc)
                    
                newindex = self.indptr[row]
                self.data[newindex+1:]   = self.data[newindex:-1]
                self.colind[newindex+1:] = self.colind[newindex:-1]
                
                self.data[newindex]   = val
                self.colind[newindex] = col
                self.indptr[row+1:] += 1
                
            elif len(indxs[0]) == 1:
                #value already present
                self.data[self.indptr[row]:self.indptr[row+1]][indxs[0]] = val
            else:
                raise IndexError, "row index occurs more than once"
            
            self._check()
        else:
            # We should allow slices here!
            raise IndexError, "invalid index"

    def rowcol(self, ind):
        col = self.colind[ind]
        row = searchsorted(self.indptr, ind+1)-1
        return (row, col)

    def tocsr(self, copy=False):
        return self.toself(copy)

    def tocoo(self):
        _cs_matrix._tocoo(self, csrtocoo, self.colind)

    def tocsc(self):
        indptr, rowind, data = csrtocsc(self.shape[0], self.shape[1], \
                                        self.indptr, self.colind, self.data)
        return csc_matrix((data, rowind, indptr), self.shape)

    
    def toarray(self):
        data = numpy.zeros(self.shape, self.data.dtype)
        csrtodense(self.shape[0], self.shape[1], self.indptr, self.colind,
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
        self.colind = self.colind[:nnz]
        self.nzmax = nnz
        self._check()

    def ensure_sorted_indices(self, inplace=False):
        """Return a copy of this matrix where the column indices are sorted
        """
        if inplace:
            sparsetools.ensure_sorted_indices(self.shape[0], self.shape[1],
                                              self.indptr, self.colind,
                                              self.data )
        else:
            return self.tocsc().tocsr()


    def copy(self):
        new = csr_matrix(self.shape, nzmax=self.nzmax, dtype=self.dtype)
        new.data = self.data.copy()
        new.colind = self.colind.copy()
        new.indptr = self.indptr.copy()
        new._check()
        return new

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
        spmatrix.__init__(self)
        if shape is None:
            self.shape = (0, 0)
        else:
            try:
                m, n = shape
            except:
                raise "shape not understood"
            else:
                self.shape = shape
        self.dtype = getdtype(dtype, A, default=float)
        if A is not None:
            if type(A) == tuple:
                # Interpret as dimensions
                try:
                    dims = A
                    (M, N) = dims
                    assert M == int(M) and M > 0
                    assert N == int(N) and N > 0
                    self.shape = (int(M), int(N))
                    return
                except (TypeError, ValueError, AssertionError):
                    raise TypeError, "dimensions must be a 2-tuple of positive"\
                            " integers"
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
            assert len(key) == 2
        except (AssertionError, TypeError):
            raise TypeError, "index must be a pair of integers or slices"
        i, j = key

        # Bounds checking
        if isintlike(i):
            if i < 0 or i >= self.shape[0]:
                raise IndexError, "index out of bounds"
        if isintlike(j):
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
                    for (ii, jj) in self:
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
            for (ii, jj) in self:
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
            if i < 0 or i >= self.shape[0] or j < 0 or j >= self.shape[1]:
                raise IndexError, "index out of bounds"
            if isintlike(value) and value == 0:
                if key in self:  # get rid of it something already there
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
            for key in other:
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
        for key in self:
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
            for key in self:
                num = searchsorted(cols_or_rows, key[1])
                if num < N:
                    newkey = (key[0], num)
                    new[newkey] = self[key]
        else:
            for key in self:
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
            for key in self:
                num = searchsorted(cols_or_rows, key[1])
                if cols_or_rows[num] == key[1]:
                    newkey = (key[0], num)
                    ext[newkey] = self[key]
                else:
                    newkey = (key[0], key[1]-num)
                    base[newkey] = self[key]
        else:
            for key in self:
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
            for key in self:
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
            for key in self:
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
        keys = [(k[1], k[0]) for k in self]
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
        for key in self:
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
        M, N = self.shape
        try:
            newM, newN = shape
            assert newM == int(newM) and newM > 0
            assert newN == int(newN) and newN > 0
        except (TypeError, ValueError, AssertionError):
            raise TypeError, "dimensions must be a 2-tuple of positive"\
                             " integers"
        if newM < M or newN < N:
            # Remove all elements outside new dimensions
            for (i, j) in self.keys():
                if i >= newM or j >= newN:
                    del self[i, j]
        self.shape = (newM, newN)



class coo_matrix(spmatrix):
    """ A sparse matrix in coordinate list format.

    COO matrices are created either as:
        A = coo_matrix(None, dims=(m, n), [dtype])
    for a zero matrix, or as:
        A = coo_matrix((obj, ij), [dims])
    where the dimensions are optional.  If supplied, we set (M, N) = dims.
    If not supplied, we infer these from the index arrays
    ij[0][:] and ij[1][:]
    
    The arguments 'obj' and 'ij' represent three arrays:
        1. obj[:]    the entries of the matrix, in any order
        2. ij[0][:]  the row indices of the matrix entries
        3. ij[1][:]  the column indices of the matrix entries
    
    So the following holds:
        A[ij[0][k], ij[1][k] = obj[k]
    """
    def __init__(self, arg1, dims=None, dtype=None):
        spmatrix.__init__(self)
        if isinstance(arg1, tuple):
            try:
                obj, ij = arg1
            except:
                raise TypeError, "invalid input format"
        elif arg1 is None:      # clumsy!  We should make ALL arguments
                                # keyword arguments instead!
            # Initialize an empty matrix.
            if not isinstance(dims, tuple) or not isintlike(dims[0]):
                raise TypeError, "dimensions not understood"
            self.shape = dims
            self.dtype = getdtype(dtype, default=float)
            self.data = array([])
            self.row = array([])
            self.col = array([])
            self._check()
            return
        else:
            raise TypeError, "invalid input format"
        
        self.dtype = getdtype(dtype, obj, default=float)

        try:
            if len(ij) != 2:
                raise TypeError
        except TypeError:
            raise TypeError, "invalid input format"
        
        if dims is None:
            M = int(amax(ij[0])) + 1
            N = int(amax(ij[1])) + 1
            self.shape = (M, N)
        else:
            # Use 2 steps to ensure dims has length 2.
            M, N = dims
            self.shape = (M, N)
        self.row = asarray(ij[0], dtype=numpy.intc)
        self.col = asarray(ij[1], dtype=numpy.intc)
        self.data = asarray(obj, dtype=self.dtype)
        self._check()

    def _check(self):
        """ Checks for consistency and stores the number of non-zeros as
        self.nnz.
        """
        nnz = len(self.data)
        if (nnz != len(self.row)) or (nnz != len(self.col)):
            raise ValueError, "row, column, and data array must all be "\
                  "the same length"
        if (self.row.dtype != numpy.intc):
            self.row = self.row.astype(numpy.intc)
        if (self.col.dtype != numpy.intc):
            self.col = self.col.astype(numpy.intc)

        # some functions pass floats
        self.shape = tuple([int(x) for x in self.shape])
        self.nnz = nnz
        self.ftype = _transtabl.get(self.dtype.char,'')

    def _normalize(self, rowfirst=False):
        if rowfirst:
            #sort by increasing rows first, columns second
            if getattr(self, '_is_normalized', None):
                #columns already sorted, use stable sort for rows
                P = numpy.argsort(self.row, kind='mergesort')        
                return self.data[P], self.row[P], self.col[P]
            else:
                #nothing already sorted
                P  = numpy.lexsort(keys=(self.col, self.row))                
                return self.data[P], self.row[P], self.col[P]
        if getattr(self, '_is_normalized', None):
            return self.data, self.row, self.col
        #sort by increasing rows first, columns second
        P  = numpy.lexsort(keys=(self.row, self.col))        
        self.data, self.row, self.col = self.data[P], self.row[P], self.col[P]
        setattr(self, '_is_normalized', 1)
        return self.data, self.row, self.col

    def rowcol(self, num):
        return (self.row[num], self.col[num])

    def getdata(self, num):
        return self.data[num]

    def tocsc(self):
        if self.nnz == 0:
            return csc_matrix(self.shape, dtype=self.dtype)
        else:
            indptr, rowind, data = cootocsc(self.shape[0], self.shape[1], \
                                            self.size, self.row, self.col, \
                                            self.data)
            return csc_matrix((data, rowind, indptr), self.shape)

    
    def tocsr(self):
        if self.nnz == 0:
            return csr_matrix(self.shape, dtype=self.dtype)
        else:
            indptr, colind, data = cootocsr(self.shape[0], self.shape[1], \
                                            self.size, self.row, self.col, \
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
        self.rows = []
        for i in xrange(M):
            self.rows.append([])
        # The non-zero values of the matrix:
        self.data = [[] for i in xrange(M)]

        if A is not None:
            for i in xrange(A.shape[0]):
                self[i, :] = A[i, :]


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
        
    def __getitem__(self, index):
        """Return the element(s) index=(i, j), where j may be a slice.
        This always returns a copy for consistency, since slices into
        Python lists return copies.
        """
        try:
            assert len(index) == 2
        except (AssertionError, TypeError):
            raise IndexError, "invalid index"
        i, j = index
        if type(i) is slice:
            raise IndexError, "lil_matrix supports slices only of a single row"
            # TODO: add support for this, like in __setitem__
        elif isintlike(i):
            i = int(i)   # Python list indices must be machine-sized ints
            if not (i>=0 and i<self.shape[0]):
                raise IndexError, "lil_matrix index out of range"
        elif operator.isSequenceType(i):
            raise NotImplementedError, "sequence indexing not yet fully supported"
        else:
            raise IndexError, "invalid index"
        row = self.rows[i]
        if type(j) is slice:
            start, stop, stride = j.indices(self.shape[1])
            if stride != 1:
                raise ValueError, "slicing with step != 1 not supported"
            if min(start, stop) < 0 or max(start, stop-1) >= self.shape[1]:
                raise IndexError, "invalid index"
            if stop <= start:
                raise ValueError, "slice width must be >= 1"
            # Look up 'start' and 'stop' in column index
            startind = bisect_left(row, start)
            stopind = bisect_left(row, stop)
            new = lil_matrix((1, stop - start), dtype=self.dtype)
            new.data = [self.data[i][startind:stopind]]
            new.rows = [[colind - start for colind in row[startind:stopind]]]
            return new
        elif operator.isSequenceType(j):
            raise NotImplementedError, "sequence indexing not yet fully supported"
        elif isintlike(j):
            j = int(j)   # Python list indices must be machine-sized ints
            if not (j>=0 and j<self.shape[1]):
                raise IndexError, "lil_matrix index out of range"
        else:
            raise IndexError, "invalid index"
        pos = bisect_left(row, j)
        if pos == len(row) or row[pos] != j:
            # Element doesn't exist (is zero)
            return self.dtype.type(0)
        else:
            return self.data[i][pos]
    
    
    def __setitem__(self, index, x):
        if isscalar(x):
            x = self.dtype.type(x)
        elif not isinstance(x, spmatrix):
            x = numpy.asarray(x, dtype=self.dtype)
        try:
            assert len(index) == 2
        except (AssertionError, TypeError):
            raise IndexError, "invalid index"
        i, j = index
        if isintlike(i):
            i = int(i)   # Python list indices must be machine-sized ints
            if not (i>=0 and i<self.shape[0]):
                raise IndexError, "lil_matrix index out of range"
        else:
            if isinstance(i, slice):
                seq = xrange(i.start or 0, i.stop or self.shape[1], i.step or 1)
            elif operator.isSequenceType(i):
                seq = i
            else:
                raise IndexError, "invalid index"
            try:
                if not len(x) == len(seq):
                    raise ValueError, "number of elements in source must be" \
                            " same as number of elements in destimation"
            except TypeError:
                # Either x or seq is not a sequence.  Note that a sparse matrix
                # is also not a sequence under this definition.
                # Currently we don't support setting to/from non-sequence types.
                # This could be enhanced, though, to allow a scalar source,
                # and/or a sparse vector.
                raise TypeError, "unsupported type for lil_matrix.__setitem__"
            else:
                # Sequence: call __setitem__ recursively, once for each row
                for i in xrange(len(seq)):
                    self[seq[i], index[1]] = x[i]
                return

        # Below here, i is an integer
        row = self.rows[i]
        if isintlike(j):
            j = int(j)   # Python list indices must be machine-sized ints
            if not (j>=0 and j<self.shape[1]):
                raise IndexError, "lil_matrix index out of range"
            # Find element to be set or removed
            pos = bisect_left(row, j)
            if x != 0:
                if pos == len(row):
                    # New element at end
                    row.append(j)
                    self.data[i].append(x)
                elif row[pos] != j:
                    # New element
                    row.insert(pos, j)
                    self.data[i].insert(pos, x)
                else:
                    # Element already exists
                    self.data[i][pos] = x
            else:
                # Remove element
                if pos < len(row) and row[pos] == j:
                    # Element already exists
                    del row[pos]
                    del self.data[i][pos]
                # otherwise no element exists -- do nothing
        else:
            if isinstance(j, slice):
                # Is there an easier way to do this?
                seq = xrange(j.start or 0, j.stop or self.shape[1], j.step or 1)
            elif operator.isSequenceType(j):
                seq = j
            else:
                raise IndexError, "invalid index"

            if isscalar(x):             # does this work with 0-dim arrays too?
                if len(row) == 0:
                    row[:] = seq
                    self.data[i] = [x for item in seq]  # make copies 
                else:
                    # add elements the slow way
                    for k, col in enumerate(seq):
                        self[i, col] = x
                return
            else:
                if isinstance(x, lil_matrix):
                    if x.shape == (1, self.shape[1]):
                        self.rows[i] = x.rows[0]
                        self.data[i] = x.data[0]  
                        # This doesn't make a copy.  Whether to copy or
                        # not with sparse matrix slicing needs a thorough
                        # review for consistency!
                    elif x.shape == (1, len(seq)):
                        # Add elements the slow way.  This could be made
                        # more efficient!
                        for k, col in enumerate(seq):
                            self[i, col] = x[0, k]
                    else:
                        # This should never happen:
                        raise ValueError, "source and destination must have" \
                                          " the same shape"
                    return
                elif isinstance(x, csr_matrix):

                    if x.shape != (1, self.shape[1]):
                        raise ValueError, "sparse matrix source must be (1 x n)"
                    self.rows[i] = x.colind.tolist()
                    self.data[i] = x.data.tolist()
                    # This should be generalized to other shapes than an entire
                    # row.

                else:
                    # Try converting to a dense array.
                    try:
                        x = asarray(x)
                    except Error, e:
                        raise TypeError, "unsupported type for" \
                                         " lil_matrix.__setitem__"
                    else:
                        if x.ndim == 2:
                            # If it's 2d, ensure it's 1 x n, then
                            # squeeze it down to 1d
                            if x.shape != (1, len(seq)):
                                raise ValueError, \
                                    "source and destination have incompatible "\
                                    "dimensions"
                            else:
                                x = x.squeeze()
                     
                    # Add elements one by one
                    for k, col in enumerate(seq):
                        self[i, col] = x[k]
                    return
                    
                return
                

    def __mul__(self, other):           # self * other
        if isscalarlike(other):
            new = self.copy()
            if other == 0:
                # Multiply by zero: return the zero matrix
                return new
            # Multiply this scalar by every element.
            new.data = [[val*other for val in rowvals] for rowvals in new.data]
            return new
        else:
            return self.dot(other)


    def copy(self):
        new = lil_matrix(self.shape, dtype=self.dtype)
        new.data = copy.deepcopy(self.data)
        new.rows = copy.deepcopy(self.rows)
        return new
    
    
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

def _isinstance( x, _class ):
    ##
    # This makes scipy.sparse.sparse.csc_matrix == __main__.csc_matrix.
    c1 = ('%s' % x.__class__).split( '.' )
    c2 = ('%s' % _class).split( '.' )
    aux = c1[-1] == c2[-1]
    return isinstance( x, _class ) or aux

def isspmatrix(x):
    return _isinstance(x, spmatrix)

issparse = isspmatrix

def isspmatrix_csr( x ):
    return _isinstance(x, csr_matrix)

def isspmatrix_csc( x ):
    return _isinstance(x, csc_matrix)

def isspmatrix_dok( x ):
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
        assert M == int(M) and N == int(N)   # raises TypeError unless integers
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

    allowed = 'fdFD'
    if newdtype.char not in allowed:
        if default is None or canCast:
            newdtype = numpy.dtype( 'd' )
        else:
            raise TypeError, "dtype must be one of 'fdFD'"
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
    
def spidentity(n, dtype='d'):
    """
    spidentity( n ) returns the identity matrix of shape (n, n) stored
    in CSC sparse matrix format.
    """
    diags = ones( (1, n), dtype = dtype )
    return spdiags( diags, 0, n, n )


def speye(n, m, k = 0, dtype = 'd'):
    """
    speye(n, m) returns a (n x m) matrix stored
    in CSC sparse matrix format, where the  k-th diagonal is all ones,
    and everything else is zeros.
    """
    diags = ones((1, n), dtype = dtype)
    return spdiags(diags, k, n, m)

def _testme():
    a = csc_matrix((arange(1, 9), \
            transpose([[0, 1, 1, 2, 2, 3, 3, 4], [0, 1, 3, 0, 2, 3, 4, 4]])))
    print "Representation of a matrix:"
    print repr(a)
    print "How a matrix prints:"
    print a
    print "Adding two matrices:"
    b = a+a
    print b
    print "Subtracting two matrices:"
    c = b - a
    print c
    print "Multiplying a sparse matrix by a dense vector:"
    d = a*[1, 2, 3, 4, 5]
    print d
    print [1, 2, 3, 4, 5]*a

    print "Inverting a sparse linear system:"
    print "The sparse matrix (constructed from diagonals):"
    a = spdiags([[1, 2, 3, 4, 5], [6, 5, 8, 9, 10]], [0, 1], 5, 5)
    b = array([1, 2, 3, 4, 5])
    
    print "(Various small tests follow ...)\n"
    print "Dictionary of keys matrix:"
    a = dok_matrix(shape=(10, 10))
    a[1, 1] = 1.
    a[1, 5] = 1.
    print a
    print "Adding it to itself:"
    print a + a

    print "Multiplying by a scalar:"
    print a * 100

    print "Dense representation:"
    print a.todense()

    print "Converting to a CSR matrix:"
    c = a.tocsr()
    print c

if __name__ == "__main__":
    _testme()

