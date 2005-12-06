""" Scipy 2D sparse matrix module.

Original code by Travis Oliphant.

Modified by Ed Schofield and Robert Cimrman.
"""

from scipy.base import *
import types
import sparsetools
import _superlu
import sys
import itertools

def resize1d(arr, newlen):
    old = len(arr)
    new = zeros((newlen,), arr.dtypechar)
    new[:old] = arr
    return new

MAXPRINT=50
ALLOCSIZE=1000

_coerce_rules = {('f', 'f'):'f', ('f', 'd'):'d', ('f', 'F'):'F',
                 ('f', 'D'):'D', ('d', 'f'):'d', ('d', 'd'):'d',
                 ('d', 'F'):'D', ('d', 'D'):'D', ('F', 'f'):'F',
                 ('F', 'd'):'D', ('F', 'F'):'F', ('F', 'D'):'D',
                 ('D', 'f'):'D', ('D', 'd'):'d', ('D', 'F'):'D',
                 ('D', 'D'):'D'}
_transtabl = {'f':'s', 'd':'d', 'F':'c', 'D':'z'}
_itranstabl = {'s':'f', 'd':'d', 'c':'F', 'z':'D'}

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

def _convert_data(data1, data2, newtype):
    if data1.dtypechar != newtype:
        data1 = data1.astype(newtype)
    if data2.dtypechar != newtype:
        data2 = data2.astype(newtype)
    return data1, data2        
    
class spmatrix:
    """ This class provides a base class for all sparse matrices.  It
    cannot be instantiated.  Most of the work is provided by subclasses.
    """
   
    __array_priority__ = 10.1
    def __init__(self, maxprint=MAXPRINT, allocsize=ALLOCSIZE):
        self.format = self.__class__.__name__[:3]
        if self.format == 'spm':
            raise ValueError, "This class is not intended" \
                  " to be instantiated directly."
        self.maxprint = maxprint
        self.allocsize = allocsize
        
    def getmaxprint(self):
        try:
            maxprint = self.maxprint
        except AttributeError:
            maxprint = MAXPRINT
        return maxprint
        
    #def typecode(self):
    #    try:
    #        typ = self.dtypechar
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
        # provides a way to print over a single index
        val = ''
        for ind in xrange(start, stop):
            val = val + '  %s\t%s\n' % (self.rowcol(ind), self.getdata(ind))
        return val
    
    def __repr__(self):
        format = self.getformat()
        return "<%dx%d sparse matrix of type '%s'\n\twith %d stored "\
               "elements (space for %d)\n\tin %s format>" % \
               (self.shape + (self.dtypechar, self.getnnz(), self.nzmax, \
                   _formats[format][1]))

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
        raise TypeError, "comparison of sparse matrices not implemented"

    def __nonzero__(self):  # Simple -- other ideas?
        return self.getnnz() > 0

    # What should len(sparse) return? For consistency with dense matrices,
    # perhaps it should be the number of rows?  For now we return the number of
    # non-zero elements.
    def __len__(self):
        return self.getnnz()

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
            return self.todense()
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
        # csc = self.tocsc()
        # new = csc.transpose()
        # new.data = conj(new.data)
        # return new
 
    def _real(self):
        csc = self.tocsc()
        csc.data = real(csc.data)
        csc.dtypechar = csc.data.dtypechar
        csc.ftype = _transtabl[csc.dtypechar]
        return csc

    def _imag(self):
        csc = self.tocsc()
        csc.data = imag(csc.data)
        csc.dtypechar = csc.data.dtypechar
        csc.ftype = _transtabl[csc.dtypechar]        
        return csc
        
    def dot(self, other):
        """ A generic interface for matrix-matrix or matrix-vector
        multiplication.  Returns A.transpose().conj() * other or
        A.transpose() * other.
        """
        M, K1 = self.shape
        try:
            K2, N = other.shape
        except (AttributeError, TypeError):
            # Not sparse or dense.  Interpret it as a sequence.
            try:
                return self.matvec(asarray(other))
            except:
                raise TypeError, "x.dot(y): y must be matrix, vector, or seq"
        except ValueError:
            # Assume it's a rank-1 array
            K2 = other.shape[0]
            N = 1
        
        if N == 1:
            return self.matvec(other)
        else:
            if K1 != K2:
                raise ValueError, "dimension mismatch"
            return self.matmat(other)

    def matmat(self, other):
        csc = self.tocsc()
        return csc.matmat(other)

    def matvec(self, other):
        csc = self.tocsc()
        return csc.matvec(other)

    def rmatvec(self, other, conjugate=True):
        """ If 'conjugate' is True:
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
        csc = self.tocsc()
        return csc.todense()

    def tocoo(self):
        csc = self.tocsc()
        return csc.tocoo()

    def copy(self):
        csc = self.tocsc()
        return csc.copy()

    def save( self, file_name, format = '%d %d %f\n' ):
        try:
            fd = open( file_name, 'w' )
        except Exception, e:
            raise e, file_name
        
        fd.write( '%d %d\n' % self.shape )
        fd.write( '%d\n' % self.size )
        for ii in xrange( self.size ):
            ir, ic = self.rowcol( ii )
            data = self.getdata( ii )
            fd.write( format % (ir, ic, data) )
        fd.close()
        
class csc_matrix(spmatrix):
    """ Compressed sparse column matrix
        This can be instantiated in several ways:
          - csc_matrix(d)
            with a dense matrix d

          - csc_matrix(s)
            with another sparse matrix s (sugar for .tocsc())

          - csc_matrix((M, N), [nzmax, dtype])
            to construct a container, where (M, N) are dimensions and
            nzmax, dtype are optional, defaulting to nzmax=100 and dtype='d'.

          - csc_matrix((data, ij), [(M, N), nzmax])
            where data, ij satisfy:
                a[ij[k, 0], ij[k, 1]] = data[k]

          - csc_matrix((data, row, ptr), [(M, N)])
            standard CSC representation
    """
    def __init__(self, arg1, dims=(None,None), nzmax=100, dtype='d', copy=False):
        spmatrix.__init__(self)
        if isdense(arg1):
            # Convert the dense matrix arg1 to CSC format
            if rank(arg1) == 2:
                s = asarray(arg1)
                if s.dtypechar not in 'fdFD':
                    # Use a double array as the source (but leave it alone)
                    s = s*1.0            
                if (rank(s) == 2):
                    M, N = s.shape
                    dtype = s.dtypechar
                    func = getattr(sparsetools, _transtabl[dtype]+'fulltocsc')
                    ierr = irow = jcol = 0
                    nnz = sum(ravel(s != 0.0))
                    a = zeros((nnz,), dtype)
                    rowa = zeros((nnz,), 'i')
                    ptra = zeros((N+1,), 'i')
                    while 1:
                        a, rowa, ptra, irow, jcol, ierr = \
                           func(s, a, rowa, ptra, irow, jcol, ierr)
                        if (ierr == 0): break
                        nnz = nnz + ALLOCSIZE
                        a = resize1d(a, nnz)
                        rowa = resize1d(rowa, nnz)
                    self.data = a
                    self.rowind = rowa
                    self.indptr = ptra
                    self.shape = (M, N)
                
                # s = dok_matrix(arg1).tocsc(nzmax)
                # self.shape = s.shape
                # self.data = s.data
                # self.rowind = s.rowind
                # self.indptr = s.indptr
            else:
                raise ValueError, "dense array does not have rank 1 or 2"
        
        elif isspmatrix(arg1):
            s = arg1
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
                func = getattr(sparsetools, s.ftype+'transp')
                self.data, self.rowind, self.indptr = \
                           func(s.shape[1], s.data, s.colind, s.indptr)
            else:
                temp = s.tocsc()
                self.data = temp.data
                self.rowind = temp.rowind
                self.indptr = temp.indptr
                self.shape = temp.shape
        elif type(arg1) == tuple:
            try:
                # Assume it's a tuple of matrix dimensions (M, N)
                (M, N) = arg1
                M = int(M)      # will raise TypeError if (data, ij)
                N = int(N)
                self.data = zeros((nzmax,), dtype)
                self.rowind = zeros((nzmax,), int)
                self.indptr = zeros((N+1,), int)
                self.shape = (M, N)
            except (ValueError, TypeError):
                try:
                    # Try interpreting it as (data, ij)
                    (s, ij) = arg1
                    assert isinstance(ij, ArrayType) and (rank(ij) == 2) and (shape(ij) == (len(s), 2))
                    temp = coo_matrix( s, ij, dims=dims, nzmax=nzmax, \
                            dtype=dtype).tocsc()
                    self.shape = temp.shape
                    self.data = temp.data
                    self.rowind = temp.rowind
                    self.indptr = temp.indptr
                except:
                    try:
                        # Try interpreting it as (data, rowind, indptr)
                        (s, rowind, indptr) = arg1
                        if copy:
                            self.data = array(s)
                            self.rowind = array(rowind)
                            self.indptr = array(indptr)
                        else:
                            self.data = asarray(s)
                            self.rowind = asarray(rowind)
                            self.indptr = asarray(indptr)
                    except:
                        raise ValueError, "unrecognized form for csc_matrix constructor"
        else:
            raise ValueError, "unrecognized form for csc_matrix constructor"

        # Read existing matrix dimensions
        try:
            (oldM, oldN) = self.shape
        except:
            oldM = oldN = None
        # Read matrix dimensions given, if any
        try:
            (M, N) = dims
        except TypeError:
            M = N = None
        M = max(oldM, M, int(amax(self.rowind)) + 1)
        N = max(0, oldN, N, len(self.indptr) - 1)
        self.shape = (M, N)

        self._check()

    def _check(self):
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
        self.nnz = nnz
        self.nzmax = nzmax
        self.dtypechar = self.data.dtypechar
        if self.dtypechar not in 'fdFD':
            self.data = 1.0 * self.data
            self.dtypechar = self.data.dtypechar
        self.ftype = _transtabl[self.dtypechar]
        

    def __radd__(self, other):
        """ Function supporting the operation: self + other.
        This does not currently work correctly for self + dense.
        Perhaps dense matrices need some hooks to support this.
        """
        if isscalar(other) or (isdense(other) and rank(other)==0):
            raise NotImplementedError, 'adding a scalar to a CSC matrix is ' \
                    'not yet supported'
        elif isspmatrix(other):
            ocs = other.tocsc()
            if (ocs.shape != self.shape):
                raise ValueError, "inconsistent shapes"
            dtypechar = _coerce_rules[(self.dtypechar, ocs.dtypechar)]
            nnz1, nnz2 = self.nnz, ocs.nnz
            data1, data2 = _convert_data(self.data[:nnz1], ocs.data[:nnz2], dtypechar)
            func = getattr(sparsetools, _transtabl[dtypechar]+'cscadd')
            c, rowc, ptrc, ierr = func(data1, self.rowind[:nnz1], self.indptr, data2, ocs.rowind[:nnz2], ocs.indptr)
            if ierr:
                raise ValueError, "ran out of space (but shouldn't have happened)"
            M, N = self.shape
            return csc_matrix((c, rowc, ptrc), dims=(M, N))
        elif isdense(other):
            # Convert this matrix to a dense matrix and add them.
            # This does not currently work.
            return self.todense() + other
        else:
            raise TypeError, "unsupported type for sparse matrix addition"

    def __add__(self, other):
        if isscalar(other) or (isdense(other) and rank(other)==0):
            raise NotImplementedError, 'adding a scalar to a CSC matrix is ' \
                    'not yet supported'
        elif isspmatrix(other):
            ocs = other.tocsc()
            if (ocs.shape != self.shape):
                raise ValueError, "inconsistent shapes"
            dtypechar = _coerce_rules[(self.dtypechar, ocs.dtypechar)]
            nnz1, nnz2 = self.nnz, ocs.nnz
            data1, data2 = _convert_data(self.data[:nnz1], ocs.data[:nnz2], dtypechar)
            func = getattr(sparsetools, _transtabl[dtypechar]+'cscadd')
            c, rowc, ptrc, ierr = func(data1, self.rowind[:nnz1], self.indptr, data2, ocs.rowind[:nnz2], ocs.indptr)
            if ierr:
                raise ValueError, "ran out of space (but shouldn't have happened)"
            M, N = self.shape
            return csc_matrix((c, rowc, ptrc), dims=(M, N))
        elif isdense(other):
            # Convert this matrix to a dense matrix and add them
            return other + self.todense()
        else:
            raise TypeError, "unsupported type for sparse matrix addition"

    def __mul__(self, other):
        """ Scalar, vector, or matrix multiplication
        """
        if isscalar(other) or (isdense(other) and rank(other)==0):
            new = self.copy()
            new.data *= other
            new.dtypechar = new.data.dtypechar
            new.ftype = _transtabl[new.dtypechar]
            return new
        else:
            return self.dot(other)
        #else:
        #    return TypeError, "unknown type for sparse matrix multiplication"

    def __rmul__(self, other):  # other * self
        if isscalar(other) or (isdense(other) and rank(other)==0):
            new = self.copy()
            new.data = other * new.data
            new.dtypechar = new.data.dtypechar
            new.ftype = _transtabl[new.dtypechar]
            return new
        else:
            other = asarray(other)
            return self.transpose().dot(other.transpose()).transpose()

    def __neg__(self):
        new = self.copy()
        new.data *= -1
        return new
        
    def __pow__(self, other):  
        """ Element-by-element power (unless other is a scalar, in which
        case return the matrix power.)
        """
        if isscalar(other) or (isdense(other) and rank(other)==0):
            new = self.copy()
            new.data = new.data ** other
            new.dtypechar = new.data.dtypechar
            new.ftype = _transtabl[new.dtypechar]
            return new
        else:
            ocs = other.tocsc()
            if (ocs.shape != self.shape):
                raise ValueError, "inconsistent shapes"
            dtypechar = _coerce_rules[(self.dtypechar, ocs.dtypechar)]
            nnz1, nnz2 = self.nnz, ocs.nnz
            data1, data2 = _convert_data(self.data[:nnz1], ocs.data[:nnz2], dtypechar)
            func = getattr(sparsetools, _transtabl[dtypechar]+'cscmul')
            c, rowc, ptrc, ierr = func(data1, self.rowind[:nnz1], self.indptr, data2, ocs.rowind[:nnz2], ocs.indptr)
            if ierr:
                raise ValueError, "ran out of space (but shouldn't have happened)"
            M, N = self.shape
            return csc_matrix((c, rowc, ptrc), dims=(M, N))

    def transpose(self, copy=False):
        M, N = self.shape
        new = csr_matrix((N, M), nzmax=self.nzmax, dtype=self.dtypechar)
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
        new = csc_matrix(self.shape, nzmax=self.nzmax, dtype=self.dtypechar)
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

    def matvec(self, other):
        if isdense(other):
            if (rank(other) != 1) or (len(other) != self.shape[1]):
                raise ValueError, "dimension mismatch"
            func = getattr(sparsetools, self.ftype+'cscmux')
            y = func(self.data, self.rowind, self.indptr, other, self.shape[0])
            return y
        elif isspmatrix(other):
            raise NotImplementedError, "use matmat() for sparse * sparse"
        else:
            raise TypeError, "need a dense vector"

    def rmatvec(self, other, conjugate=True):
        if isdense(other):
            if (rank(other) != 1) or (len(other) != self.shape[0]):
                raise ValueError, "dimension mismatch"
            func = getattr(sparsetools, self.ftype+'csrmux')
            if conjugate:
                cd = conj(self.data)
            else:
                cd = self.data
            y = func(cd, self.rowind, self.indptr, other)
            return y
        elif isspmatrix(other):
            raise NotImplementedError, "use matmat() for sparse * sparse"
        else:
            raise TypeError, "need a dense vector"

    def matmat(self, other):
        if isspmatrix(other):
            M, K1 = self.shape
            K2, N = other.shape
            if (K1 != K2):
                raise ValueError, "shape mismatch error"
            a, rowa, ptra = self.data, self.rowind, self.indptr
            if isinstance(other, csr_matrix):
                other._check()
                dtypechar = _coerce_rules[(self.dtypechar, other.dtypechar)]
                ftype = _transtabl[dtypechar]            
                func = getattr(sparsetools, ftype+'cscmucsr')
                b = other.data
                rowb = other.colind
                ptrb = other.indptr
            elif isinstance(other, csc_matrix):
                other._check()
                dtypechar = _coerce_rules[(self.dtypechar, other.dtypechar)]
                ftype = _transtabl[dtypechar]                        
                func = getattr(sparsetools, ftype+'cscmucsc')
                b = other.data
                rowb = other.rowind
                ptrb = other.indptr
            else:
                other = other.tocsc()
                dtypechar = _coerce_rules[(self.dtypechar, other.dtypechar)]
                ftype = _transtabl[dtypechar]                        
                func = getattr(sparsetools, ftype+'cscmucsc')
                b = other.data
                rowb = other.rowind
                ptrb = other.indptr
            a, b = _convert_data(a, b, dtypechar)
            newshape = (M, N)
            ptrc = zeros((N+1,), 'i')
            nnzc = 2*max(ptra[-1], ptrb[-1])
            c = zeros((nnzc,), dtypechar)
            rowc = zeros((nnzc,), 'i')
            ierr = irow = kcol = 0
            while 1:
                c, rowc, ptrc, irow, kcol, ierr = func(M, a, rowa, ptra, b, rowb, ptrb, c, rowc, ptrc, irow, kcol, ierr)
                if (ierr==0): break
                # otherwise we were too small and must resize
                #  calculations continue where they left off...
                percent_to_go = 1- (1.0*kcol) / N
                newnnzc = int(ceil((1+percent_to_go)*nnzc))
                c = resize1d(c, newnnzc)
                rowc = resize1d(rowc, newnnzc)
                nnzc = newnnzc
            return csc_matrix((c, rowc, ptrc), dims=(M, N))
        elif isdense(other):
            # This is SLOW!  We need a more efficient implementation
            # of sparse * dense matrix multiplication!
            return self.matmat(csc_matrix(other))
        else:
            raise TypeError, "need a dense or sparse matrix"
            

    def __getitem__(self, key):
        if isinstance(key, types.TupleType):
            row = key[0]
            col = key[1]
            func = getattr(sparsetools, self.ftype+'cscgetel')
            M, N = self.shape
            if not (0<=row<M) or not (0<=col<N):
                raise KeyError, "index out of bounds"
            ind, val = func(self.data, self.rowind, self.indptr, row, col)
            return val
        #elif isinstance(key, type(3)):
        elif type(key) == int:
            return self.data[key]
        else:
            raise NotImplementedError

    def __setitem__(self, key, val):
        if isinstance(key, types.TupleType):
            row = key[0]
            col = key[1]
            func = getattr(sparsetools, self.ftype+'cscsetel')
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
            nzmax = self.nzmax
            if (nzmax < self.nnz+1):  # need more room 
                alloc = max(1, self.allocsize)
                self.data = resize1d(self.data, nzmax + alloc)
                self.rowind = resize1d(self.rowind, nzmax + alloc)
            func(self.data, self.rowind, self.indptr, row, col, val)
            self._check()
        elif isinstance(key, types.IntType):
            if (key < self.nnz):
                self.data[key] = val
            else:
                raise KeyError, "key out of bounds"
        else:
            raise NotImplementedError
            
    def rowcol(self, ind):
        row = self.rowind[ind]
        col = searchsorted(self.indptr, ind+1)-1
        return (row, col)

    def getdata(self, ind):
        return self.data[ind]
    
    def tocsc(self, copy=False):
        if copy:
            new = self.copy()
        else:
            new = self
        return new

    def tocoo(self):
        func = getattr(sparsetools, self.ftype+"csctocoo")
        data, row, col = func(self.data, self.rowind, self.indptr)
        return coo_matrix(data, (row, col), dims=self.shape)

    def tocsr(self):
        return self.tocoo().tocsr()

    def todense(self):
        func = getattr(sparsetools, self.ftype+'csctofull')
        return func(self.shape[0], self.data, self.rowind, self.indptr)

    def prune(self):
        """ Remove empty space after all non-zero elements.
        """
        nnz = self.indptr[-1]
        if self.nzmax <= nnz:
            if self.nzmax < nnz:
                raise RunTimeError, "should never have nnz > nzmax"            
            return
        self.nnz = nnz
        self.data = self.data[:nnz]
        self.rowind = self.rowind[:nnz]
        self.nzmax = nnz
        self._check()

    def copy(self):
        new = csc_matrix(self.shape, nzmax=self.nzmax, dtype=self.dtypechar)
        new.data = self.data.copy()
        new.rowind = self.rowind.copy()
        new.indptr = self.indptr.copy()
        new._check()
        return new
    
 
class csr_matrix(spmatrix):
    """ Compressed sparse row matrix
        This can be instantiated in several ways:
          - csr_matrix(d)
            with a dense matrix d

          - csr_matrix(s)
             with another sparse matrix s (sugar for .tocsr())

          - csr_matrix((M, N), [nzmax, dtype])
             to construct a container, where (M, N) are dimensions and
             nzmax, dtype are optional, defaulting to nzmax=100 and dtype='d'.

          - csr_matrix((data, ij), [(M, N), nzmax])
             where data, ij satisfy:
                a[ij[k, 0], ij[k, 1]] = data[k]

          - csr_matrix((data, col, ptr), [(M, N)])
             standard CSR representation
    """
    def __init__(self, arg1, dims=(None,None), nzmax=100, dtype='d', copy=False):
        spmatrix.__init__(self)
        if isdense(arg1):
            # Convert the dense matrix arg1 to CSR format
            if rank(arg1) == 2:
                s = asarray(arg1)
                ocsc = csc_matrix(transpose(s))
                self.colind = ocsc.rowind
                self.indptr = ocsc.indptr
                self.data = ocsc.data
                self.shape = (ocsc.shape[1], ocsc.shape[0])

        elif isspmatrix(arg1):
            s = arg1
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
            elif isinstance(s, csc_matrix):
                self.shape = s.shape
                func = getattr(sparsetools, s.ftype+'transp')
                self.data, self.colind, self.indptr = \
                           func(s.shape[1], s.data, s.rowind, s.indptr)
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
            try:
                # Assume it's a tuple of matrix dimensions (M, N)
                (M, N) = arg1
                M = int(M)      # will raise TypeError if (data, ij)
                N = int(N)
                self.data = zeros((nzmax,), dtype)
                self.colind = zeros((nzmax,), int)
                self.indptr = zeros((M+1,), int)
                self.shape = (M, N)
            except (ValueError, TypeError):
                try:
                    # Try interpreting it as (data, ij)
                    (s, ij) = arg1
                    assert isinstance(ij, ArrayType) and (rank(ij) == 2) and (shape(ij) == (len(s), 2))
                    ijnew = ij.copy()
                    ijnew[:, 0] = ij[:, 1]
                    ijnew[:, 1] = ij[:, 0]
                    temp = coo_matrix(s, ijnew, dims=dims, nzmax=nzmax,
                                      dtype=dtype).tocsr()
                    self.shape = temp.shape
                    self.data = temp.data
                    self.colind = temp.colind
                    self.indptr = temp.indptr
                except:
                    try:
                        # Try interpreting it as (data, colind, indptr)
                        (s, colind, indptr) = arg1
                        if copy:
                            self.data = array(s)
                            self.colind = array(colind)
                            self.indptr = array(indptr)
                        else:
                            self.data = asarray(s)
                            self.colind = asarray(colind)
                            self.indptr = asarray(indptr)
                    except:
                        raise ValueError, "unrecognized form for csr_matrix constructor"
        else:
            raise ValueError, "unrecognized form for csr_matrix constructor"

        # Read existing matrix dimensions
        try:
            (oldM, oldN) = self.shape
        except:
            oldM = oldN = None
        # Read matrix dimensions given, if any
        try:
            (M, N) = dims
        except TypeError:
            M = N = None
        M = max(0, oldM, M, len(self.indptr) - 1)
        N = max(oldN, N, int(amax(self.colind)) + 1)
        self.shape = (M, N)
        self._check()
        
    def _check(self):
        
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
        self.nnz = nnz
        self.nzmax = nzmax
        self.dtypechar = self.data.dtypechar
        if self.dtypechar not in 'fdFD':
            self.data = self.data + 0.0
            self.dtypechar = self.data.dtypechar
            
        self.ftype = _transtabl[self.dtypechar]

        
    def __add__(self, other):
        # First check if argument is a scalar
        if isscalar(other) or (isdense(other) and rank(other)==0):
            # Now we would add this scalar to every element.
            raise NotImplementedError, 'adding a scalar to a sparse matrix ' \
                    'is not yet supported'
        elif isspmatrix(other):
            ocs = other.tocsr()
            if (ocs.shape != self.shape):
                raise ValueError, "inconsistent shapes"

            dtypechar = _coerce_rules[(self.dtypechar, ocs.dtypechar)]
            data1, data2 = _convert_data(self.data, ocs.data, dtypechar)
            func = getattr(sparsetools, _transtabl[dtypechar]+'cscadd')
            c, colc, ptrc, ierr = func(data1, self.colind, self.indptr, data2, ocs.colind, ocs.indptr)
            if ierr:
                raise ValueError, "ran out of space (but shouldn't have happened)"
            M, N = self.shape
            return csr_matrix((c, colc, ptrc), dims=(M, N))
        elif isdense(other):
            # Convert this matrix to a dense matrix and add them.
            # This does not currently work.
            return self.todense() + other
        else:
            raise TypeError, "unsupported type for sparse matrix addition"
    
    def __mul__(self, other):
        """ Scalar, vector, or matrix multiplication
        """
        if isscalar(other) or (isdense(other) and rank(other)==0):
            new = self.copy()
            new.data = other * new.data         # allows type conversion
            new.dtypechar = new.data.dtypechar
            new.ftype = _transtabl[new.dtypechar]
            return new
        else:
            return self.dot(other)
    
    def __rmul__(self, other):  # other * self
        if isscalar(other) or (isdense(other) and rank(other)==0):
            new = self.copy()
            new.data = other * new.data         # allows type conversion
            new.dtypechar = new.data.dtypechar
            new.ftype = _transtabl[new.dtypechar]
            return new
        else:
            other = asarray(other)
            return self.transpose().dot(other.transpose()).transpose()
    
    def __neg__(self):
        new = self.copy()
        new.data *= -1
        return new
        
    def __pow__(self, other):  
        """ Element-by-element power (unless other is a scalar, in which
        case return the matrix power.)
        """
        if isscalar(other) or (isdense(other) and rank(other)==0):
            new = self.copy()
            new.data = new.data ** other
            new.dtypechar = new.data.dtypechar
            new.ftype = _transtabl[new.dtypechar]
            return new
        elif isspmatrix(other):
            ocs = other.tocsr()
            if (ocs.shape != self.shape):
                raise ValueError, "inconsistent shapes"
            dtypechar = _coerce_rules[(self.dtypechar, ocs.dtypechar)]
            data1, data2 = _convert_data(self.data, ocs.data, dtypechar)
            func = getattr(sparsetools, _transtabl[dtypechar]+'cscmul')
            c, colc, ptrc, ierr = func(data1, self.colind, self.indptr, data2, ocs.colind, ocs.indptr)
            if ierr:
                raise ValueError, "ran out of space (but shouldn't have happened)"
            M, N = self.shape
            return csr_matrix((c, colc, ptrc), dims=(M, N))
        else:
            raise TypeError, "unsupported type for sparse matrix power"

    def transpose(self, copy=False):
        M, N = self.shape
        new = csc_matrix((N, M), nzmax=self.nzmax, dtype=self.dtypechar)
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

    def matvec(self, other):
        if (rank(other) != 1) or (len(other) != self.shape[1]):
            raise ValueError, "dimension mismatch"
        func = getattr(sparsetools, self.ftype+'csrmux')
        y = func(self.data, self.colind, self.indptr, other)
        return y

    def rmatvec(self, other, conjugate=True):
        if (rank(other) != 1) or (len(other) != self.shape[0]):
            raise ValueError, "dimension mismatch"
        func = getattr(sparsetools, self.ftype+'cscmux')
        if conjugate:
            cd = conj(self.data)
        else:
            cd = self.data
        y = func(cd, self.colind, self.indptr, other, self.shape[1])
        return y

    def matmat(self, other):
        if isspmatrix(other):
            M, K1 = self.shape
            K2, N = other.shape
            a, rowa, ptra = self.data, self.colind, self.indptr
            if (K1 != K2):
                raise ValueError, "shape mismatch error"
            if isinstance(other, csc_matrix):            
                other._check()
                dtypechar = _coerce_rules[(self.dtypechar, other.dtypechar)]
                ftype = _transtabl[dtypechar]            
                func = getattr(sparsetools, ftype+'csrmucsc')
                b = other.data
                colb = other.rowind
                ptrb = other.indptr
                out = 'csc'
                firstarg = ()
            elif isinstance(other, csr_matrix):
                other._check()
                dtypechar = _coerce_rules[(self.dtypechar, other.dtypechar)]
                ftype = _transtabl[dtypechar]            
                func = getattr(sparsetools, ftype+'cscmucsc')
                b, colb, ptrb = a, rowa, ptra
                a, rowa, ptra = other.data, other.colind, other.indptr
                out = 'csr'
                firstarg = (N,)
            else:
                other = other.tocsc()
                dtypechar = _coerce_rules[(self.dtypechar, other.dtypechar)]
                ftype = _transtabl[dtypechar]            
                func = getattr(sparsetools, ftype+'csrmucsc')
                b = other.data
                colb = other.rowind
                ptrb = other.indptr
                out = 'csc'
                firstarg = ()
            a, b = _convert_data(a, b, dtypechar)            
            newshape = (M, N)
            if out == 'csr':
                ptrc = zeros((M+1,), 'i')
            else:
                ptrc = zeros((N+1,), 'i')
            nnzc = 2*max(ptra[-1], ptrb[-1])
            c = zeros((nnzc,), dtypechar)
            rowc = zeros((nnzc,), 'i')
            ierr = irow = kcol = 0
            while 1:
                args = firstarg+(a, rowa, ptra, b, colb, ptrb, c, rowc, ptrc, irow,
                                 kcol, ierr)
                c, rowc, ptrc, irow, kcol, ierr = func(*args)
                if (ierr==0): break
                # otherwise we were too small and must resize
                percent_to_go = 1- (1.0*kcol) / N
                newnnzc = int(ceil((1+percent_to_go)*nnzc))
                c = resize1d(c, newnnzc)
                rowc = resize1d(rowc, newnnzc)
                nnzc = newnnzc

            if out == 'csr':
                # Note: 'rowc' is deliberate
                return csr_matrix((c, rowc, ptrc), dims=(M, N))
            else:
                return csc_matrix((c, rowc, ptrc), dims=(M, N))
        elif isdense(other):
            # This is SLOW!  We need a more efficient implementation
            # of sparse * dense matrix multiplication!
            return self.matmat(csc_matrix(other))
        else:
            raise TypeError, "need a dense or sparse matrix"

    def __getitem__(self, key):
        if isinstance(key, types.TupleType):
            row = key[0]
            col = key[1]
            func = getattr(sparsetools, self.ftype+'cscgetel')
            M, N = self.shape
            if (row < 0):
                row = M + row
            if (col < 0):
                col = N + col
            if (row >= M ) or (col >= N) or (row < 0) or (col < 0):
                raise IndexError, "index out of bounds"
            ind, val = func(self.data, self.colind, self.indptr, col, row)
            return val
        #elif isinstance(key, type(3)):
        elif type(key) == int:
            return self.data[key]
        else:
            raise NotImplementedError


    def __setitem__(self, key, val):
        if isinstance(key, types.TupleType):
            row = key[0]
            col = key[1]
            func = getattr(sparsetools, self.ftype+'cscsetel')
            M, N = self.shape
            if (row < 0):
                row = M + row
            if (col < 0):
                col = N + col
            if (row < 0) or (col < 0):
                raise KeyError, "index out of bounds"
            if (row >= M):
                self.indptr = resize1d(self.indptr, row+2)
                self.indptr[M+1:] = self.indptr[M]
                M = row+1
            elif (row < 0):
                row = M - row
            if (col >= N):
                N = col+1
            elif (col < 0):
                col = N - col
            self.shape = (M, N)
            nzmax = self.nzmax
            if (nzmax < self.nnz+1):  # need more room 
                alloc = max(1, self.allocsize)
                self.data = resize1d(self.data, nzmax + alloc)
                self.colind = resize1d(self.colind, nzmax + alloc)
            func(self.data, self.colind, self.indptr, col, row, val)
            self._check()
        elif isinstance(key, types.IntType):
            if (key < self.nnz):
                self.data[key] = val
            else:
                raise KeyError, "key out of bounds"
        else:
            raise NotImplementedError
            
    def rowcol(self, ind):
        col = self.colind[ind]
        row = searchsorted(self.indptr, ind+1)-1
        return (row, col)

    def getdata(self, ind):
        return self.data[ind]

    def tocsr(self, copy=False):
        if copy:
            return self.copy()
        else:
            return self

    def tocoo(self):
        func = getattr(sparsetools, self.ftype+"csctocoo")
        data, col, row = func(self.data, self.colind, self.indptr)
        return coo_matrix(data, (row, col), dims=self.shape)

    def tocsc(self):
        return self.tocoo().tocsc()

    def todense(self):
        func = getattr(sparsetools, self.ftype+'csctofull')
        s = func(self.shape[1], self.data, self.colind, self.indptr)
        return transpose(s)

    def prune(self):
        """ Eliminate non-zero entries, leaving space for at least
        newnzmax elements.
        """
        nnz = self.indptr[-1]
        if self.nzmax <= nnz:
            if self.nzmax < nnz:
                raise RunTimeError, "should never have nnz > nzmax"
            return
        self.data = self.data[:nnz]
        self.colind = self.colind[:nnz]
        self.nzmax = nnz
        self._check()

    def copy(self):
        new = csr_matrix(self.shape, nzmax=self.nzmax, dtype=self.dtypechar)
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
    """ A dictionary of keys based matrix.  This is relatively efficient
    for constructing sparse matrices for conversion to other sparse
    matrix types.
    
    It does type checking on input by default.  To disable this type
    checking and speed up element accesses slightly, set self._validate
    to False.
    """
    def __init__(self, A=None):
        """ Create a new dictionary-of-keys sparse matrix.  An optional
        argument A is accepted, which initializes the dok_matrix with it.
        This can be a tuple of dimensions (m, n) or a (dense) array
        to copy.
        """
        dict.__init__(self)
        spmatrix.__init__(self)
        self.shape = (0, 0)
        # If _validate is True, ensure __setitem__ keys are integer tuples
        self._validate = True
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
            if isspmatrix(A):
                # For sparse matrices, this is too inefficient; we need 
                # something else.
                raise NotImplementedError, "initializing a dok_matrix with " \
                        "a sparse matrix is not yet supported"
            elif isdense(A):
                A = asarray(A)
                if rank(A) == 2:
                    M, N = A.shape
                    self.shape = (M, N)
                    for i in range(M):
                        for j in range(N):
                            if A[i, j] != 0:
                                self[i, j] = A[i, j]
                elif rank(A) == 1:
                    M = A.shape[0]
                    self.shape = (M, 1)
                    for i in range(M):
                        if A[i] != 0:
                            self[i, 0] = A[i]
                else:
                    raise TypeError, "array for initialization must have rank 2"
            else:
                raise TypeError, "argument should be a tuple of dimensions " \
                        "or a sparse or dense matrix"

    def getnnz(self):
        return dict.__len__(self)
    
    def __len__(self):
        return dict.__len__(self)
    
    def __str__(self):
        val = ''
        nnz = len(self)
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

    def __repr__(self):
        nnz = len(self)
        format = self.getformat()
        return "<%dx%d sparse matrix with %d stored "\
               "elements in %s format>" % \
               (self.shape + (nnz, _formats[format][1]))

    def __getitem__(self, key):
        if self._validate:
            # Sanity checks: key must be a pair of integers
            if not isinstance(key, tuple) or len(key) != 2:
                raise TypeError, "key must be a tuple of two integers"
            if type(key[0]) != int or type(key[1]) != int:
                raise TypeError, "key must be a tuple of two integers"

        return self.get(key, 0)

    def __setitem__(self, key, value):
        if self._validate:
            # Sanity checks: key must be a pair of integers
            try:
                # Cast to integers and compare.  This way the key
                # can be a pair of rank-0 arrays.
                i, j = int(key[0]), int(key[1])
                assert i == key[0] and j == key[1]
            except:
                raise TypeError, "key must be a tuple of two integers"
             
            if (value == 0):
                if self.has_key(key):  # get rid of it something already there
                    del self[key]              
                return                 # do nothing
            dict.__setitem__(self, key, value)
            newrows = max(self.shape[0], int(key[0])+1)
            newcols = max(self.shape[1], int(key[1])+1)
            self.shape = (newrows, newcols)
        else:
            # Faster version without sanity checks
            if (value == 0):
                if self.has_key(key):  # get rid of it something already there
                    del self[key]              
                return                 # do nothing
            dict.__setitem__(self, key, value)
    
    def __add__(self, other):
        # First check if argument is a scalar
        if isscalar(other) or (isdense(other) and rank(other)==0):
            new = dok_matrix()
            # Add this scalar to every element.
            M, N = self.shape
            for i in range(M):
                for j in range(N):
                    aij = self.get((i, j), 0) + other
                    if aij != 0:
                        new[i, j] = aij
            #new.dtypechar = self.dtypechar
        elif isinstance(other, dok_matrix):
            new = dok_matrix()
            new.update(self)
            new.shape = self.shape
            for key in other.keys():
                new[key] += other[key]
        elif isspmatrix(other):
            csc = self.tocsc()
            new = csc + other
        else:
            # Perhaps it's a dense matrix?
            new = self.todense() + other
        return new

    def __radd__(self, other):
        # First check if argument is a scalar
        if isscalar(other) or (isdense(other) and rank(other)==0):
            new = dok_matrix()
            # Add this scalar to every element.
            M, N = self.shape
            for i in range(M):
                for j in range(N):
                    aij = self.get((i, j), 0) + other
                    if aij != 0:
                        new[i, j] = aij
            #new.dtypechar = self.dtypechar
        elif isinstance(other, dok_matrix):
            new = dok_matrix()
            new.update(self)
            new.shape = self.shape
            for key in other.keys():
                new[key] += other[key]
        elif isspmatrix(other):
            csc = self.tocsc()
            new = csc + other
        else:
            # Perhaps it's a dense matrix?
            new = other + self.todense()
        return new

    def __neg__(self):
        new = dok_matrix()
        for key in self.keys():
            new[key] = -self[key]
        return new

    def __mul__(self, other):           # self * other
        if isscalar(other) or (isdense(other) and rank(other)==0):
            new = dok_matrix()
            # Multiply this scalar by every element.
            for (key, val) in self.items():
                new[key] = val * other
            #new.dtypechar = self.dtypechar
            return new
        else:
            return self.dot(other)

    def __rmul__(self, other):          # other * self
        if isscalar(other) or (isdense(other) and rank(other)==0):
            new = dok_matrix()
            # Multiply this scalar by every element.
            for (key, val) in self.items():
                new[key] = other * val
            #new.dtypechar = self.dtypechar
            return new
        else:
            other = asarray(other)
            return self.transpose().dot(other.transpose()).transpose()

    # What should len(sparse) return? For consistency with dense matrices,
    # perhaps it should be the number of rows?  For now it returns the number
    # of non-zeros.

    def transpose(self):
        """ Return the transpose
        """
        newshape = (self.shape[1], self.shape[0])
        new = dok_matrix(newshape)
        for key in self.keys():
            new[key[1], key[0]] = self[key]
        return new

    def conjtransp(self):
        """ Return the conjugate transpose
        """
        new = dok_matrix()
        for key in self.keys():
            new[key[1], key[0]] = conj(self[key])
        return new

    def copy(self):
        new = dok_matrix()
        new.update(self)
        new.shape = self.shape
        return new
        
    def take(self, cols_or_rows, columns=1):
        # Extract columns or rows as indictated from matrix
        # assume cols_or_rows is sorted
        new = dok_matrix()
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
        # similar to take but returns two array, the extracted
        #  columns plus the resulting array
        #  assumes cols_or_rows is sorted
        base = dok_matrix()
        ext = dok_matrix()
        indx = int((columns == 1))
        N = len(cols_or_rows)
        if indx:
            for key in self.keys():
                num = searchsorted(cols_or_rows, key[1])
                if cols_or_rows[num]==key[1]:
                    newkey = (key[0], num)
                    ext[newkey] = self[key]
                else:
                    newkey = (key[0], key[1]-num)
                    base[newkey] = self[key]
        else:
            for key in self.keys():
                num = searchsorted(cols_or_rows, key[0])
                if cols_or_rows[num]==key[0]:
                    newkey = (num, key[1])
                    ext[newkey] = self[key]
                else:
                    newkey = (key[0]-num, key[1])
                    base[newkey] = self[key]            
        return base, ext


    def matvec(self, other):
        other = asarray(other)
        if other.shape[0] != self.shape[1]:
            raise ValueError, "dimensions do not match"
        new = [0]*self.shape[0]
        for key in self.keys():
            new[int(key[0])] += self[key] * other[int(key[1]), ...]
        return array(new)        

    def rmatvec(self, other, conjugate=True):
        other = asarray(other)
	if other.shape[-1] != self.shape[0]:
	    raise ValueError, "dimensions do not match"
	new = [0]*self.shape[1]
	for key in self.keys():
            new[int(key[1])] += other[..., int(key[0])] * conj(self[key])
	return array(new)

    def setdiag(self, values, k=0):
        M, N = self.shape
        m = len(values)
        for i in range(min(M, N-k)):
            self[i, i+k] = values[i]
        return

    def tocsr(self, nzmax=None):
        """ Return Compressed Sparse Row format arrays for this matrix
        """
        keys = self.keys()
        keys.sort()
        nnz = len(keys)
        nzmax = max(nnz, nzmax)
        data = [0]*nzmax
        colind = [0]*nzmax
        # Empty rows will leave row_ptr dangling.  We assign row_ptr[i] 
        # for each empty row i to point off the end.  Is this sufficient??
        row_ptr = [nnz]*(self.shape[0]+1)
        current_row = -1
        k = 0
        for key in keys:
            ikey0 = int(key[0])
            ikey1 = int(key[1])
            if ikey0 != current_row:
                N = ikey0-current_row
                row_ptr[current_row+1:ikey0+1] = [k]*N
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
        #  Sort based on columns
        # This works, but is very slow for matrices with many non-zero
        # elements (requiring a function call for every element)
        #keys.sort(csc_cmp)

        # Faster sort: Schwartzian transform
        keys = [(k[1], k[0]) for k in self.keys()]
        keys.sort()
        keys = [(k[1], k[0]) for k in keys]
        
        nnz = len(keys)
        nzmax = max(nnz, nzmax)
        data = [0]*nzmax
        rowind = [0]*nzmax
        # Empty columns will leave col_ptr dangling.  We assign col_ptr[j] 
        # for each empty column j to point off the end.  Is this sufficient??
        col_ptr = [nnz]*(self.shape[1]+1)
        current_col = -1
        k = 0
        for key in keys:
            ikey0 = int(key[0])
            ikey1 = int(key[1])
            if ikey1 != current_col:
                N = ikey1-current_col
                col_ptr[current_col+1:ikey1+1] = [k]*N
                current_col = ikey1
            data[k] = self[key]
            rowind[k] = ikey0
            k += 1
        data = array(data)
        rowind = array(rowind)
        col_ptr = array(col_ptr)
        return csc_matrix((data, rowind, col_ptr), dims=self.shape, nzmax=nzmax)

    def todense(self, dtype=None):
        if dtype is None:
            dtype = 'd'
        new = zeros(self.shape, dtype=dtype)
        for key in self.keys():
            ikey0 = int(key[0])
            ikey1 = int(key[1])
            new[ikey0, ikey1] = self[key]
        if amax(ravel(abs(new.imag))) == 0:
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
            for (i,j) in self.keys():
                if i >= newM or j >= newN:
                    del self[i,j]
        self.shape = (newM, newN)


# dictionary of dictionaries based matrix
class dod_matrix(spmatrix):
    pass

# linked list matrix
class lnk_matrix(spmatrix):
    pass


class coo_matrix(spmatrix):
    """ A sparse matrix in coordinate list format.

    COO matrices are instantiated as follows:
        A = coo_matrix(obj, ij, [dims])
    The dimensions are optional.  If supplied, we set (M, N) = dims.
    If not supplied, we infer these from the index arrays
    ij[:][0] and ij[:][1]

    The arguments 'obj' and 'ij' represent three arrays:
        1. obj[:]: the entries of the matrix, in any order
        2. ij[:][0]: the row indices of the matrix entries
        3. ij[:][1]: the column indices of the matrix entries

    So the following holds:
        A[ij[k][0], ij[k][1]] = obj[k]
    """
    def __init__(self, obj, ij_in, dims=None, nzmax=None, dtype=None):
        spmatrix.__init__(self)
        try:
            # Assume the first calling convention
            #            assert len(ij) == 2
            if len(ij_in) != 2:
                if isdense( ij_in ) and (ij_in.shape[1] == 2):
                    ij = (ij_in[:,0], ij_in[:,1])
                else:
                    raise AssertionError
            else:
                ij = ij_in
            if dims is None:
                M = int(amax(ij[0]))
                N = int(amax(ij[1]))
                self.shape = (M, N)
            else:
                # Use 2 steps to ensure dims has length 2.
                M, N = dims
                self.shape = (M, N)
            self.row = asarray(ij[0])
            self.col = asarray(ij[1])
            self.data = asarray(obj, dtype=dtype)
            self.dtypechar = self.data.dtypechar
            if nzmax is None:
                nzmax = len(self.data)
            self.nzmax = nzmax
            self._check()
        except Exception, e:
            raise e, "invalid input format"

    def _check(self):
        """ Checks for consistency and stores the number of non-zeros as
        self.nnz.
        """
        nnz = len(self.data)
        if (nnz != len(self.row)) or (nnz != len(self.col)):
            raise ValueError, "row, column, and data array must all be "\
                  "the same length"
        if (self.nzmax < nnz):
            raise ValueError, "nzmax must be >= nnz"
        self.nnz = nnz
        self.ftype = _transtabl[self.dtypechar]

    def _normalize(self, rowfirst=False):
        if rowfirst:
            l = zip(self.row, self.col, self.data)
            l.sort()
            row, col, data = list(itertools.izip(*l))
            return data, row, col
        if getattr(self, '_is_normalized', None):
            return self.data, self.row, self.col
        l = zip(self.col, self.row, self.data)
        l.sort()
        col, row, data = list(itertools.izip(*l))
        self.col = asarray(col, 'i')
        self.row = asarray(row, 'i')
        self.data = array(data, self.dtypechar)
        setattr(self, '_is_normalized', 1)
        return self.data, self.row, self.col

    def rowcol(self, num):
        return (self.row[num], self.col[num])

    def getdata(self, num):
        return self.data[num]

    def tocsc(self):
        func = getattr(sparsetools, self.ftype+"cootocsc")
        data, row, col = self._normalize()
        a, rowa, ptra, ierr = func(self.shape[1], data, row, col)
        if ierr:
            raise RuntimeError, "error in conversion"
        return csc_matrix((a, rowa, ptra), dims=self.shape)

    def tocsr(self):
        func = getattr(sparsetools, self.ftype+"cootocsc")
        data, row, col = self._normalize(rowfirst=True)
        a, cola, ptra, ierr = func(self.shape[0], data, col, row)
        if ierr:
            raise RuntimeError, "error in conversion"
        return csr_matrix((a, cola, ptra), dims=self.shape)

    def tocoo(self, copy=False):
        if copy:
            return self.copy()
        else:
            return self
        
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
               
def isspmatrix(x):
    return isinstance(x, spmatrix)

def isspmatrix_csr( x ):
    return isinstance(x, csr_matrix)

def isspmatrix_csc( x ):
    return isinstance(x, csc_matrix)

def isspmatrix_dok( x ):
    return isinstance(x, dok_matrix)

def isspmatrix_dod( x ):
    return isinstance(x, dod_matrix)

def isspmatrix_lnk( x ):
    return isinstance(x, lnk_matrix)

def isspmatrix_coo( x ):
    return isinstance(x, coo_matrix)

def isdense(x):
    # What's the best way to check for this?  The following fails on import:
    # import numerictypes
    # return numerictypes.isdtype(x)
    temp = zeros(1)
    return type(x) == type(temp)

def _spdiags_tosub(diag_num, a, b):
    part1 = where(less(diag_num, a), abs(diag_num-a), 0)
    part2 = where(greater(diag_num, b), abs(diag_num-b), 0)
    return part1+part2
                      
def spdiags(diags, offsets, M, N):
    """Return a sparse matrix in CSR format given its diagonals.
    
    B = spdiags(diags, offsets, M, N)

    Inputs:
        diags  --  rows contain diagonal values
        offsets -- diagonals to set (0 is main)
        M, N    -- sparse matrix returned is M X N
    """
    diags = array(transpose(diags), copy=True)
    if diags.dtypechar not in 'fdFD':
        diags = diags.astype('d')
    offsets = array(offsets, copy=False)
    mtype = diags.dtypechar
    assert(len(offsets) == diags.shape[1])
    # set correct diagonal to csr conversion routine for this type
    diagfunc = eval('sparsetools.'+_transtabl[mtype]+'diatocsr')
    a, rowa, ptra, ierr = diagfunc(M, N, diags, offsets)
    if ierr:
        raise ValueError, "ran out of memory (shouldn't have happened)"
    return csc_matrix((a, rowa, ptra), dims=(M, N))

def solve(A, b, permc_spec=2):
    if not hasattr(A, 'tocsr') and not hasattr(A, 'tocsc'):
        raise ValueError, "sparse matrix must be able to return CSC format--"\
              "A.tocsc()--or CSR format--A.tocsr()"
    if not hasattr(A, 'shape'):
        raise ValueError, "sparse matrix must be able to return shape (rows, cols) = A.shape"
    M, N = A.shape
    if (M != N):
        raise ValueError, "matrix must be square"    
    if hasattr(A, 'tocsc') and not isspmatrix_csr( A ):
        mat = A.tocsc()
        ftype, lastel, data, index0, index1 = \
               mat.ftype, mat.nnz, mat.data, mat.rowind, mat.indptr
        csc = 1
    else:
        mat = A.tocsr()
        ftype, lastel, data, index0, index1 = \
               mat.ftype, mat.nnz, mat.data, mat.colind, mat.indptr
        csc = 0
    gssv = eval('_superlu.' + ftype + 'gssv')
    return gssv(N, lastel, data, index0, index1, b, csc, permc_spec)[0]
    

def lu_factor(A, permc_spec=2, diag_pivot_thresh=1.0,
              drop_tol=0.0, relax=1, panel_size=10):
    M, N = A.shape
    if (M != N):
        raise ValueError, "can only factor square matrices"
    csc = A.tocsc()
    gstrf = eval('_superlu.' + csc.ftype + 'gstrf')
    return gstrf(N, csc.nnz, csc.data, csc.rowind, csc.indptr, permc_spec,
                 diag_pivot_thresh, drop_tol, relax, panel_size)
        

if __name__ == "__main__":
    a = csc_matrix((arange(1, 9), transpose([[0, 1, 1, 2, 2, 3, 3, 4], [0, 1, 3, 0, 2, 3, 4, 4]])))
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
    print a
    print "Solve: single precision complex:"
    a = a.astype('F')
    x = solve(a, b)
    print x
    print "Error: ", a*x-b

    print "Solve: double precision complex:"
    a = a.astype('D')
    x = solve(a, b)
    print x
    print "Error: ", a*x-b

    print "Solve: double precision:"
    a = a.astype('d')
    x = solve(a, b)
    print x
    print "Error: ", a*x-b

    print "Solve: single precision:"
    a = a.astype('f')
    x = solve(a, b.astype('f'))
    print x
    print "Error: ", a*x-b

    print "(Various small tests follow ...)\n"
    print "Dictionary of keys matrix:"
    a = dok_matrix()
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

    print "Adding a constant:"
    c += 5
    print c
