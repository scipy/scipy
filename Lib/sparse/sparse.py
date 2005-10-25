""" Scipy 2D sparse matrix module.

Original code by Travis Oliphant.

Modified by Ed Schofield and Robert Cimrman.
"""

from scipy.base import *
import types
import sparsetools
import _superlu
import sys

def resize1d(arr, newlen):
    old = len(arr)
    new = zeros((newlen,),arr.dtypechar)
    new[:old] = arr
    return new

MAXPRINT=50
ALLOCSIZE=1000

_coerce_rules = {('f','f'):'f', ('f','d'):'d', ('f','F'):'F',
                 ('f','D'):'D', ('d','f'):'d', ('d','d'):'d',
                 ('d','F'):'D', ('d','D'):'D', ('F','f'):'F',
                 ('F','d'):'D', ('F','F'):'F', ('F','D'):'D',
                 ('D','f'):'D', ('D','d'):'d', ('D','F'):'D',
                 ('D','D'):'D'}
_transtabl = {'f':'s','d':'d','F':'c','D':'z'}
_itranstabl = {'s':'f','d':'d','c':'F','z':'D'}

# The formats that we might potentially understand.
_formats = {'csc':[0,"Compressed Sparse Column"],
            'csr':[1,"Compressed Sparse Row"],
            'dok':[2,"Dictionary Of Keys"],
            'lil':[3,"LInked List"],
            'dod':[4,"Dictionary of Dictionaries"],
            'sss':[5,"Symmetric Sparse Skyline"],
            'coo':[6,"COOrdinate"],
            'lba':[7,"Linpack BAnded"],
            'egd':[8,"Ellpack-itpack Generalized Diagonal"],
            'dia':[9,"DIAgonal"],
            'bsr':[10,"Block Sparse Row"],
            'msr':[11,"Modified compressed Sparse Row"],
            'bsc':[12,"Block Sparse Column"],
            'msc':[13,"Modified compressed Sparse Column"],
            'ssk':[14,"Symmetric SKyline"],
            'nsk':[15,"Nonsymmetric SKyline"],
            'jad':[16,"JAgged Diagonal"],
            'uss':[17,"Unsymmetric Sparse Skyline"],
            'vbr':[18,"Variable Block Row"],
            'und':[19,"Undefined"]
            }

def _convert_data(data1,data2,newtype):
    if data1.dtypechar != newtype:
        data1 = data1.astype(newtype)
    if data2.dtypechar != newtype:
        data2 = data2.astype(newtype)
    return data1, data2        
    
class spmatrix:
    """ This class provides a base class for all sparse matrices.  It
    cannot be instantiated.  Most of the work is provided by subclasses.
    """
    
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
    #        typ = self._dtypechar
    #    except AttributeError:
    #        typ = None
    #    return typ

    def getnnz(self):
        try:
            nnz = self.nnz
        except AttributeError:
            nnz = 0
        return nnz

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

    def listprint(self,start,stop):
        # provides a way to print over a single index
        val = ''
        for ind in xrange(start,stop):
            val = val + '  %s\t%s\n' % (self.rowcol(ind),self.getdata(ind))
        return val
    
    def __repr__(self):
        typ = self._dtypechar
        nnz = self.getnnz()
        format = self.getformat()
        nzmax = self.getnzmax()
        return "<%dx%d sparse matrix of type '%s' with %d stored "\
               "elements (space for %d) in\n\t%s format>" % \
               (self.shape + (typ, nnz, nzmax, _formats[format][1]))

    def __str__(self):
        nnz = self.getnnz()
        maxprint = self.getmaxprint()
        val = ''
        if nnz > maxprint:
            val = val + self.listprint(0,maxprint/2)
            val = val + "  :\t:\n"
            val = val + self.listprint(nnz-maxprint/2,nnz)
        else:
            val = val + self.listprint(0,nnz)
        return val[:-1]

    def __cmp__(self,other):
        raise TypeError, "Comparison of sparse matrices is not implemented."

    def __nonzero__(self):  # Simple -- other ideas?
        return self.getnnz() > 0

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

    def __sub__(self, other):
        csc = self.tocsc()
        return csc - other

    def __rsub__(self, other): # other - self
        csc = self.tocsc()
        return csc.__rsub__(other)


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

    def __getattr__(self,attr):
        if attr == 'A':
            return self.todense()
        elif attr == 'T':
            return self.transpose()
        elif attr == 'H':
            return self.conjtranspose()
        elif attr == 'real':
            return self._real()
        elif attr == 'imag':
            return self._imag()
        else:
            raise AttributeError, attr + " not found."
        
    def transpose(self):
        csc = self.tocsc()
        return csc.transpose()

    def conjtranspose(self):
        csc = self.tocsc()
        res = csc.transpose()
        res.data = conj(res.data)
        return res

    def _real(self):
        csc = self.tocsc()
        csc.data = real(csc.data)
        csc._dtypechar = csc.data.dtypechar
        csc.ftype = _transtabl[csc._dtypechar]
        return csc

    def _imag(self):
        csc = self.tocsc()
        csc.data = imag(csc.data)
        csc._dtypechar = csc.data.dtypechar
        csc.ftype = _transtabl[csc._dtypechar]        
        return csc
        
    def matrixmultiply(self, other):
        """ A generic interface for matrix-matrix or matrix-vector
        multiplication.
        """
        csc = self.tocsc()
        return csc.matrixmultiply(other)
        
    def matmat(self, other):
        csc = self.tocsc()
        return csc.matmat(other)

    def matvec(self, vec):
        csc = self.tocsc()
        return csc.matvec(vec)

    # implements A.H * x or A.T * x depending on conj
    def rmatvec(self, vec, conj=1):
        csc = self.tocsc()
        return csc.rmatvec(vec, conj=conj)

    def todense(self):
        csc = self.tocsc()
        return csc.todense()

    def tocoo(self):
        csc = self.tocsc()
        return csc.tocoo()

    def copy(self):
        csc = self.tocsc()
        return csc.copy()


 
class csc_matrix(spmatrix):
    """ Compressed sparse column matrix
        This can be instantiated in two ways:
          - csc_matrix(s)
            with another sparse matrix s (sugar for .tocsc())

          - csc_matrix((M,N), [nzmax, dtype])
            to construct a container, where (M,N) are dimensions and
            nzmax, dtype are optional, defaulting to nzmax=100 and dtype='d'.

          - csc_matrix(data, ij, [(M,N), nzmax])
            where data, ij satisfy:
                a[ij[k,0],ij[k,1]] = data[k]

          - csc_matrix(data, (row, ptr), [(M,N)])
            ??
    """
    # Perhaps we should split up the different calling conventions to this
    # __init__ function somehow.  This could improve robustness and remove
    # some code duplication.
    def __init__(self, arg1, arg2=None, arg3=None, nzmax=100, dtype='d', copy=False):
        spmatrix.__init__(self)
        
        if isspmatrix(arg1):
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
                func = getattr(sparsetools,s.ftype+'transp')
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
                assert len(arg1) == 2 and type(arg1[0]) == int and type(arg1[1]) == int
            except AssertionError:
                raise TypeError, "matrix dimensions must be a tuple of two integers"
            (M, N) = arg1
            self.data = zeros((nzmax,), dtype)
            self.rowind = zeros((nzmax,),'i')
            self.indptr = zeros((N+1,),'i')
            self.shape = (M, N)
        elif isinstance(arg1, ArrayType) or type(arg1) == list:
            s = asarray(arg1)
            if s.dtypechar not in 'fdFD':
                # Use a double array as the source (but leave it alone)
                s = s*1.0            
            if (rank(s) == 2):  # converting from a full array
                M, N = s.shape
                dtype = s.dtypechar
                func = getattr(sparsetools,_transtabl[dtype]+'fulltocsc')
                ierr = irow = jcol = 0
                nnz = sum(ravel(s != 0.0))
                a = zeros((nnz,), dtype)
                rowa = zeros((nnz,),'i')
                ptra = zeros((N+1,),'i')
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
                self.shape = (M,N)
            elif isinstance(arg2, ArrayType) and (rank(arg2) == 2) and (shape(arg2) == (len(s),2)):
                ij = arg2
                if arg3 != None:
                    try:
                        M, N = arg3
                    except TypeError:
                        raise TypeError, "argument 3 must be a pair (M, N) of dimensions"
                else:
                    M = N = None
                temp = coo_matrix(s, ij, M=M, N=N, nzmax=nzmax, dtype=dtype)
                temp = temp.tocsc()
                self.shape = temp.shape
                self.data = temp.data
                self.rowind = temp.rowind
                self.indptr = temp.indptr
            elif type(arg2) == tuple and len(arg2)==2:
                self.data = asarray(s)
                self.rowind = arg2[0]
                self.indptr = arg2[1]
                if arg3 != None:
                    try:
                        M, N = arg3
                    except TypeError:
                        raise TypeError, "argument 3 must be a pair (M, N) of dimensions"
                else:
                    M = N = None
                if M is None:
                    try:                        
                        # we cast this to an int so type checking works
                        M = int(amax(self.rowind)) + 1
                    except ValueError:
                        M = 0
                if N is None:
                    N = len(self.indptr) - 1
                    if N == -1: N = 0
                self.shape = (M, N)
            else:
                raise ValueError, "Unrecognized form for csc_matrix constructor."
        else:
            raise ValueError, "Unrecognized form for csc_matrix constructor."

        self._check()

    def _check(self):
        M,N = self.shape
        nnz = self.indptr[-1]
        nzmax = len(self.rowind)

        if (rank(self.data) != 1) or (rank(self.rowind) != 1) or \
           (rank(self.indptr) != 1):
            raise ValueError, "Data, rowind, and indptr arrays "\
                  "should be rank 1."
        if (len(self.data) != nzmax):
            raise ValueError, "Data and row list should have same length"
        if (len(self.indptr) != N+1):
            raise ValueError, "Index pointer should be of of size N+1"
        if (nzmax < nnz):
            raise ValueError, "Nzmax must not be less than nnz."
        if (nnz>0) and (max(self.rowind[:nnz]) >= M):
            raise ValueError, "Row-values must be < M."
        if (self.indptr[-1] > len(self.rowind)):
            raise ValueError, \
                  "Last value of index list should be less than "\
                  "the size of data list"
        self.nnz = nnz
        self.nzmax = nzmax
        self._dtypechar = self.data.dtypechar
        if self._dtypechar not in 'fdFD':
            self.data = 1.0 * self.data
            self._dtypechar = self.data.dtypechar
        self.ftype = _transtabl[self._dtypechar]
        

    def __add__(self, other):
        if isscalar(other):
            raise NotImplementedError('adding a scalar to a sparse matrix is not yet supported')
        elif isspmatrix(other):
            ocs = other.tocsc()
            if (ocs.shape != self.shape):
                raise ValueError, "Inconsistent shapes."
            dtypechar = _coerce_rules[(self._dtypechar, ocs._dtypechar)]
            nnz1, nnz2 = self.nnz, other.nnz
            data1, data2 = _convert_data(self.data[:nnz1], ocs.data[:nnz2], dtypechar)
            func = getattr(sparsetools,_transtabl[dtypechar]+'cscadd')
            c,rowc,ptrc,ierr = func(data1, self.rowind[:nnz1], self.indptr, data2, ocs.rowind[:nnz2], ocs.indptr)
            if ierr:
                raise ValueError, "Ran out of space (but shouldn't have happened)."
            M, N = self.shape
            return csc_matrix(c, (rowc, ptrc), (M, N))

    def __mul__(self, other):  # implement matrix multiplication and matrix-vector multiplication
        if isspmatrix(other):
            return self.matmat(other)
        elif isscalar(other):
            new = self.copy()
            new.data *= other
            new._dtypechar = new.data.dtypechar
            new.ftype = _transtabl[new._dtypechar]
            return new
        else:
            return self.matvec(other)

    def __rmul__(self, other):  # other * self
        if isspmatrix(other):
            ocs = other.tocsc()
            return ocs.matmat(self)
        elif isscalar(other):
            new = self.copy()
            new.data = other * new.data
            new._dtypechar = new.data.dtypechar
            new.ftype = _transtabl[new._dtypechar]
            return new
        else:
            return transpose(self.rmatvec(transpose(other),conj=0))

    def __neg__(self):
        new = self.copy()
        new.data *= -1
        return new
        
    def __sub__(self, other):
        if isscalar(other):
            raise NotImplementedError('adding a scalar to a sparse matrix is not yet supported')
        elif isspmatrix(other):
            ocs = other.tocsc()
            if (ocs.shape != self.shape):
                raise ValueError, "Inconsistent shapes."
            dtypechar = _coerce_rules[(self._dtypechar,ocs._dtypechar)]
            data1, data2 = _convert_data(self.data, ocs.data, dtypechar)
            func = getattr(sparsetools,_transtabl[dtypechar]+'cscadd')
            c,rowc,ptrc,ierr = func(data1,self.rowind,self.indptr,-data2,ocs.rowind,ocs.indptr)
            if ierr:
                raise ValueError, "Ran out of space (but shouldn't have happened)."
            M, N = self.shape
            return csc_matrix(c, (rowc, ptrc), (M, N))


    def __rsub__(self, other):  # implement other - self
        if isscalar(other):
            raise NotImplementedError('adding a scalar to a sparse matrix is not yet supported')
        elif isspmatrix(other):
            ocs = other.tocsc()
            if (ocs.shape != self.shape):
                raise ValueError, "Inconsistent shapes."
            dtypechar = _coerce_rules[(self._dtypechar,ocs._dtypechar)]
            data1, data2 = _convert_data(self.data, ocs.data, dtypechar)
            func = getattr(sparsetools,_transtabl[dtypechar]+'cscadd')
            c,rowc,ptrc,ierr = func(-data1,self.rowind,self.indptr,data2,ocs.rowind,ocs.indptr)
            if ierr:
                raise ValueError, "Ran out of space (but shouldn't have happened)."
            M, N = self.shape
            return csc_matrix(c, (rowc, ptrc), (M, N))
        
    def __pow__(self, other):  
        """ Element-by-element power (unless other is a scalar, in which
        case return the matrix power.)
        """
        if isscalar(other):
            new = self.copy()
            new.data = new.data ** other
            new._dtypechar = new.data.dtypechar
            new.ftype = _transtabl[new._dtypechar]
            return new
        else:
            ocs = other.tocsc()
            if (ocs.shape != self.shape):
                raise ValueError, "Inconsistent shapes."
            dtypechar = _coerce_rules[(self._dtypechar,ocs._dtypechar)]
            nnz1, nnz2 = self.nnz, ocs.nnz
            data1, data2 = _convert_data(self.data[:nnz1], ocs.data[:nnz2], dtypechar)
            func = getattr(sparsetools,_transtabl[dtypechar]+'cscmul')
            c,rowc,ptrc,ierr = func(data1,self.rowind[:nnz1],self.indptr,data2,ocs.rowind[:nnz2],ocs.indptr)
            if ierr:
                raise ValueError, "Ran out of space (but shouldn't have happened)."
            M, N = self.shape
            return csc_matrix(c, (rowc, ptrc), (M, N))

    def transpose(self, copy=False):
        M,N = self.shape
        new = csr_matrix((N,M), nzmax=0, dtype=self._dtypechar)
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
    
    def matvec(self, x):
        if (rank(x) != 1) or (len(x) != self.shape[1]):
            raise ValueError, "Dimension mismatch"
        self._check()  # be sure everything is as it should be
        func = getattr(sparsetools,self.ftype+'cscmux')
        y = func(self.data, self.rowind, self.indptr, x, self.shape[0])
        return y

    def rmatvec(self, x, conj=1):
        if (rank(x) != 1) or (len(x) != self.shape[0]):
            raise ValueError, "Dimension mismatch"
        self._check()  # be sure everything is as it should be4
        func = getattr(sparsetools,self.ftype+'csrmux')
        if conj: cd = conj(self.data)
        else: cd = self.data
        y = func(cd, self.rowind, self.indptr, x)
        return y

    def matrixmultiply(self, other):
        """ A generic interface for matrix-matrix or matrix-vector
        multiplication.
        """
        M, K1 = self.shape
        K2, N = other.shape
        if N == 1:
            return self.matvec(other)
        else:
            return self.matmat(other)

    def matmat(self, bmat):
        self._check()
        M,K1 = self.shape
        K2,N = bmat.shape
        if (K1 != K2):
            raise ValueError, "Shape mismatch error."
        a, rowa, ptra = self.data, self.rowind, self.indptr
        if isinstance(bmat,csr_matrix):
            bmat._check()
            dtypechar = _coerce_rules[(self._dtypechar,bmat._dtypechar)]
            ftype = _transtabl[dtypechar]            
            func = getattr(sparsetools,ftype+'cscmucsr')
            b = bmat.data
            rowb = bmat.colind
            ptrb = bmat.indptr
        elif isinstance(bmat,csc_matrix):
            bmat._check()
            dtypechar = _coerce_rules[(self._dtypechar,bmat._dtypechar)]
            ftype = _transtabl[dtypechar]                        
            func = getattr(sparsetools,ftype+'cscmucsc')
            b = bmat.data
            rowb = bmat.rowind
            ptrb = bmat.indptr
        else:
            bmat = bmat.tocsc()
            dtypechar = _coerce_rules[(self._dtypechar,bmat._dtypechar)]
            ftype = _transtabl[dtypechar]                        
            func = getattr(sparsetools,ftype+'cscmucsc')
            b = bmat.data
            rowb = bmat.rowind
            ptrb = bmat.indptr
        a, b = _convert_data(a, b, dtypechar)
        newshape = (M,N)
        ptrc = zeros((N+1,),'i')
        nnzc = 2*max(ptra[-1],ptrb[-1])
        c = zeros((nnzc,),dtypechar)
        rowc = zeros((nnzc,),'i')
        ierr = irow = kcol = 0
        while 1:
            c, rowc, ptrc, irow, kcol, ierr = func(M,a,rowa,ptra,b,rowb,ptrb,c,rowc,ptrc,irow,kcol, ierr)
            if (ierr==0): break
            # otherwise we were too small and must resize
            #  calculations continue where they left off...
            percent_to_go = 1- (1.0*kcol) / N
            newnnzc = int(ceil((1+percent_to_go)*nnzc))
            c = resize1d(c,newnnzc)
            rowc = resize1d(rowc,newnnzc)
            nnzc = newnnzc

        return csc_matrix(c, (rowc, ptrc), (M, N))
            

    def __getitem__(self, key):
        if isinstance(key,types.TupleType):
            row = key[0]
            col = key[1]
            func = getattr(sparsetools,self.ftype+'cscgetel')
            M, N = self.shape
            if not (0<=row<M) or not (0<=col<N):
                raise KeyError, "Index out of bounds."
            ind, val = func(self.data, self.rowind, self.indptr, row, col)
            return val
        #elif isinstance(key,type(3)):
        elif type(key) == int:
            return self.data[key]
        else:
            raise NotImplementedError

    def __setitem__(self, key, val):
        if isinstance(key,types.TupleType):
            row = key[0]
            col = key[1]
            func = getattr(sparsetools,self.ftype+'cscsetel')
            M, N = self.shape
            if (row < 0):
                row = M + row
            if (col < 0):
                col = N + col
            if (row < 0) or (col < 0):
                raise IndexError, "Index out of bounds."
            if (col >= N):
                self.indptr = resize1d(self.indptr, col+2)
                self.indptr[N+1:] = self.indptr[N]
                N = col+1
            if (row >= M):
                M = row+1
            self.shape = (M,N)
            nzmax = self.nzmax
            if (nzmax < self.nnz+1):  # need more room 
                alloc = max(1,self.allocsize)
                self.data = resize1d(self.data, nzmax + alloc)
                self.rowind = resize1d(self.rowind, nzmax + alloc)
            func(self.data, self.rowind, self.indptr, row, col, val)
            self._check()
        elif isinstance(key, types.IntType):
            if (key < self.nnz):
                self.data[key] = val
            else:
                raise KeyError, "Key out of bounds."
        else:
            raise NotImplementedError
            
    def rowcol(self, ind):
        row = self.rowind[ind]
        col = searchsorted(self.indptr,ind+1)-1
        return (row, col)

    def getdata(self, ind):
        return self.data[ind]
    
    def tocsc(self,copy=False):
        if copy:
            new = self.copy()
        else:
            new = self
        return new

    def tocoo(self,copy=False):
        func = getattr(sparsetools,self.ftype+"csctocoo")
        data,row,col = func(self.data, self.rowind,self.indptr)
        return coo_matrix(data, (row, col), M=self.shape[0], N=self.shape[1])

    def tocsr(self,copy=False):
        return self.tocoo().tocsr()

    def todense(self):
        func = getattr(sparsetools, self.ftype+'csctofull')
        return func(self.shape[0],self.data,self.rowind,self.indptr)

    #  should add full option to eliminate non-zero entries.
    def prune(self):
        nnz = self.indptr[-1]
        if self.nzmax <= nnz:
            if self.nzmax < nnz:
                raise RunTimeError, "Should never have nnz > nzmax"            
            return
        self.nnz = nnz
        self.data = self.data[:nnz]
        self.rowind = self.rowind[:nnz]
        self.nzmax = nnz
        self._check()

    def copy(self):
        dtype = self._dtypechar
        new = csc_matrix(self.shape, nzmax=0, dtype=dtype)
        new.data = self.data.copy()
        new.rowind = self.rowind.copy()
        new.indptr = self.indptr.copy()
        new._check()
        return new
    
# compressed sparse row matrix
# 
class csr_matrix(spmatrix):
    """ Compressed sparse row matrix
        This can be instantiated in two ways:
          - csr_matrix(s)
            with another sparse matrix s (sugar for .tocsr())

          - csr_matrix((M,N), [nzmax, dtype])
            to construct a container, where (M,N) are dimensions and
            nzmax, dtype are optional, defaulting to nzmax=100 and dtype='d'.

          - csr_matrix(data, ij, [(M,N),nzmax])
            where data, ij satisfy:
                a[ij[k,0],ij[k,1]] = data[k]

          - csr_matrix(data, (row, ptr), [(M,N)])
            ??
    """
    # Perhaps we should split up the different calling conventions to this
    # __init__ function somehow.  This could improve robustness and remove
    # some code duplication.
    def __init__(self, arg1, arg2=None, arg3=None, nzmax=100, dtype='d', copy=False):
        spmatrix.__init__(self)
        if isspmatrix(arg1):
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
                func = getattr(sparsetools,s.ftype+'transp')
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
                assert len(arg1) == 2 and type(arg1[0]) == int and type(arg1[1]) == int
            except AssertionError:
                raise TypeError, "matrix dimensions must be a tuple of two integers"
            (M, N) = arg1
            self.data = zeros((nzmax,), dtype)
            self.colind = zeros((nzmax,), 'i')
            self.indptr = zeros((M+1,), 'i')
            self.shape = (M, N)
        elif isinstance(arg1, ArrayType) or type(arg1) == list:
            s = asarray(arg1)
            if (rank(s) == 2):  # converting from a full array
                ocsc = csc_matrix(transpose(s))
                self.colind = ocsc.rowind
                self.indptr = ocsc.indptr
                self.data = ocsc.data
                self.shape = (ocsc.shape[1], ocsc.shape[0])

            elif isinstance(arg2, ArrayType) and (rank(arg2) == 2) and (shape(arg2) == (len(s),2)):
                ij = arg2
                ijnew = ij.copy()
                ijnew[:,0] = ij[:,1]
                ijnew[:,1] = ij[:,0]
                temp = coo_matrix(s,ijnew,M=M,N=N,nzmax=nzmax,
                                  dtype=dtype)
                temp = temp.tocsc()
                self.shape = temp.shape
                self.data = temp.data
                self.colind = temp.colind
                self.indptr = temp.indptr
            elif type(arg2) == tuple and len(arg2)==2:
                self.data = asarray(s)
                self.colind = arg2[0]
                self.indptr = arg2[1]
                if arg3 != None:
                    try:
                        M, N = arg3
                    except TypeError:
                        raise TypeError, "argument 3 must be a pair (M, N) of dimensions"
                else:
                    M = N = None
                if N is None:
                    try:
                        # we cast this to an int so type checking works
                        N = int(amax(self.colind)) + 1
                    except ValueError:
                        N = 0
                if M is None:
                    M = len(self.indptr) - 1
                    if M == -1: M = 0
                self.shape = (M,N)
            else:
                raise ValueError, "Unrecognized form for csr_matrix constructor."
        else:
            raise ValueError, "Unrecognized form for csr_matrix constructor."

        self._check()

    def _check(self):
        M, N = self.shape
        nnz = self.indptr[-1]
        nzmax = len(self.colind)
        if (rank(self.data) != 1) or (rank(self.colind) != 1) or \
           (rank(self.indptr) != 1):
            raise ValueError, "Data, colind, and indptr arrays "\
                  "should be rank 1."
        if (len(self.data) != nzmax):
            raise ValueError, "Data and row list should have same length"
        if (len(self.indptr) != M+1):
            raise ValueError, "Index pointer should be of length #rows + 1"
        if (nnz>0) and (max(self.colind[:nnz]) >= N):
            raise ValueError, "Column-values must be < N."
        if (nnz > nzmax):
            raise ValueError, \
                  "Last value of index list should be less than "\
                  "the size of data list"
        self.nnz = nnz
        self.nzmax = nzmax
        self._dtypechar = self.data.dtypechar
        if self._dtypechar not in 'fdFD':
            self.data = self.data + 0.0
            self._dtypechar = self.data.dtypechar
            
        self.ftype = _transtabl[self._dtypechar]

        
    def __add__(self, other):
        # First check if argument is a scalar
        if isscalar(other):
            # Now we would add this scalar to every element.
            raise NotImplementedError('adding a scalar to a sparse matrix is not yet supported')
        ocs = other.tocsr()
        if (ocs.shape != self.shape):
            raise ValueError, "inconsistent shapes."

        dtypechar = _coerce_rules[(self._dtypechar, ocs._dtypechar)]
        data1, data2 = _convert_data(self.data, ocs.data, dtypechar)
        func = getattr(sparsetools,_transtabl[dtypechar]+'cscadd')
        c,colc,ptrc,ierr = func(data1,self.colind,self.indptr,data2,other.colind,other.indptr)
        if ierr:
            raise ValueError, "Ran out of space (but shouldn't have happened)."
        M, N = self.shape
        return csr_matrix(c, (colc, ptrc), (M, N))


    def __mul__(self, other):  # implement matrix multiplication and matrix-vector multiplication
        if isspmatrix(other):
            return self.matmat(other)
        elif isscalar(other):
            new = self.copy()
            new.data *= other
            new._dtypechar = new.data.dtypechar
            new.ftype = _transtabl[new._dtypechar]
            return new
        else:
            return self.matvec(other)

    def __rmul__(self, other):  # other * self
        if isspmatrix(other):
            ocs = other.tocsc()
            return occ.matmat(self)
        elif isscalar(other):
            new = self.copy()
            new.data = other * new.data
            new._dtypechar = new.data.dtypechar
            new.ftype = _transtabl[new._dtypechar]
            return new
        else:
            return transpose(self.rmatvec(transpose(other),conj=0))

    def __neg__(self):
        new = self.copy()
        new.data *= -1
        return new
        
    def __sub__(self, other):
        # First check if argument is a scalar
        if isscalar(other):
            # Now we would add this scalar to every element.
            raise NotImplementedError('subtracting a scalar from a sparse matrix is not yet supported')
        elif isspmatrix(other):
            ocs = other.tocsr()
            if (ocs.shape != self.shape):
                raise ValueError, "Inconsistent shapes."
            dtypechar = _coerce_rules[(self._dtypechar, ocs._dtypechar)]
            data1, data2 = _convert_data(self.data, ocs.data, dtypechar)
            func = getattr(sparsetools,_transtabl[dtypechar]+'cscadd')
            c,colc,ptrc,ierr = func(data1,self.colind,self.indptr,-data2,other.colind,other.indptr)
            if ierr:
                raise ValueError, "Ran out of space (but shouldn't have happened)."
            M, N = self.shape
            return csr_matrix(c, (colc, ptrc), (M, N))


    def __rsub__(self, other):  # implement other - self
        ocs = other.tocsr()
        if (ocs.shape != self.shape):
            raise ValueError, "Inconsistent shapes."
        dtypechar = _coerce_rules[(self._dtypechar, ocs._dtypechar)]
        data1, data2 = _convert_data(self.data, ocs.data, dtypechar)
        func = getattr(sparsetools,_transtabl[dtypechar]+'cscadd')
        c,colc,ptrc,ierr = func(-data1,self.colind,self.indptr,data2,other.colind,other.indptr)
        if ierr:
            raise ValueError, "Ran out of space (but shouldn't have happened)."
        M, N = self.shape
        return csc_matrix(c, (colc, ptrc), (M, N))
        
    def __pow__(self, other):  
        """ Element-by-element power (unless other is a scalar, in which
        case return the matrix power.)
        """
        if isscalar(other):
            new = self.copy()
            new.data = new.data ** other
            new._dtypechar = new.data.dtypechar
            new.ftype = _transtabl[new._dtypechar]
            return new
        else:
            ocs = other.tocsr()
            if (ocs.shape != self.shape):
                raise ValueError, "Inconsistent shapes."
            dtypechar = _coerce_rules[(self._dtypechar, ocs._dtypechar)]
            data1, data2 = _convert_data(self.data, ocs.data, dtypechar)
            func = getattr(sparsetools,_transtabl[dtypechar]+'cscmul')
            c,colc,ptrc,ierr = func(data1,self.colind,self.indptr,data2, ocs.colind, ocs.indptr)
            if ierr:
                raise ValueError, "Ran out of space (but shouldn't have happened)."
            M, N = self.shape
            return csr_matrix(c, (colc, ptrc), (M, N))

    def transpose(self, copy=False):
        M,N = self.shape
        new = csc_matrix(N,M,nzmax=0,dtypechar=self._dtypechar)
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

    def matvec(self, x):
        if (rank(x) != 1) or (len(x) != self.shape[1]):
            raise ValueError, "Dimension mismatch"
        self._check()  # be sure everything is as it should be
        func = getattr(sparsetools,self.ftype+'csrmux')
        y = func(self.data, self.colind, self.indptr, x)
        return y

    def rmatvec(self, x, conj=1):
        if (rank(x) != 1) or (len(x) != self.shape[0]):
            raise ValueError, "Dimension mismatch"
        self._check()  # be sure everything is as it should be4
        func = getattr(sparsetools,self.ftype+'cscmux')
        if conj: cd = conj(self.data)
        else: cd = self.data
        y = func(cd, self.colind, self.indptr, x, self.shape[1])
        return y

    def matrixmultiply(self, other):
        """ A generic interface for matrix-matrix or matrix-vector
        multiplication.
        """
        M, K1 = self.shape
        K2, N = other.shape
        if N == 1:
            return self.matvec(other)
        else:
            return self.matmat(other)

    def matmat(self, bmat):
        self._check()
        M,K1 = self.shape
        K2,N = bmat.shape
        a, rowa, ptra = self.data, self.colind, self.indptr
        if (K1 != K2):
            raise ValueError, "Shape mismatch error."
        if isinstance(bmat,csc_matrix):            
            bmat._check()
            dtypechar = _coerce_rules[(self._dtypechar,bmat._dtypechar)]
            ftype = _transtabl[dtypechar]            
            func = getattr(sparsetools,ftype+'csrmucsc')
            b = bmat.data
            colb = bmat.rowind
            ptrb = bmat.indptr
            out = 'csc'
            firstarg = ()
        elif isinstance(bmat,csr_matrix):
            bmat._check()
            dtypechar = _coerce_rules[(self._dtypechar,bmat._dtypechar)]
            ftype = _transtabl[dtypechar]            
            func = getattr(sparsetools,ftype+'cscmucsc')
            b, colb, ptrb = a, rowa, ptra
            a, rowa, ptra = bmat.data, bmat.colind, bmat.indptr
            out = 'csr'
            firstarg = (N,)
        else:
            bmat = bmat.tocsc()
            dtypechar = _coerce_rules[(self._dtypechar,bmat._dtypechar)]
            ftype = _transtabl[dtypechar]            
            func = getattr(sparsetools,ftype+'csrmucsc')
            b = bmat.data
            colb = bmat.colind
            ptrb = bmat.indptr
            out = 'csc'
            firstarg = ()
        a, b = _convert_data(a, b, dtypechar)            
        newshape = (M,N)
        if out == 'csr':
            ptrc = zeros((M+1,),'i')
        else:
            ptrc = zeros((N+1,),'i')
        nnzc = 2*max(ptra[-1],ptrb[-1])
        c = zeros((nnzc,),dtypechar)
        rowc = zeros((nnzc,),'i')
        ierr = irow = kcol = 0
        while 1:
            args = firstarg+(a,rowa,ptra,b,colb,ptrb,c,rowc,ptrc,irow,
                             kcol, ierr)
            c, rowc, ptrc, irow, kcol, ierr = func(*args)
            if (ierr==0): break
            # otherwise we were too small and must resize
            percent_to_go = 1- (1.0*kcol) / N
            newnnzc = int(ceil((1+percent_to_go)*nnzc))
            c = resize1d(c,newnnzc)
            rowc = resize1d(rowc,newnnzc)
            nnzc = newnnzc

        if out == 'csr':
            # FIXME
            # Is this correct??  Does rowc here mean colc?
            return csr_matrix(c, (rowc, ptrc), (M, N))
        else:
            return csc_matrix(c, (rowc, ptrc), (M, N))

    def __getitem__(self, key):
        if isinstance(key,types.TupleType):
            row = key[0]
            col = key[1]
            func = getattr(sparsetools,self.ftype+'cscgetel')
            M, N = self.shape
            if (row < 0):
                row = M + row
            if (col < 0):
                col = N + col
            if (row >= M ) or (col >= N) or (row < 0) or (col < 0):
                raise IndexError, "Index out of bounds."
            ind, val = func(self.data, self.colind, self.indptr, col, row)
            return val
        #elif isinstance(key,type(3)):
        elif type(key) == int:
            return self.data[key]
        else:
            raise NotImplementedError


    def __setitem__(self, key, val):
        if isinstance(key,types.TupleType):
            row = key[0]
            col = key[1]
            func = getattr(sparsetools,self.ftype+'cscsetel')
            M, N = self.shape
            if (row < 0):
                row = M + row
            if (col < 0):
                col = N + col
            if (row < 0) or (col < 0):
                raise KeyError, "Index out of bounds."
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
            self.shape = (M,N)
            nzmax = self.nzmax
            if (nzmax < self.nnz+1):  # need more room 
                alloc = max(1,self.allocsize)
                self.data = resize1d(self.data, nzmax + alloc)
                self.colind = resize1d(self.colind, nzmax + alloc)
            func(self.data, self.colind, self.indptr, col, row, val)
            self._check()
        elif isinstance(key, types.IntType):
            if (key < self.nnz):
                self.data[key] = val
            else:
                raise KeyError, "Key out of bounds."
        else:
            raise NotImplementedError
            
    def rowcol(self, ind):
        col = self.colind[ind]
        row = searchsorted(self.indptr,ind+1)-1
        return (row, col)

    def getdata(self, ind):
        return self.data[ind]

    def tocsr(self,copy=False):
        if copy:
            new = self.copy()
        else:
            new = self
        return new

    def tocoo(self,copy=False):
        func = getattr(sparsetools,self.ftype+"csctocoo")
        data,col,row = func(self.data, self.colind,self.indptr)
        return coo_matrix(data, (row, col), M=self.shape[0], N=self.shape[1])

    def tocsc(self,copy=False):
        return self.tocoo().tocsc()

    def todense(self):
        func = getattr(sparsetools, self.ftype+'csctofull')
        s = func(self.shape[1],self.data,self.colind,self.indptr)
        return transpose(s)

    #  should add full option to eliminate non-zero entries.
    def prune(self):
        nnz = self.indptr[-1]
        if self.nzmax <= nnz:
            if self.nzmax < nnz:
                raise RunTimeError, "Should never have nnz > nzmax"            
            return
        self.nnz = nnz
        self.data = self.data[:nnz]
        self.colind = self.colind[:nnz]
        self.nzmax = nnz
        self._check()

    def copy(self):
        new = csr_matrix(self.shape, nzmax=0, dtype=self._dtypechar)
        new.data = self.data.copy()
        new.colind = self.colind.copy()
        new.indptr = self.indptr.copy()
        new._check()
        return new

#   A simple "dictionary-based" sparse matrix.
#   keys must be 2-tuples of any object.  The object must have __int__
#    defined to return integer
#   should also define __cmp__ as well to compare the integer-based keys

def csc_cmp(x,y):
    if (x == y): return 0
    elif (x[1] == y[1]):
        if (x[0] > y[0]): return 1
        elif (x[0] == y[0]): return 0
        else: return -1
    elif (x[1] > y[1]): return 1
    else: return -1
            
# dictionary of keys based matrix
class dok_matrix(spmatrix, dict):
    """ A dictionary of keys based matrix.  This is slow: does type
        checking on input and uses dicts, but it is efficient for constructing
        sparse matrices for conversion to other sparse matrix types.
    """
    def __init__(self, A=None):
        """ Create a new dictionary-of-keys sparse matrix.  An optional
        argument A is accepted, which initializes the dok_matrix with it.
        (For now this only supports dense matrices.)
        """
        dict.__init__(self)
        spmatrix.__init__(self)
        self.shape = (0,0)
        self.nnz = 0

        if A is not None:
            if isspmatrix(A):
                # For sparse matrices, this is too inefficient; we need 
                # something else.
                raise NotImplementedError, "initializing a dok_matrix with a sparse matrix is not supported yet"
            elif isdense(A):
                A = asarray(A)
                N,M = A.shape
                for n in range(N):
                    for m in range(M):
                        if A[n,m] != 0:
                            self[n,m] = A[n,m]
            else:
                raise TypeError, "argument should be a sparse or dense matrix"

    def __str__(self):
        val = ''
        nnz = len(self.keys())
        self.nnz = nnz
        keys = self.keys()
        keys.sort()
        if nnz > self.maxprint:
            for k in xrange(self.maxprint / 2):
                key = keys[k]
                val += "  %s\t%s\n" % (str(key),str(self[key]))
            val = val + "   :    \t  :\n"
            for k in xrange(nnz-self.maxprint/2,nnz):
                key = keys[k]
                val += "  %s\t%s\n" % (str(key),str(self[key]))
        else:
            for k in xrange(nnz):
                key = keys[k]
                val += "  %s\t%s\n" % (str(key),str(self[key]))
        return val[:-1]

    def __repr__(self):
        nnz = self.getnnz()
        format = self.getformat()
        return "<%dx%d sparse matrix with %d stored "\
               "elements in %s format>" % \
               (self.shape + (nnz, _formats[format][1]))

    def __getitem__(self, key):
        # Sanity checks: key must be a pair of integers
        if not isinstance(key, tuple) or len(key) != 2:
            raise TypeError, "key must be a tuple of two integers"
        if type(key[0]) != int or type(key[1]) != int:
            raise TypeError, "key must be a tuple of two integers"

        return self.get(key,0)

    def __setitem__(self, key, value):
        # Sanity checks: key must be a pair of integers
        if not isinstance(key, tuple) or len(key) != 2:
            raise TypeError, "key must be a tuple of two integers"
        if type(key[0]) != int or type(key[1]) != int:
            raise TypeError, "key must be a tuple of two integers"

        if (value == 0):
            if self.has_key(key):  # get rid of it something already there
                del self[key]              # otherwise do nothing.
            return
        if not self.has_key(key):
            self.nnz += 1
        dict.__setitem__(self, key, value)
        newrows = max(self.shape[0], int(key[0])+1)
        newcols = max(self.shape[1], int(key[1])+1)
        self.shape = (newrows, newcols)
    
    def __add__(self, other):
        # First check if argument is a scalar
        if isscalar(other):
            # Now we would add this scalar to every element.
            raise NotImplementedError('adding a scalar to a sparse matrix is not yet supported')
        elif isinstance(other, dok_matrix):
            res = dok_matrix()
            res.update(self)
            res.shape = self.shape
            res.nnz = self.nnz
            for key in other.keys():
                res[key] += other[key]
        else:
            csc = self.tocsc()
            res = csc + other
        return res

    def __sub__(self, other):
        # First check if argument is a scalar
        if isscalar(other):
            # Now we would add this scalar to every element.
            raise NotImplementedError('subtracting a scalar from a sparse matrix is not yet supported')
        elif isinstance(other, dok_matrix):
            res = dok_matrix()
            res.update(self)
            res.shape = self.shape
            res.nnz = self.nnz
            for key in other.keys():
                res[key] -= other[key]
        else:
            csc = self.tocsc()
            res = csc - other
        return res
    
    def __neg__(self):
        res = dok_matrix()
        for key in self.keys():
            res[key] = -self[key]
        return res

    def __mul__(self, other):
        if isspmatrix(other):
            return self.matmat(other)
        other = asarray(other)
        if rank(other) > 0:
            return self.matvec(other)
        res = dok_matrix()
        for key in self.keys():
            res[key] = other * self[key]
        return res

    def __len__(self):
        return len(self.keys())    

    def transpose(self):
        """ Return the transpose
        """
        new = dok_matrix()
        for key in self.keys():
            new[key[1],key[0]] = self[key]
        return new

    def conjtransp(self):
        """ Return the conjugate transpose
        """
        new = dok_matrix()
        for key in self.keys():
            new[key[1],key[0]] = conj(self[key])
        return new

    def copy(self):
        new = dok_matrix()
        new.update(self)
        new.nnz = self.nnz
        new.shape = self.shape
        return new
        
    def take(self, cols_or_rows, columns=1):
        # Extract columns or rows as indictated from matrix
        # assume cols_or_rows is sorted
        res = dok_matrix()
        indx = int((columns == 1))
        N = len(cols_or_rows)
        if indx: # columns
            for key in self.keys():
                num = searchsorted(cols_or_rows,key[1])
                if num < N:
                    newkey = (key[0],num)
                    res[newkey] = self[key]
        else:
            for key in self.keys():
                num = searchsorted(cols_or_rows,key[0])
                if num < N:
                    newkey = (num,key[1])
                    res[newkey] = self[key]            
        return res

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
                num = searchsorted(cols_or_rows,key[1])
                if cols_or_rows[num]==key[1]:
                    newkey = (key[0],num)
                    ext[newkey] = self[key]
                else:
                    newkey = (key[0],key[1]-num)
                    base[newkey] = self[key]
        else:
            for key in self.keys():
                num = searchsorted(cols_or_rows,key[0])
                if cols_or_rows[num]==key[0]:
                    newkey = (num,key[1])
                    ext[newkey] = self[key]
                else:
                    newkey = (key[0]-num,key[1])
                    base[newkey] = self[key]            
        return base, ext


    def matvec(self, other):
        other = asarray(other)
        if other.shape[0] != self.shape[1]:
            raise ValueError, "Dimensions do not match."
        res = [0]*self.shape[0]
        for key in self.keys():
            res[int(key[0])] += self[key] * other[int(key[1]),...]
        return array(res)        

    def rmatvec(self, other):
        other = asarray(other)
	if other.shape[-1] != self.shape[0]:
	    raise ValueError, "Dimensions do not match."
	res = [0]*self.shape[1]
	for key in self.keys():
            res[int(key[1])] += other[..., int(key[0])] * conj(self[key])
	return array(res)

    def setdiag(self, values, k=0):
        N = len(values)
        for n in range(N):
            self[n,n+k] = values[n]
        return

    def tocsr(self):
        """ Return Compressed Sparse Row format arrays for this matrix
        """
        keys = self.keys()
        keys.sort()
        nnz = self.nnz
        assert nnz == len(keys)
        data = [0]*nnz
        colind = [0]*nnz
        row_ptr = [0]*(self.shape[0]+1)
        current_row = 0
        k = 0
        for key in keys:
            ikey0 = int(key[0])
            ikey1 = int(key[1])
            if ikey0 != current_row:
                N = ikey0-current_row
                row_ptr[current_row+1:ikey0+1] = [k]*N
                current_row = ikey0
            data[k] = self[key]
            colind[k] = ikey1
            k += 1
        row_ptr[-1] = nnz
        data = array(data)
        colind = array(colind)
        row_ptr = array(row_ptr)
        return csr_matrix(data, (colind, row_ptr))

    def tocsc(self):
        """ Return Compressed Sparse Column format arrays for this matrix
        """
        keys = self.keys()
        #  Sort based on columns
        keys.sort(csc_cmp)
        nnz = self.nnz
        assert nnz == len(keys)
        data = [0]*nnz
        rowind = [0]*nnz
        col_ptr = [0]*(self.shape[1]+1)
        current_col = 0
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
        col_ptr[-1] = nnz
        data = array(data)
        rowind = array(rowind)
        col_ptr = array(col_ptr)
        return csc_matrix(data, (rowind, col_ptr))

    def todense(self,dtypechar=None):
        if dtypechar is None:
            dtypechar = 'd'
        new = zeros(self.shape, dtypechar)
        for key in self.keys():
            ikey0 = int(key[0])
            ikey1 = int(key[1])
            new[ikey0,ikey1] = self[key]
        if max(ravel(abs(new.imag))) == 0:
            new = new.real
        return new
    

# dictionary of dictionaries based matrix
class dod_matrix(spmatrix):
    pass

# linked list matrix
class lnk_matrix(spmatrix):
    pass

# coordinate lists format
#   a[ij[k][0],ij[k][1]] = obj[k]
# 
class coo_matrix(spmatrix):
    def __init__(self, obj, ij, M=None, N=None, nzmax=None, dtype=None):
        spmatrix.__init__(self)
        if type(ij) is type(()) and len(ij)==2:
            if M is None:
                # we cast this to an int so type checking works
                M = int(amax(ij[0]))
            if N is None:
                N = int(amax(ij[1]))
            self.row = asarray(ij[0],'i')
            self.col = asarray(ij[1],'i')
        else:
            aij = asarray(ij,'i')
            if M is None:
                M = int(amax(aij[:,0]))
            if N is None:
                N = int(amax(aij[:,1]))
            self.row = aij[:,0]
            self.col = aij[:,1]
        aobj = asarray(obj,dtype=dtype)
        self.shape = (M,N)
        if nzmax is None:
            nzmax = len(aobj)
        self.nzmax = nzmax
        self.data = aobj
        self._dtypechar = aobj.dtypechar
        self._check()

    def _check(self):
        nnz = len(self.data)
        if (nnz != len(self.row)) or (nnz != len(self.col)):
            raise ValueError, "Row, column, and data array must all be "\
                  "the same length."
        if (self.nzmax < nnz):
            raise ValueError, "nzmax must be >= nnz"
        self.nnz = nnz
        self.ftype = _transtabl[self._dtypechar]

    def _normalize(self, rowfirst = 0):
        if rowfirst:
            import itertools
            l = zip(self.row,self.col,self.data)
            l.sort()
            row,col,data = list(itertools.izip(*l))
            return data, row, col
        if getattr(self,'_is_normalized',None):
            return self.data, self.row, self.col
        import itertools
        l = zip(self.col,self.row,self.data)
        l.sort()
        col,row,data = list(itertools.izip(*l))
        self.col = asarray(col,'i')
        self.row = asarray(row,'i')
        self.data = array(data,self._dtypechar)
        setattr(self,'_is_normalized',1)
        return self.data, self.row, self.col

    def rowcol(self, num):
        return (self.row[num], self.col[num])

    def getdata(self, num):
        return self.data[num]

    def tocsc(self):
        func = getattr(sparsetools,self.ftype+"cootocsc")
        data, row, col = self._normalize()
        a, rowa, ptra, ierr = func(self.shape[1], data, row, col)
        if ierr:
            raise RuntimeError, "Error in conversion."
        return csc_matrix(a, (rowa, ptra), self.shape)

    def tocsr(self):
        func = getattr(sparsetools,self.ftype+"cootocsc")
        data,row,col = self._normalize(rowfirst=1)
        a, cola, ptra, ierr = func(self.shape[0], data, col, row)
        if ierr:
            raise RuntimeError, "Error in conversion."
        return csr_matrix(a, (cola, ptra), self.shape)

    def tocoo(self,copy=False):
        if copy:
            new = self.copy()
        else:
            new = self
        return new
        
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

def isdense(x):
    # What's the best way to check for this?  The following fails on import:
    # import numerictypes
    # return numerictypes.isdtype(x)
    temp = zeros(1)
    return type(x) == type(temp)

def _spdiags_tosub(diag_num,a,b):
    part1 = where(less(diag_num,a),abs(diag_num-a),0)
    part2 = where(greater(diag_num,b),abs(diag_num-b),0)
    return part1+part2
                      
def spdiags(diags, offsets, M, N):
    """Return a sparse matrix in CSR format given its diagonals.
    
    B = spdiags(diags, offsets, M, N)

    Inputs:
        diags  --  rows contain diagonal values
        offsets -- diagonals to set (0 is main)
        M, N    -- sparse matrix returned is M X N
    """
    diags = array(transpose(diags),copy=True)
    if diags.dtypechar not in 'fdFD':
        diags = diags.astype('d')
    offsets = array(offsets,copy=False)
    mtype = diags.dtypechar
    assert(len(offsets) == diags.shape[1])
    # set correct diagonal to csr conversion routine for this type
    diagfunc = eval('sparsetools.'+_transtabl[mtype]+'diatocsr')
    a, rowa, ptra, ierr = diagfunc(M,N,diags,offsets)
    if ierr:
        raise ValueError, "Ran out of memory (shouldn't have happened)"
    return csc_matrix(a,(rowa,ptra), (M, N))

def solve(A,b,permc_spec=2):
    if not hasattr(A, 'tocsr') and not hasattr(A, 'tocsc'):
        raise ValueError, "Sparse matrix must be able to return CSC format--"\
              "A.tocsc()--or CSR format--A.tocsr()"
    if not hasattr(A,'shape'):
        raise ValueError, "Sparse matrix must be able to return shape (rows,cols) = A.shape"
    M, N = A.shape
    if (M != N):
        raise ValueError, "Matrix must be square."    
    if hasattr(A, 'tocsc'):
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
    return gssv(N,lastel,data,index0,index1,b,csc,permc_spec)[0]
    

def lu_factor(A, permc_spec=2, diag_pivot_thresh=1.0,
              drop_tol=0.0, relax=1, panel_size=10):
    M, N = A.shape
    if (M != N):
        raise ValueError, "Can only factor square matrices."
    csc = A.tocsc()
    gstrf = eval('_superlu.' + csc.ftype + 'gstrf')
    return gstrf(N,csc.nnz,csc.data,csc.rowind,csc.indptr,permc_spec,
                 diag_pivot_thresh, drop_tol, relax, panel_size)
        

if __name__ == "__main__":
    a = csc_matrix(arange(1,9),transpose([[0,1,1,2,2,3,3,4],[0,1,3,0,2,3,4,4]]))
    print "Representation of a matrix:"
    print repr(a)
    print "How a matrix prints."
    print a
    print "Adding two matrices."
    b = a+a
    print b
    print "Subtracting two matrices."
    c = b - a
    print c
    print "Multiplying a sparse matrix by a dense vector."
    d = a*[1,2,3,4,5]
    print d
    print [1,2,3,4,5]*a

    print "Inverting a sparse linear system."
    print "The sparse matrix (constructed from diagonals)."
    a = spdiags([[1,2,3,4,5],[6,5,8,9,10]],[0,1],5,5)
    b = array([1,2,3,4,5])
    print a
    print "Solve: single precision complex."
    a = a.astype('F')
    x = solve(a,b)
    print x
    print "Error: ",a*x-b

    print "Solve: double precision complex."
    a = a.astype('D')
    x = solve(a,b)
    print x
    print "Error: ",a*x-b

    print "Solve: double precision."
    a = a.astype('d')
    x = solve(a,b)
    print x
    print "Error: ",a*x-b

    print "Solve: single precision."
    a = a.astype('f')
    x = solve(a,b.astype('f'))
    print x
    print "Error: ",a*x-b

    print "(Various small tests follow ...)\n"
    print "Dictionary of keys matrix:"
    a = dok_matrix()
    a[1,1] = 1.
    a[1,5] = 1.
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
