from scipy_base import *
from scipy_base.fastumath import *
#from scipy_base import isscalar, rank, shape, resize, ArrayType, transpose
import types
import sparsetools
import _superlu
import sys

if sys.version[:3] < '2.3':
    True = 1
    False = 0

def resize1d(arr, newlen):
    old = len(arr)
    new = zeros((newlen,),arr.typecode())
    new[:old] = arr
    return new


MAXPRINT=50
ALLOCSIZE=1000
# The formats that we might potentially understand.

_coerce_rules = {('f','f'):'f', ('f','d'):'d', ('f','F'):'F',
                 ('f','D'):'D', ('d','f'):'d', ('d','d'):'d',
                 ('d','F'):'D', ('d','D'):'D', ('F','f'):'F',
                 ('F','d'):'D', ('F','F'):'F', ('F','D'):'D',
                 ('D','f'):'D', ('D','d'):'d', ('D','F'):'D',
                 ('D','D'):'D'}
_transtabl = {'f':'s','d':'d','F':'c','D':'z'}
_itranstabl = {'s':'f','d':'d','c':'F','z':'D'}
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
    if data1.typecode() != newtype:
        data1 = data1.astype(newtype)
    if data2.typecode() != newtype:
        data2 = data2.astype(newtype)
    return data1, data2        
    
# This class provides a base class for all sparse matrices
#   most of the work is provided by subclasses

class spmatrix:
    def __init__(self, format, maxprint=MAXPRINT, allocsize=ALLOCSIZE):
        self.format = format
        self.maxprint = maxprint
        self.allocsize = allocsize
        
    def getmaxprint(self):
        try:
            maxprint = self.maxprint
        except AttributeError:
            maxprint = MAXPRINT
        return maxprint
        
    def gettypecode(self):
        try:
            typecode = self.typecode
        except AttributeError:
            typecode = None
        return typecode

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
        typecode = self.gettypecode()
        nnz = self.getnnz()
        format = self.getformat()
        nzmax = self.getnzmax()
        return "<%dx%d sparse matrix of type '%s' with %d stored "\
               "elements (space for %d) in\n\t%s format>" % \
               (self.shape + (typecode, nnz, nzmax, _formats[format][1]))

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
        res = csc + other
        return res

    def __sub__(self, other):
        csc = self.tocsc()
        res = csc - other
        return res

    def __rsub__(self, other): # other - self
        csc = self.tocsc()
        res = csc.__rsub__(other)
        return res


    def __mul__(self, other):
        csc = self.tocsc()
        res = csc * other
        return res

    def __rmul__(self, other):
        csc = self.tocsc()
        res = csc.__rmul__(other)
        return res
        
    def __neg__(self):
        csc = self.tocsc()
        res = -csc
        return res

    def matmat(self, other):
        csc = self.tocsc()
        res = csc.matmat(other)
        return res

    def matvec(self, vec):
        csc = self.tocsc()
        res = csc.matvec(vec)
        return res

    def todense(self):
        csc = self.tocsc()
        return csc.todense()

    def tocoo(self):
        csc = self.tocsc()
        return csc.tocoo()

# compressed sparse column matrix
#  This can be instantiated in many ways
#    - with another sparse matrix (sugar for .tocsc())
#    - with M,N,nzmax,typecode  to construct a container
#    - with data, ij, {M,N,nzmax}
#           a[ij[k,0],ij[k,1]] = data[k]
#    - with data, (row, ptr)
# 

class csc_matrix(spmatrix):
    def __init__(self,s,ij=None,M=None,N=None,nzmax=100,typecode=Float,copy=0):
        spmatrix.__init__(self, 'csc')
        if isinstance(s,spmatrix):
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
        elif isinstance(s,type(3)):
            M=s
            N=ij
            self.data = zeros((nzmax,),typecode)
            self.rowind = zeros((nzmax,),'i')
            self.indptr = zeros((N+1,),'i')
            self.shape = (M,N)
        elif (isinstance(s,ArrayType) or \
              isinstance(s,type([]))):
            s = asarray(s)
            if (rank(s) == 2):  # converting from a full array
                M, N = s.shape
                s = asarray(s)
                if s.typecode() not in 'fdFD':
                    s = s*1.0
                typecode = s.typecode()
                func = getattr(sparsetools,_transtabl[typecode]+'fulltocsc')
                ierr = irow = jcol = 0
                nnz = sum(ravel(s != 0.0))
                a = zeros((nnz,),typecode)
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
            elif isinstance(ij, ArrayType) and (rank(ij) == 2) and (shape(ij) == (len(s),2)):
                temp = coo_matrix(s,ij,M=M,N=N,nzmax=nzmax,typecode=typecode)
                temp = temp.tocsc()
                self.data = temp.data
                self.rowind = temp.rowind
                self.indptr = temp.indptr
                self.shape = temp.shape
            elif isinstance(ij, types.TupleType) and (len(ij)==2):
                self.data = asarray(s)
                self.rowind = ij[0]
                self.indptr = ij[1]
                if M is None:
                    try:                        
                        M = amax(self.rowind) + 1
                    except ValueError:
                        M = 0
                if N is None:
                    N = len(self.indptr) - 1
                    if N == -1: N = 0
                self.shape = (M,N)
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
        self.typecode = self.data.typecode()
        if self.typecode not in 'fdFD':
            self.data = self.data.astype('d')
            self.typecode = 'd'
        self.ftype = _transtabl[self.typecode]
        

    def __add__(self, other):
        ocs = csc_matrix(other)
        if (ocs.shape != self.shape):
            raise ValueError, "Inconsistent shapes."
        typecode = _coerce_rules[(self.typecode,other.typecode)]
        nnz1, nnz2 = self.nnz, other.nnz
        data1, data2 = _convert_data(self.data[:nnz1], ocs.data[:nnz2], typecode)
        func = getattr(sparsetools,_transtabl[typecode]+'cscadd')
        c,rowc,ptrc,ierr = func(data1,self.rowind[:nnz1],self.indptr,data2,ocs.rowind[:nnz2],ocs.indptr)
        if ierr:
            raise ValueError, "Ran out of space (but shouldn't have happened)."
        M, N = self.shape
        return csc_matrix(c,(rowc,ptrc),M=M,N=N)

    def __mul__(self, other):  # implement matrix multiplication and matrix-vector multiplication
        if isspmatrix(other):
            return self.matmat(other)
        elif isscalar(other):
            new = self.copy()
            new.data = new.data * other
            new.typecode = new.data.typecode()
            new.ftype = _transtabl[new.typecode]
            return new
        else:
            return self.matvec(other)

    def __rmul__(self, other):  # other * self
        if isspmatrix(other):
            ocs = csc_matrix(other)
            return ocs.matmat(self)
        elif isscalar(other):
            new = self.copy()
            new.data = other * new.data
            new.typecode = new.data.typecode()
            new.ftype = _transtabl[new.typecode]
            return new
        else:
            return self.rmatvec(other)

    def __neg__(self):
        new = self.copy()
        new.data = -new.data
        return new
        
    def __sub__(self, other):
        ocs = csc_matrix(other)
        if (ocs.shape != self.shape):
            raise ValueError, "Inconsistent shapes."
        typecode = _coerce_rules[(self.typecode,other.typecode)]
        data1, data2 = _convert_data(self.data, other.data, typecode)
        func = getattr(sparsetools,_transtabl[typecode]+'cscadd')
        c,rowc,ptrc,ierr = func(data1,self.rowind,self.indptr,-data2,other.rowind,other.indptr)
        if ierr:
            raise ValueError, "Ran out of space (but shouldn't have happened)."
        M, N = self.shape
        return csc_matrix(c,(rowc,ptrc),M=M,N=N)


    def __rsub__(self, other):  # implement other - self
        ocs = csc_matrix(other)
        if (ocs.shape != self.shape):
            raise ValueError, "Inconsistent shapes."
        typecode = _coerce_rules[(self.typecode,other.typecode)]
        data1, data2 = _convert_data(self.data, other.data, typecode)
        func = getattr(sparsetools,_transtabl[typecode]+'cscadd')
        c,rowc,ptrc,ierr = func(-data1,self.rowind,self.indptr,data2,other.rowind,other.indptr)
        if ierr:
            raise ValueError, "Ran out of space (but shouldn't have happened)."
        M, N = self.shape
        return csc_matrix(c,(rowc,ptrc),M=M,N=N)
        

    # element-by-element multiplication (unless other is an
    #    integer and then matrix power)
    def __pow__(self, other):  
        if isinstance(other, type(3)):
            raise NotImplementedError
        elif isscalar(other):
            new = self.copy()
            new.data = new.data * other
            new.typecode = new.data.typecode()
            new.ftype = _transtabl[new.typecode]
            return new
        else:
            ocs = csc_matrix(other)
            if (ocs.shape != self.shape):
                raise ValueError, "Inconsistent shapes."
            typecode = _coerce_rules[(self.typecode,other.typecode)]
            nnz1, nnz2 = self.nnz, other.nnz
            data1, data2 = _convert_data(self.data[:nnz1], other.data[:nnz2], typecode)
            func = getattr(sparsetools,_transtabl[typecode]+'cscmul')
            c,rowc,ptrc,ierr = func(data1,self.rowind[:nnz1],self.indptr,data2,other.rowind[:nnz2],other.indptr)
            if ierr:
                raise ValueError, "Ran out of space (but shouldn't have happened)."
            M, N = self.shape
            return csc_matrix(c,(rowc,ptrc),M=M,N=N)

    def transp(self, copy=0):
        M,N = self.shape
        new = csr_matrix(N,M,nzmax=0,typecode=self.typecode)
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

    def rmatvec(self, x):
        if (rank(x) != 1) or (len(x) != self.shape[0]):
            raise ValueError, "Dimension mismatch"
        self._check()  # be sure everything is as it should be4
        func = getattr(sparsetools,self.ftype+'csrmux')
        y = func(self.data, self.rowind, self.indptr, x)
        return y

    def matmat(self, bmat):
        self._check()
        M,K1 = self.shape
        K2,N = bmat.shape
        if (K1 != K2):
            raise ValueError, "Shape mismatch error."
        a, rowa, ptra = self.data, self.rowind, self.indptr
        if isinstance(bmat,csr_matrix):
            bmat._check()
            typecode = _coerce_rules[(self.typecode,bmat.typecode)]
            ftype = _transtabl[typecode]            
            func = getattr(sparsetools,ftype+'cscmucsr')
            b = bmat.data
            rowb = bmat.colind
            ptrb = bmat.indptr
        elif isinstance(bmat,csc_matrix):
            bmat._check()
            typecode = _coerce_rules[(self.typecode,bmat.typecode)]
            ftype = _transtabl[typecode]                        
            func = getattr(sparsetools,ftype+'cscmucsc')
            b = bmat.data
            rowb = bmat.rowind
            ptrb = bmat.indptr
        else:
            bmat = bmat.tocsc()
            typecode = _coerce_rules[(self.typecode,bmat.typecode)]
            ftype = _transtabl[typecode]                        
            func = getattr(sparsetools,ftype+'cscmucsc')
            b = bmat.data
            rowb = bmat.rowind
            ptrb = bmat.indptr
        a, b = _convert_data(a, b, typecode)
        newshape = (M,N)
        ptrc = zeros((N+1,),'i')
        nnzc = 2*max(ptra[-1],ptrb[-1])
        c = zeros((nnzc,),typecode)
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

        return csc_matrix(c, (rowc, ptrc), M=M, N=N)
            

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
        elif isinstance(key,type(3)):
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
    
    def tocsc(self,copy=0):
        if copy:
            new = self.copy()
        else:
            new = self
        return new

    def tocoo(self,copy=0):
        func = getattr(sparsetools,self.ftype+"csctocoo")
        data,row,col = func(self.data, self.rowind,self.indptr)
        return coo_matrix(data, (row, col), M=self.shape[0], N=self.shape[1])

    def tocsr(self,copy=0):
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
        M, N = self.shape
        typecode = self.typecode
        new = csc_matrix(M, N, nzmax=0, typecode=typecode)
        new.data = self.data.copy()
        new.rowind = self.rowind.copy()
        new.indptr = self.indptr.copy()
        new._check()
        return new
    
# compressed sparse row matrix
# 
class csr_matrix(spmatrix):
    def __init__(self,s,ij=None,M=None,N=None,nzmax=100,typecode=Float,copy=0):
        spmatrix.__init__(self, 'csr')
        if isinstance(s,spmatrix):
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
                self.rowind = temp.rowind
                self.indptr = temp.indptr
                self.shape = temp.shape
        elif isinstance(s,type(3)):
            M=s
            N=ij
            self.data = zeros((nzmax,),typecode)
            self.colind = zeros((nzmax,),'i')
            self.indptr = zeros((M+1,),'i')
            self.shape = (M,N)
        elif (isinstance(s,ArrayType) or \
              isinstance(s,type([]))):
            s = asarray(s)
            if (rank(s) == 2):  # converting from a full array
                ocsc = csc_matrix(transpose(s))
                self.shape = ocsc.shape[1], ocsc.shape[0]
                self.colind = ocsc.rowind
                self.indptr = ocsc.indptr
                self.data = ocsc.data
            elif isinstance(ij, ArrayType) and (rank(ij) == 2) and (shape(ij) == (len(s),2)):
                ijnew = ij.copy()
                ijnew[:,0] = ij[:,1]
                ijnew[:,1] = ij[:,0]
                temp = coo_matrix(s,ijnew,M=M,N=N,nzmax=nzmax,
                                  typecode=typecode)
                temp = temp.tocsc()
                self.data = temp.data
                self.colind = temp.colind
                self.indptr = temp.indptr
                self.shape = temp.shape
            elif isinstance(ij, types.TupleType) and (len(ij)==2):
                self.data = asarray(s)
                self.colind = ij[0]
                self.indptr = ij[1]
                if N is None:
                    N = max(self.colind)
                if M is None:
                    M = len(self.indptr) - 1
                self.shape = (M,N)
            else:
                raise ValueError, "Unrecognized form for csr_matrix constructor."
        else:
            raise ValueError, "Unrecognized form for csr_matrix constructor."

        self._check()


    def _check(self):
        M,N = self.shape
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
        self.typecode = self.data.typecode()
        if self.typecode not in 'fdFD':
            self.typecode = 'd'
            self.data = self.data.astype('d')
        self.ftype = _transtabl[self.typecode]

        
    def __add__(self, other):
        ocs = csr_matrix(other)
        if (ocs.shape != self.shape):
            raise ValueError, "Inconsistent shapes."
        typecode = _coerce_rules[(self.typecode,other.typecode)]
        data1, data2 = _convert_data(self.data, other.data, typecode)
        func = getattr(sparsetools,_transtabl[typecode]+'cscadd')
        c,colc,ptrc,ierr = func(data1,self.colind,self.indptr,data2,other.colind,other.indptr)
        if ierr:
            raise ValueError, "Ran out of space (but shouldn't have happened)."
        M, N = self.shape
        return csr_matrix(c,(colc,ptrc),M=M,N=N)


    def __mul__(self, other):  # implement matrix multiplication and matrix-vector multiplication
        if isspmatrix(other):
            return self.matmat(other)
        elif isscalar(other):
            new = self.copy()
            new.data = new.data * other
            new.typecode = new.data.typecode()
            new.ftype = _transtabl[new.typecode]
            return new
        else:
            return self.matvec(other)

    def __rmul__(self, other):  # other * self
        if isspmatrix(other):
            ocs = csr_matrix(other)
            return occ.matmat(self)
        elif isscalar(other):
            new = self.copy()
            new.data = other * new.data
            new.typecode = new.data.typecode()
            new.ftype = _transtabl[new.typecode]
            return new
        else:
            return self.rmatvec(other)

    def __neg__(self):
        new = self.copy()
        new.data = -new.data
        return new
        
    def __sub__(self, other):
        ocs = csr_matrix(other)
        if (ocs.shape != self.shape):
            raise ValueError, "Inconsistent shapes."
        typecode = _coerce_rules[(self.typecode,other.typecode)]
        data1, data2 = _convert_data(self.data, other.data, typecode)
        func = getattr(sparsetools,_transtabl[typecode]+'cscadd')
        c,colc,ptrc,ierr = func(data1,self.colind,self.indptr,-data2,other.colind,other.indptr)
        if ierr:
            raise ValueError, "Ran out of space (but shouldn't have happened)."
        M, N = self.shape
        return csr_matrix(c,(colc,ptrc),M=M,N=N)


    def __rsub__(self, other):  # implement other - self
        ocs = csr_matrix(other)
        if (ocs.shape != self.shape):
            raise ValueError, "Inconsistent shapes."
        typecode = _coerce_rules[(self.typecode,other.typecode)]
        data1, data2 = _convert_data(self.data, other.data, typecode)
        func = getattr(sparsetools,_transtabl[typecode]+'cscadd')
        c,colc,ptrc,ierr = func(-data1,self.colind,self.indptr,data2,other.colind,other.indptr)
        if ierr:
            raise ValueError, "Ran out of space (but shouldn't have happened)."
        M, N = self.shape
        return csr_matrix(c,(colc,ptrc),M=M,N=N)
        
    # element-by-element multiplication (unless other is an
    #    integer and then matrix power)
    def __pow__(self, other):  
        if isinstance(other, type(3)):
            raise NotImplementedError
        elif isscalar(other):
            new = self.copy()
            new.data = new.data * other
            new.typecode = new.data.typecode()
            new.ftype = _transtabl[new.typecode]
            return new
        else:
            ocs = csr_matrix(other)
            if (ocs.shape != self.shape):
                raise ValueError, "Inconsistent shapes."
            typecode = _coerce_rules[(self.typecode,other.typecode)]
            data1, data2 = _convert_data(self.data, other.data, typecode)
            func = getattr(sparsetools,_transtabl[typecode]+'cscmul')
            c,colc,ptrc,ierr = func(data1,self.colind,self.indptr,data2,other.colind,other.indptr)
            if ierr:
                raise ValueError, "Ran out of space (but shouldn't have happened)."
            M, N = self.shape
            return csr_matrix(c,(colc,ptrc),M=M,N=N)

    def transp(self, copy=0):
        M,N = self.shape
        new = csc_matrix(N,M,nzmax=0,typecode=self.typecode)
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

    def rmatvec(self, x):
        if (rank(x) != 1) or (len(x) != self.shape[0]):
            raise ValueError, "Dimension mismatch"
        self._check()  # be sure everything is as it should be4
        func = getattr(sparsetools,self.ftype+'cscmux')
        y = func(self.data, self.colind, self.indptr, x, self.shape[1])
        return y

    def matmat(self, bmat):
        self._check()
        M,K1 = self.shape
        K2,N = bmat.shape
        a, rowa, ptra = self.data, self.colind, self.indptr
        if (K1 != K2):
            raise ValueError, "Shape mismatch error."
        if isinstance(bmat,csc_matrix):            
            bmat._check()
            typecode = _coerce_rules[(self.typecode,bmat.typecode)]
            ftype = _transtabl[typecode]            
            func = getattr(sparsetools,ftype+'csrmucsc')
            b = bmat.data
            colb = bmat.rowind
            ptrb = bmat.indptr
            out = 'csc'
            firstarg = ()
        elif isinstance(bmat,csr_matrix):
            bmat._check()
            typecode = _coerce_rules[(self.typecode,bmat.typecode)]
            ftype = _transtabl[typecode]            
            func = getattr(sparsetools,ftype+'cscmucsc')
            b, colb, ptrb = a, rowa, ptra
            a, rowa, ptra = bmat.data, bmat.colind, bmat.indptr
            out = 'csr'
            firstarg = (N,)
        else:
            bmat = bmat.tocsc()
            typecode = _coerce_rules[(self.typecode,bmat.typecode)]
            ftype = _transtabl[typecode]            
            func = getattr(sparsetools,ftype+'csrmucsc')
            b = bmat.data
            colb = bmat.colind
            ptrb = bmat.indptr
            out = 'csc'
            firstarg = ()
        a, b = _convert_data(a, b, typecode)            
        newshape = (M,N)
        if out == 'csr':
            ptrc = zeros((M+1,),'i')
        else:
            ptrc = zeros((N+1,),'i')
        nnzc = 2*max(ptra[-1],ptrb[-1])
        c = zeros((nnzc,),typecode)
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

        outinit = eval('%s_matrix' % out)
        return outinit(c, (rowc, ptrc), M=M, N=N)


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
        elif isinstance(key,type(3)):
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

    def tocsr(self,copy=0):
        if copy:
            new = self.copy()
        else:
            new = self
        return new

    def tocoo(self,copy=0):
        func = getattr(sparsetools,self.ftype+"csctocoo")
        data,col,row = func(self.data, self.colind,self.indptr)
        return coo_matrix(data, (row, col), M=self.shape[0], N=self.shape[1])

    def tocsc(self,copy=0):
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
        M, N = self.shape
        typecode = self.typecode
        new = csr_matrix(M, N, nzmax=0, typecode=typecode)
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
    def __init__(self,A=None):
        dict.__init__(self)
        spmatrix.__init__(self,'dok')
        self.shape = (0,0)
        self.nnz = 0
        if A is not None:
            A = asarray(A)
            N,M = A.shape
            for n in range(N):
                for m in range(M):
                    if A[n,m] != 0:
                        self[n,m] = A[n,m]

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
        return self.get(key,0)

    def __setitem__(self, key, value):
        if (value == 0):
            return
        if not isinstance(key, tuple) or len(key) != 2:
            raise KeyError, "Key must be a 2-tuple"
        if (self.get(key) == None): new = 1
        else: new = 0
        dict.__setitem__(self, key, value)
        newrows = max(self.shape[0], int(key[0])+1)
        newcols = max(self.shape[1], int(key[1])+1)
        self.shape = (newrows, newcols)
        if new:
            self.nnz += 1
    
    def __add__(self, other):
        res = dok_matrix()
        res.update(self)
        res.shape = self.shape
        for key in other.keys():
            try:
                res[key] += other[key]
            except KeyError:
                res[key] = other[key]
        return res

    def __sub__(self, other):
        res = dok_matrix()
        res.update(self)
        res.shape = self.shape
        for key in other.keys():
            try:
                res[key] -= other[key]
            except KeyError:
                res[key] = -other[key]
        return res
    
    def __neg__(self):
        res = dok_matrix()
        for key in self.keys():
            res[key] = -self[key]
        return res

    def __mul__(self, other):
        if isinstance(other, dok_matrix):
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

    def transp(self):
        # Transpose (return the transposed)
        new = dok_matrix()
        for key in self.keys():
            new[key[1],key[0]] = self[key]
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
        keys = self.keys()
        res = [0]*self.shape[0]
        for key in keys:
            res[int(key[0])] += self[key] * other[int(key[1]),...]
        return array(res)        

    def setdiag(self, values, k=0):
        N = len(values)
        for n in range(N):
            self[n,n+k] = values[n]
        return

    def tocsr(self):
        # Return Compressed Sparse Row format arrays for this matrix
        keys = self.keys()
        keys.sort()
        nnz = len(keys)
        data = [0]*nnz
        colind = [0]*nnz
        row_ptr = [0]*(self.shape[0]+1)
        current_row = -1
        k = 0
        for key in keys:
            ikey0 = int(key[0])
            ikey1 = int(key[1])
            if ikey0 != current_row:
                current_row = ikey0
                row_ptr[ikey0] = k
            data[k] = self[key]
            colind[k] = ikey1
            k += 1
        row_ptr[-1] = nnz
        data = array(data)
        colind = array(colind)
        row_ptr = array(row_ptr)
        return csr_matrix(data,(colind, row_ptr))

    def tocsc(self):
        # Return Compressed Sparse Column format arrays for this matrix
        keys = self.keys()
        keys.sort(csc_cmp)
        nnz = len(keys)
        data = [None]*nnz
        colind = [None]*nnz
        col_ptr = [None]*(self.shape[1]+1)
        current_col = -1
        k = 0
        for key in keys:
            ikey0 = int(key[0])
            ikey1 = int(key[1])
            if ikey1 != current_col:
                current_col = ikey1
                col_ptr[ikey1] = k
            data[k] = self[key]
            colind[k] = ikey0
            k += 1
        col_ptr[-1] = nnz
        data = array(data)
        colind = array(colind)
        col_ptr = array(col_ptr)
        return csc_matrix(data, (colind, col_ptr))

    def todense(self,typecode=None):
        if typecode is None:
            typecode = 'd'
        new = zeros(self.shape,typecode)
        for key in self.keys():
            ikey0 = int(key[0])
            ikey1 = int(key[1])
            new[ikey0,ikey1] = self[key]
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
    def __init__(self, obj, ij, M=None, N=None, nzmax=None, typecode=None):
        spmatrix.__init__(self, 'coo')
        if type(ij) is type(()) and len(ij)==2:
            if M is None:
                M = amax(ij[0])
            if N is None:
                N = amax(ij[1])
            self.row = asarray(ij[0],'i')
            self.col = asarray(ij[1],'i')
        else:
            aij = asarray(ij,'i')
            if M is None:
                M = amax(aij[:,0])
            if N is None:
                N = amax(aij[:,1])
            self.row = aij[:,0]
            self.col = aij[:,1]
        aobj = asarray(obj)
        self.shape = (M,N)
        if nzmax is None:
            nzmax = len(aobj)
        self.nzmax = nzmax
        self.data = aobj
        self.typecode = aobj.typecode()
        self._check()

    def _check(self):
        nnz = len(self.data)
        if (nnz != len(self.row)) or (nnz != len(self.col)):
            raise ValueError, "Row, column, and data array must all be "\
                  "the same length."
        if (self.nzmax < nnz):
            raise ValueError, "nzmax must be >= nnz"
        self.nnz = nnz
        self.ftype = _transtabl[self.typecode]

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
        self.data = array(data,self.typecode)
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
        return csc_matrix(a, (rowa, ptra), M=self.shape[0], N=self.shape[1])

    def tocsr(self):
        func = getattr(sparsetools,self.ftype+"cootocsc")
        data,row,col = self._normalize(rowfirst=1)
        a, cola, ptra, ierr = func(self.shape[0], data, col, row)
        if ierr:
            raise RuntimeError, "Error in conversion."
        return csr_matrix(a, (cola, ptra), M=self.shape[0],N=self.shape[1])

    def tocoo(self,copy=0):
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

def _spdiags_tosub(diag_num,a,b):
    part1 = where(less(diag_num,a),abs(diag_num-a),0)
    part2 = where(greater(diag_num,b),abs(diag_num-b),0)
    return part1+part2
                      
def spdiags(diags,offsets,m,n):
    """Return a sparse matrix given it's diagonals.
    
    B = spdiags(diags, offsets, M, N)

    Inputs:
        diags  --  rows contain diagonal values
        offsets -- diagonals to set (0 is main)
        M, N    -- sparse matrix returned is M X N
    """
    diags = array(transpose(diags),copy=1)
    if diags.typecode() not in 'fdFD':
        diags = diags.astype('d')
    offsets = array(offsets,copy=0)
    mtype = diags.typecode()
    assert(len(offsets) == diags.shape[1])
    # set correct diagonal to csr conversion routine for this type
    diagfunc = eval('sparsetools.'+_transtabl[mtype]+'diatocsr')
    a, rowa, ptra, ierr = diagfunc(m,n,diags,offsets)
    if ierr:
        raise ValueError, "Ran out of memory (shouldn't have happened)"
    return csc_matrix(a,(rowa,ptra),M=m,N=n)

def solve(A,b,permc_spec=2):
    if not hasattr(A, 'tocsr') and not hasattr(A, 'tocsc'):
        raise ValueError, "Sparse matrix must be able to return CSC format--"\
              "A.tocsc()--or CSR format--A.tocsr()"
    if not hasattr(A,'shape'):
        raise ValueError, "Sparse matrix must be able to return shape (rows,cols) = A.shape"
    M,N = A.shape
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
    M,N = A.shape
    if (M != N):
        raise ValueError, "Can only factor square matrices."
    csc = A.tocsc()
    gstrf = eval('_superlu.' + csc.ftype + 'gstrf')
    return gstrf(N,csc.nnz,csc.data,csc.rowind,csc.indptr,permc_spec,
                 diag_pivot_thresh, drop_tol, relax, panel_size)
        

if __name__ == "__main__":
    a = spmatrix(arange(1,9),[0,1,1,2,2,3,3,4],[0,1,3,0,2,3,4,4])
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
