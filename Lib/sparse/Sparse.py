from Numeric import *
from scipy_base.fastumath import *
from scipy_base import isscalar, rank, shape, resize, ArrayType
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
ALLOCSIZE = 1000
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
    def __init__(self, format, maxprint=MAXPRINT):
        self.format = format
        self.maxprint = maxprint
        
    def getmaxprint(self):
        try:
            maxprint = self.maxprint
        except AttributeError:
            maxprint = MAXPRINT
        return maxprint
        
    def getnumtype(self):
        try:
            numtype = self.numtype
        except AttributeError:
            numtype = None
        return numtype

    def getnnz(self):
        try:
            nnz = self.nnz
        except AttributeError:
            nnz = 0
        return nnz

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
        numtype = self.getnumtype()
        nnz = self.getnnz()
        format = self.getformat()
        return "<%dx%d sparse matrix of type '%s' with %d non-zero "\
               "elements in %s format>" % \
               (self.shape + (numtype, nnz, _formats[format][1]))

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
    # thus, a new sparse matrix format just needs to define
    #  a tocsc method and then take a csc_matrix in its constructor
    #  to get functionality (though not optimality)

    def __add__(self, other):
        format = self.getformat()
        csc = self.tocsc()
        res = csc + other
        return eval('%s_matrix'%format)(res)

    def __sub__(self, other):
        format = self.getformat()
        csc = self.tocsc()
        res = csc - other
        return eval('%s_matrix'%format)(res)

    def __rsub__(self, other): # other - self
        format = self.getformat()
        csc = self.tocsc()
        res = csc.__rsub__(other)
        return eval('%s_matrix'%format)(res)


    def __mul__(self, other):
        format = self.getformat()
        csc = self.tocsc()
        res = csc * other
        return eval('%s_matrix'%format)(res)

    def __rmul__(self, other):
        format = self.getformat()
        csc = self.tocsc()
        res = csc.__rmul__(other)
        return eval('%s_matrix'%format)(res)
        
    def __neg__(self):
        format = self.getformat()
        csc = self.tocsc()
        res = -csc
        return eval('%s_matrix'%format)(res)

    def matvec(self, vec):
        format = self.getformat()
        csc = self.tocsc()
        res = csc.matvec(vec)
        return res

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
            temp = s.tocsc(copy=copy)
            self.data = temp.data
            self.rowind = temp.rowind
            self.indptr = temp.indptr
            self.shape = temp.shape
            self.numtype = temp.numtype
            self.nnz = temp.nnz
            self.nzmax = nzmax
        elif isinstance(s,type(3)):
            M=s
            N=ij
            self.data = zeros((nzmax,),typecode)
            self.rowind = zeros((nzmax,),'i')
            self.indptr = zeros((N+1,),'i')
            self.shape = (M,N)
        elif (isinstance(s,ArrayType) or \
              isinstance(s,type([]))):
            if (rank(s) == 2):  # converting from a full array
                M, N = s.shape
                s = asarray(s)
                if s.typecode() not in 'fdFD':
                    s = s*1.0
                numtype = s.typecode()
                func = getattr(sparsetools,_transtabl[numtype]+'fulltocsc')
                ierr = 1
                irow = jcol = 0
                nnz = sum(ravel(s != 0.0))
                a = zeros((nnz,),numtype)
                rowa = zeros((nnz,),'i')
                ptra = zeros((N+1,),'i')
                while 1:
                    a, ptra, rowa, irow, jcol, ierr = func(s, a, ptra, rowa, irow, jcol)
                    if (ierr == 0): break
                    nnz = nnz + ALLOCSIZE
                    a = resize1d(a, nnz)
                    rowa = resize1d(rowa, nnz)

                self.data = a
                self.rowind = rowa
                self.indptr = ptra
                self.shape = (M,N)
                    
            elif (rank(ij) == 2) and (shape(ij) == (len(s),2)):
                temp = coo_matrix(s,ij,M=M,N=N,nzmax=nzmax,typecode=typecode)
                temp = temp.tocsc()
                self.data = temp.data
                self.rowind = temp.rowind
                self.indptr = temp.indptr
                self.shape = temp.shape
            else:
                self.data = asarray(s)
                self.rowind = ij[0]
                self.indptr = ij[1]
                if M is None:
                    M = max(self.rowind)
                if N is None:
                    N = len(self.indptr) - 1
                self.shape = (M,N)
        else:
            raise ValueError, "Unrecognized form for csc_matrix constructor."

        self._check()


    def _check(self):
        M,N = self.shape
        if (rank(self.data) != 1) or (rank(self.rowind) != 1) or \
           (rank(self.rowind) != 1):
            raise ValueError, "Data, row, and indptr arrays should be rank 1."
        if (len(self.data) != len(self.rowind)):
            raise ValueError, "Data and row list should have same length"
        if (len(self.indptr) != N+1):
            raise ValueError, "Index pointer should be of of size N+1"
        if (max(self.rowind) >= M):
            raise ValueError, "Row-values must be < M."
        if (self.indptr[-1] > len(self.rowind)):
            raise ValueError, \
                  "Last value of index list should be less than "\
                  "the size of data list"
        self.nnz = self.indptr[-1]
        self.nzmax = len(self.rowind)
        self.numtype = self.data.typecode()
        self.ftype = _transtabl[self.numtype]
        if self.numtype not in 'fdFD':
            raise ValueError, "Only floating point sparse matrix types allowed"
        

    def __add__(self, other):
        ocs = csc_matrix(other)
        if (ocs.shape != self.shape):
            raise ValueError, "Inconsistent shapes."
        numtype = _coerce_rules[(self.numtype,other.numtype)]
        data1, data2 = _convert_data(self.data, other.data, numtype)
        func = getattr(sparsetools,_transtabl[numtype]+'cscadd')
        c,rowc,ptrc,ierr = func(data1,self.rowind,self.indptr,data2,other.rowind,other.indptr)
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
            new.numtype = new.data.typecode()
            new.ftype = _transtabl[new.numtype]
            return new
        else:
            return self.matvec(other)

    def __neg__(self):
        new = self.copy()
        new.data = -new.data
        return new
        
    def __sub__(self, other):
        ocs = csc_matrix(other)
        if (ocs.shape != self.shape):
            raise ValueError, "Inconsistent shapes."
        numtype = _coerce_rules[(self.numtype,other.numtype)]
        data1, data2 = _convert_data(self.data, other.data, numtype)
        func = getattr(sparsetools,_transtabl[numtype]+'cscadd')
        c,rowc,ptrc,ierr = func(data1,self.rowind,self.indptr,-data2,other.rowind,other.indptr)
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
            new.numtype = new.data.typecode()
            new.ftype = _transtabl[new.numtype]
            return new
        else:
            ocs = csc_matrix(other)
            if (ocs.shape != self.shape):
                raise ValueError, "Inconsistent shapes."
            numtype = _coerce_rules[(self.numtype,other.numtype)]
            data1, data2 = _convert_data(self.data, other.data, numtype)
            func = getattr(sparsetools,_transtabl[numtype]+'cscmul')
            c,rowc,ptrc,ierr = func(data1,self.rowind,self.indptr,data2,other.rowind,other.indptr)
            if ierr:
                raise ValueError, "Ran out of space (but shouldn't have happened)."
            M, N = self.shape
            return csc_matrix(c,(rowc,ptrc),M=M,N=N)

    def matvec(self, x):
        if (rank(x) != 1) or (len(x) != self.shape[1]):
            raise ValueError, "Dimnsion mismatch"
        self._check()  # be sure everything is as it should be
        func = getattr(sparsetools,self.ftype+'cscmux')
        y = func(self.data, self.rowind, self.indptr, x, self.shape[0])
        return y

    def matmat(self, bmat):
        self._check()
        M,K1 = self.shape
        K2,N = bmat.shape
        if (K1 != K2):
            raise ValueError, "Shape mismatch error."
        if isinstance(x,csr_matrix):
            bmat._check()
            func = getattr(sparsetools,bmat.ftype+'cscmucsr')
            b = bmat.data
            colb = bmat.colind
            ptrb = bmat.indptr
        elif isintance(x,csc_matrix):
            bmat._check()
            func = getattr(sparsetools,bmat.ftype+'cscmucsc')
            b = bmat.data
            colb = bmat.rowind
            ptrb = bmat.indptr
        else:
            bmat = bmat.tocsc()
            func = getattr(sparsetools,bmat.ftype+'cscmucsc')
            b = bmat.data
            colb = bmat.rowind
            ptrb = bmat.indptr
        newshape = (M,N)
        ptrc = zeros((N+1,),'i')
        nnzc = 2*max(len(self.data)+len(b))
        c = zeros((nnzc,),typecode)
        rowc = zeros((nnzc,),'i')
        irow = kcol = 0
        while 1:
            c, rowc, ptrc, irow, kcol, ierr = func(M,a,rowa,ptra,b,colb,ptrb,c,rowc,ptrc,irow,kcol)
            if (ierr==0): break
            # otherwise we were too small and must resize
            percent_to_go = 1- (1.0*kcol) / N
            newnnzc = int(ceil((1+percent_to_go)*nnzc))
            c = resize1d(c,newnnzc)
            newrowc = resize1d(rowc,newnnzc)
            nnzc = newnnzc
            
    def prune():
        nnz = self.indptr[-1]
        self.nnz = nnz
        self.data = self.data[:nnz]
        self.rowind = self.rowind[:nnz]
        self.nzmax = nnz
        self._check()

    def copy():
        M, N = self.shape
        numtype = self.numtype
        new = csc_matrix(M, N, nzmax=0, typecode=numtype)
        new.data = self.data.copy()
        new.rowind = self.rowind.copy()
        new.indptr = self.indptr.copy()
        new._check()
        return new
    
# compressed sparse row matrix
# 
class csr_matrix(spmatrix):
    def __init__(self):
        pass


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
        self.numtype = None
        if A is not None:
            A = asarray(A)
            N,M = A.shape
            for n in range(N):
                for m in range(M):
                    if A[n,m] != 0:
                        self[n,m] = A[n,m]

    def __getitem__(self, key):
        return self.get(key,0)

    def __setitem__(self, key, value):
        if (value == 0):
            return
        if not isinstance(key, tuple) or len(key) != 2:
            raise KeyError, "Key must be a 2-tuple"
        dict.__setitem__(self, key, value)
        newrows = max(self.shape[0], int(key[0])+1)
        newcols = max(self.shape[1], int(key[1])+1)
        self.shape = (newrows, newcols)
    
    def __add__(self, other):
        res = dictmatrix()
        res.update(self)
        res.shape = self.shape
        for key in other.keys():
            try:
                res[key] += other[key]
            except KeyError:
                res[key] = other[key]
        return res

    def __sub__(self, other):
        res = dictmatrix()
        res.update(self)
        res.shape = self.shape
        for key in other.keys():
            try:
                res[key] -= other[key]
            except KeyError:
                res[key] = -other[key]
        return res
    
    def __neg__(self):
        res = dictmatrix()
        for key in self.keys():
            res[key] = -self[key]
        return res

    def __mul__(self, other):
        if isinstance(other, dictmatrix):
            return self.matmat(other)
        other = asarray(other)
        if rank(other) > 0:
            return self.matvec(other)
        res = dictmatrix()
        for key in self.keys():
            res[key] = other * self[key]
        return res

    def __len__(self):
        return len(self.keys())

    def transp(self):
        # Transpose (return the transposed)
        new = dictmatrix()
        for key in self.keys():
            new[key[1],key[0]] = self[key]
        return new

    def matmat(self, other):
        res = dictmatrix()
        spself = spmatrix(self)
        spother = spmatrix(other)
        spres = spself * spother
        return spres.todict()                            

    def take(self, cols_or_rows, columns=1):
        # Extract columns or rows as indictated from matrix
        # assume cols_or_rows is sorted
        res = dictmatrix()
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
        base = dictmatrix()
        ext = dictmatrix()
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

    def getCSR(self):
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
        ptype = data.typecode()
        if ptype not in ['d','D','f','F']:
            data = data.astype('d')
            ptype = 'd'
        return _transtabl[ptype], nnz, data, colind, row_ptr

    def getCSC(self):
        # Return Compressed Sparse Column format arrays for this matrix
        keys = self.keys()
        keys.sort(csc_cmp)
        nnz = len(keys)
        data = [None]*nnz
        rowind = [None]*nnz
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
            rowind[k] = ikey0
            k += 1
        col_ptr[-1] = nnz
        data = array(data)
        rowind = array(rowind)
        col_ptr = array(col_ptr)
        ptype = data.typecode()
        if ptype not in ['d','D','f','F']:
            data = data.astype('d')
            ptype = 'd'
        return _transtabl[ptype], nnz, data, rowind, col_ptr

    def dense(self,typecode=None):
        if typecode is None:
            typecode = self.type
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
class coo_matrix(spmatrix):
    pass

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
    diagfunc = eval('_sparsekit.'+_transtabl[mtype]+'diacsr')
    # construct empty sparse matrix and pass it's main parameters to
    #  the diagonal to csr conversion routine.
    nzmax = diags.shape[0]*diags.shape[1]
    s = spmatrix(m,n,nzmax,typecode=mtype)
    diagfunc(array(m), array(n), array(0),  diags,
             offsets, s.data, s.index[0], s.index[1],
             array(diags.shape[1]),array(diags.shape[0]))

    # compute how-many elements were actually filled
    s.lastel = min([m,n])*len(offsets) - 1 - \
               sum(_spdiags_tosub(offsets, a=min([n-m,0]), b=max([n-m,0])))
    return s

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
               A.ftype, A.nnz, A.data, A.rowind, A.indptr
        csc = 1
    else:
        mat = A.tocsr()
        ftype, lastel, data, index0, index1 = \
               A.ftype, A.nnz, A.data, A.colind, A.indptr
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
    return gstrf(N,csc.nnz,csc.data,csc.rowind,csc.colptr,permc_spec,
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
