from Numeric import *
from scipy_base.fastumath import *
from scipy_base import isscalar
import types
import _sparsekit
import _sparseutil
import _superlu

MAXPRINT=50
# The formats that SPARSEKIT's convert programs understand.

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
    def __init__(self, format):
        self.format = format
        
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
            format = self.storage
        except AttributeError:
            format = 'UND'
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
               (self.shape + (numtype, nzmax, _formats[format][1]))

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
        csc = self.toscs()
        return eval('csc.as%s' % format)()
 
    # default operations use the CSC format as a base
    # thus, a new sparse matrix format just needs to define
    #  an toscs method and then get CSCmatrix to define an
    #  tofmt type to get functionality.

    def __add__(self, other):
        format = self.getformat()
        csc = self.toscs()
        res = csc + other
        return eval('res.to%s' % format)

    def __sub__(self, other):
        format = self.getformat()
        csc = self.toscs()
        res = csc - other
        return eval('res.to%s' % format)

    def __rsub__(self, other): # other - self
        format = self.getformat()
        csc = self.toscs()
        res = csc.__rsub__(other)
        return eval('res.to%s' % format)

    def __mul__(self, other):
        format = self.getformat()
        csc = self.toscs()
        res = csc * other
        return eval('res.to%s' % format)

    def __rmul__(self, other):
        format = self.getformat()
        csc = self.toscs()
        res = csc.__rmul__(other)
        return eval('res.to%s' % format)
        
    def __neg__(self):
        format = self.getformat()
        csc = self.toscs()
        res = -csc
        return eval('res.to%s' % format)
        
# compressed sparse column matrix
#  This can be instantiated in many ways
#    - with another sparse matrix (sugar for .tocsc())
#    - with M,N,nzmax,typecode  to construct a container
#    - with data, ij, {M,N,nzmax}
#           a[ij[k,0],ij[k,1]] = data[k]
# 

class csc_matrix(spmatrix):
    def __init__(self,s,ij=None,M=None,N=None,nzmax=100,typecode=Float):       
        spmatrix.__init__(self, 'CSC')
        if isinstance(s,spmatrix):
            temp = s.tocsc()
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
            self.rowind = zeros((nzmax,))
            self.indptr = zeros((N+1,))
            self.shape = (M,N)
            self.numtype = typecode
            self.nnz = 0
            self.nzmax = nzmax
        elif isinstance(s,scipy_base.ArrayType) or \
             isinstance(s,type([])):
            temp = coo_matrix(s,ij,M=M,N=N,nzmax=nzmax,typecode=typecode)
            temp = temp.tocsc()
            self.data = temp.data
            self.rowind = temp.rowind
            self.indptr = temp.indptr
            self.shape = temp.shape
            self.numtype = temp.numtype
            self.nnz = temp.nnz
            self.nzmax = nzmax

    def __add__(self, other):
        if isspmatrix(other):
            ocs = other.tocsc()
            numtype = _coerce_rules[(self.numtype,other.numtype)]
            data1, data2 = _coerce_data(self.data1, other.data2, numtype)
            func = eval('sparsetools.'+_transtabl[numtype]+'cscadd')
    
# compressed sparse row matrix
# 
class csr_matrix(spmatrix):
    def __init__(self):
        pass
    

# dictionary of keys based matrix
class dok_matrix(spmatrix):
    pass

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
class sss_matrix(spmatrix):
    pass

# diagonal (banded) matrix
class dia_matrix(spmatrix):
    pass

# ellpack-itpack generalized diagonal
class egd_matrix(spmatrix):
    pass

# block sparse row
# modified compressed sparse row
# block sparse column
# modified compressed sparse column
# symmetric skyline
# nonsymmetric skyline
# jagged diagonal
# unsymmetric sparse skyline
# variable block row



# So far only CSR format is supported internally.




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
        

class dictmatrix(dict):
    def __init__(self,A=None):
        dict.__init__(self)
        self.shape = (0,0)
        self.storage = 'dict'
        self.type = None
        self.maxprint = MAXPRINT
        if A is not None:
            A = asarray(A)
            N,M = A.shape
            for n in range(N):
                for m in range(M):
                    if A[n,m] != 0:
                        self[n,m] = A[n,m]

    def __repr__(self):
        return "<%dx%d dictmatrix of type '%s' with %d non-zero elements>" % (self.shape + (self.type, len(self.keys())))

    def __str__(self):
        val = ''
        nnz = len(self.keys())
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

# A sparse matrix class.  A sparse matrix can be initialized as:

# a = spmatrix(M,N,nzmax,typecode=Float)
#   Create an MxN matrix with room for nzmax non-zero elements of
#     type typecode
# a = spmatrix(data,row,col{,M,N,nzmax})
#   Create a sparse matrix with non-zero elements data
#     using a[row[i],col[i]] = data[i]
               
class spmatrix:
    def __init__(self,s,i=None,j=None,M=None,N=None,nzmax=None,
                 typecode=Float):
        if isinstance(s, dictmatrix):
            ftype, nnz, data, index0, index1 = s.getCSR()
            self.ftype = ftype
            self.ptype = _itranstabl[ftype]
            self.lastel = nnz-1
            self.data = data
            self.index = [index0+1, index1+1]
            M, N = s.shape
            nzmax = nnz
        elif type(s) in [types.ListType, ArrayType]:
            s = array(s,copy=0,typecode=typecode)
            if s.typecode() not in 'fdFD':  # only support these 4 types.
                s = s.astype('d')
            sz = len(s)
            i = array(i,typecode='l',copy=0)
            j = array(j,typecode='l',copy=0)
            if nzmax is None:
                nzmax = sz
                if M is None:
                    M = max(i)+1
                if N is None:
                    N = max(j)+1
            self.ptype = s.typecode()
            self.data = zeros((nzmax,),s.typecode())
            self.ftype = _transtabl[self.ptype]
            self.index = [zeros((nzmax,)),zeros((M+1,))]
            convfunc = eval('_sparsekit.'+self.ftype+'coocsr')
            convfunc(array(M),array(nzmax),s,i+1,j+1,self.data,self.index[0],self.index[1])
            self.lastel = len(s)-1
        elif type(s) is types.IntType:
            M = int(s)
            N = int(i)
            if j is None:
                j = 0
            nzmax = int(j)
            self.ptype = typecode
            self.ftype = _transtabl[self.ptype]
            self.data = zeros((nzmax,),typecode)
            self.index = [zeros((nzmax,)),zeros((M+1,))]
            self.lastel = 0
        elif isspmatrix(s) and s.storage=='CSR':  # make a copy
            for attr in dir(s):
                if attr not in ['data','index']:
                    setattr(self,attr,getattr(s,attr))
            self.data = array(s.data,copy=1)
            self.index = [array(s.index[0],copy=1),array(s.index[1],copy=1)]
            return
        else:
            raise TypeError, "Unsupported type %s" % type(s)

        self.storage = 'CSR'
        self.shape = (M,N)
        self.nzmax = nzmax
        self.maxprint = MAXPRINT

    def rowcol(self,key):
        if key > self.lastel:
            raise ValueError, "There are only %d nonzero entries." % (self.lastel+1)
        row = searchsorted(self.index[1]-1,key+1)-1
        col = self.index[0][key]-1
        return (row,col)

    def listprint(self,start,stop):
        val = ''
        for ind in xrange(start,stop):
            val = val + '  %s\t%s\n' % (self.rowcol(ind),self.data[ind])
        return val

    def __repr__(self):
        return "<%dx%d spmatrix of type '%s' with %d elements in %s>" % (self.shape + (self.ptype, self.nzmax, _formats[self.storage][1]))

    def __str__(self):
        val = ''
        if self.nzmax > self.maxprint:
            val = val + self.listprint(0,self.maxprint/2)
            val = val + "  :\t:\n"
            val = val + self.listprint(self.lastel-self.maxprint/2,self.lastel+1)
        else:
            val = val + self.listprint(0,self.lastel+1)
        return val[:-1]

    def __cmp__(self,other):
        raise TypeError, "Comparison of sparse matrices is not implemented."

    def __nonzero__(self):  # Simple -- other ideas?
        return self.lastel > 0

    def __len__(self):
        return self.lastel + 1

    def __getitem__(self,key):  
        if key is None:
            return 0.0
        if type(key) == types.IntType:
            return (self.data[key], self.rowcol(key))
        elif type(key) == types.TupleType:
            if len(key) == 2 and (type(key[0]),type(key[1])) == (types.IntType,)*2:
                getelm = eval('_sparsekit.'+self.ftype+'getelm')
                add = array(0)
                val = getelm(array(key[0]+1),array(key[1]+1),self.data,self.index[0],self.index[1],add,array(0))
                if add[0] > 0:
                    return (self.data[add-1],add-1)
                else:
                    return (0.0, None)
            elif len(key) == 2 and (type(key[0]),type(key[1])) == (types.SliceType,)*2:
                if (key[0].step is not None and key[0].step != 1) or \
                   (key[1].step is not None and key[1].step != 1):
                    print "Ignoring step value in slice."
                assert self.storage=='CSR'
                submat = eval('_sparsekit.'+self.ftype+'submat')
                nr = array(0)
                nc = array(0)
                ao = array(self.data,copy=1)
                jao = array(self.index[0],copy=1)
                iao = zeros((self.lastel+1,))
                submat(array(self.shape[0]),array(1),array(key[0].start+1),array(key[0].stop),array(key[1].start+1),array(key[1].stop),self.data,self.index[0],self.index[1],nr,nc,ao,jao,iao)
                nels = max(iao)-1
                # Eliminate "extra memory"
                ao = array(ao[:nels],copy=1)
                jao = array(jao[:nels],copy=1)
                iao = array(iao[:nr[0]],copy=1)
                b = spmatrix(nr[0],nc[0],nels)
                b.lastel = nels-1
                b.data = ao
                b.index = [jao,iao]
                return b                                   
        raise TypeError, "Cannot access sparse matrix that way."


    def astype(self,newtype):
        if newtype == self.ptype:
            return self
        else:
            b = spmatrix(self)
            b.data = b.data.astype(newtype)
            b.ptype = newtype
            b.ftype = _transtabl[newtype]
            return b

    def __add__(self,other):
        if not isspmatrix(other):
            raise TypeError, "Both matrices must be sparse."
        spadd = eval('_sparsekit.'+self.ftype+'aplb')
        assert self.shape == other.shape
        assert self.storage == 'CSR'
        if other.ftype != self.ftype:
            other = other.astype(self.ptype)
        new = spmatrix(self.shape[0],self.shape[1],min((self.nzmax + other.nzmax,product(self.shape))),typecode=self.ptype)
        ierr = array(0)
        iw = zeros((self.shape[1],))
        spadd(1,self.data,self.index[0],self.index[1],other.data,other.index[0],other.index[1],new.data,new.index[0],new.index[1],array(new.nzmax),iw,ierr,self.shape[0],self.shape[1])
        nels = max(new.index[1])-1
        new.data = array(new.data[:nels],copy=1)
        new.index[0] = array(new.index[0][:nels],copy=1)
        new.lastel = nels - 1
        new.nzmax = nels
        return new

    def __neg__(self):
        new = spmatrix(self.shape[0],self.shape[1],self.nzmax)
        new.data = -self.data
        new.index = self.index
        new.ptype = self.ptype
        new.ftype = self.ftype
        new.lastel = self.lastel
        return new

    def __sub__(self,other):
        if not isspmatrix(other):
            raise TypeError, "Right operand must also be sparse."
        return self + (-other)

    def __mul__(self,other):
        if isspmatrix(other):
            assert other.shape[0] == self.shape[1]
            assert self.storage == 'CSR'
            new_nz = self.nzmax + other.nzmax
            new = spmatrix(self.shape[0],other.shape[1],new_nz,typecode=self.ptype)
            mult = eval('_sparsekit.'+self.ftype+'amub')
            iw = zeros((self.shape[1],))
            ierr = array(0)
            while 1:  # mult returns error if array wasn't big enough
                mult(array(1),
                     self.data, self.index[0], self.index[1], other.data,
                     other.index[0], other.index[1], new.data, new.index[0],
                     new.index[1], array(new.nzmax), iw, ierr,
                     array(self.shape[0]),array(other.shape[1]))
                if (ierr[0] == 0 or new.nzmax > 5*self.nzmax):
                    break
                # make output array bigger for the next try
                new.expand(int(self.shape[0]/float(ierr[0])*new.nzmax + 1))
            if (ierr[0] != 0):
                raise ValueError, "Could not find a good size for sparse output: ierr = %d" % ierr[0]
            new.cleanup()
            return new
        
        elif type(other) in [ArrayType, types.ListType]:
            assert self.storage == 'CSR'
            other = array(other,copy=0).astype(self.ptype)
            assert len(other.shape)==1 and len(other) == self.shape[1]
            matvec = eval('_sparsekit.'+self.ftype+'amux')
            y = zeros((self.shape[0]),self.ptype)
            matvec(array(self.shape[0]),other,y,self.data,self.index[0],self.index[1])
            return y
        
        elif type(other) in [types.IntType, types.FloatType, types.ComplexType]:
            new = spmatrix(self)           # make a copy
            new.data = other*new.data
            new.ptype = new.data.typecode()
            new.ftype = _transtabl[new.ptype]
            return new

    def __rmul__(self,other):
        return self*other

    def cleanup(self):  # eliminate unused entries from the end
        assert self.storage=='CSR'
        self.nzmax = searchsorted(equal(self.index[0],0),1)
        self.lastel = self.nzmax-1
        self.data = self.data[:self.nzmax]
        self.index[0] = self.index[0][:self.nzmax]

    def expand(self, new_nzmax):  # create more space
        assert self.storage=='CSR'
        tmp = self.data
        tmpb = self.index[0]
        self.data = zeros((new_nzmax,),self.ptype)
        self.index[0] = zeros((new_nzmax,),self.ptype)
        self.data[:self.nzmax] = tmp
        self.index[0][:self.nzmax] = tmpb
        self.nzmax = new_nzmax

    def transp(self,inplace=0):
        assert self.storage=='CSR'
        M,N = self.shape
        if inplace == 0:
            new = spmatrix(self) # make a copy
        else:
            new = self   # make a reference
        transp = eval('_sparsekit.'+self.ftype+'transp')
        iwk = zeros((len(self.index[0]),))
        ierr = array(0)
        M,N = array(M),array(N)
        transp(M,N,new.data,new.index[0],new.index[1],iwk,ierr)
        if ierr[0] != 0:
            raise ValueError, "Error during transpose: %d"  % ierr[0]
        new.shape = (N[0],M[0])
        return new

    def conj(self,inplace=0):
        if inplace == 0:
            new = spmatrix(self)
        else:
            new = self
        new.data = conjugate(self.data)
        return new

    def getCSR(self):
        return self.ftype, self.lastel+1, self.data, self.index[0]-1, self.index[1]-1

    def getCSC(self):
        B = self.transp()
        return B.ftype, B.lastel+1, B.data, B.index[0]-1, B.index[1]-1

    def todict(self):
        res = dictmatrix()
        for k in range(self.nzmax):
            row,col = self.rowcol(k)
            res[row,col] = self.data[k]
        return res
            
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
    if not hasattr(A, 'getCSR') and not hasattr(A, 'getCSC'):
        raise ValueError, "Sparse matrix must be able to return CSC format--"\
              "A.getCSC()--or CSR format--A.getCSR()"
    if not hasattr(A,'shape'):
        raise ValueError, "Sparse matrix must be able to return shape (rows,cols) = A.shape"
    if hasattr(A, 'getCSC'):
        ftype, lastel, data, index0, index1 = A.getCSC()
        csc = 1
    else:
        ftype, lastel, data, index0, index1 = A.getCSR()
        csc = 0
    M,N = A.shape
    if (M != N):
        raise ValueError, "Matrix must be square."
    gssv = eval('_superlu.' + ftype + 'gssv')
    return gssv(N,lastel,data,index0,index1,b,csc,permc_spec)[0]
    

def lu_factor(A, permc_spec=2, diag_pivot_thresh=1.0,
              drop_tol=0.0, relax=1, panel_size=10):
    ftype, nnz, data, rowind, colptr = A.getCSC()
    M,N = A.shape
    if (M != N):
        raise ValueError, "Can only factor square matrices."
    gstrf = eval('_superlu.' + ftype + 'gstrf')
    return gstrf(N,nnz,data,rowind,colptr,permc_spec,
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
