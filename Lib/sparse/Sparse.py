from Numeric import *
from scipy_base.fastumath import *
import types
import _sparsekit
import _superlu

MAXPRINT=50
# The formats that SPARSEKIT's convert programs understand.

_formats = {'BND':[0,"Linpack Banded format"],
            'CSR':[1,"Compressed Sparse Row format"],
            'CSC':[2,"Compressed Sparse Column format"],
            'COO':[3,"Coordinate format"],
            'ELL':[4,"Ellpack-Itpack genralized diagonal format"],
            'DIA':[5,"Diagonal format"],
            'BSR':[6,"Block Sparse Row format"],
            'MSR':[7,"Modified Compressed Sparse Row format"],
            'SSK':[8,"Symmetric Skyline format"],
            'NSK':[9,"Nonsymmetric Skyline format"],
            'LNK':[10,"Linked list storage format"],
            'JAD':[11,"Jagged Diagonal format"],
            'SSS':[12,"The Symmetric Sparse Skyline format"],
            'USS':[13,"The Unsymmetric Sparse Skyline format"],
            'VBR':[14,"Variable Block Row format"]}

# So far only CSR format is supported internally.

_transtabl = {'f':'s','d':'d','F':'c','D':'z'}

# A sparse matrix class.  A sparse matrix can be initialized as:

# a = spmatrix(M,N,nzmax,typecode=Float)
#   Create an MxN matrix with room for nzmax non-zero elements of
#     type typecode
# a = spmatrix(data,row,col{,M,N,nzmax})
#   Create a sparse matrix with non-zero elements data
#     using a[row[i],col[i]] = data[i]


#   A simple "dictionary-based" sparse matrix.
#   keys must be 2-tuples of any object.  The object must have __int__ defined to return integer
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
    def __init__(self):
        dict.__init__(self)
        self.shape = (0,0)
        self.storage = 'dict'
        self.data = {}

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
        for key in other.keys():
            try:
                res[key] += other[key]
            except KeyError:
                res[key] = other[key]
        return res

    def __sub__(self, other):
        res = dictmatrix()
        res.update(self)
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
        if (isinstance(other, list) or (isinstance(other, ArrayType) and rank(other)==1) and \
            len(other)==self.shape[1]):
            return self.matvec(other)
        
        res = dictmatrix()
        for key in self.keys():
            res[key] = other * self[key]
        return res

    def __len__(self):
        return len(self.keys())

    def take(self, cols_or_rows, columns=1):
        # Extract columns or rows as indictated from matrix
        res = dictmatrix()
        indx = int((columns == 1))
        for key in self.keys():
            if key[indx] in cols_or_rows:
                res[key] = self[key]
        return res

    def matvec(self, other):
        if len(other) != self.shape[1]:
            raise ValueError, "Dimensions do not match."
        keys = self.keys()
        res = [0]*self.shape[0]
        for key in keys:
            res[int(key[0])] += self[key] * other[int(key[1])]
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
        data = [None]*nnz
        colind = [None]*nnz
        row_ptr = [None]*(self.shape[0]+1)
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
        ftype = data.typecode()
        if ftype not in ['d','D','f','F']:
            data = data*1.0
            ftype = 'd'
        return ftype, nnz, data, colind, row_ptr

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
        ftype = data.typecode()
        if ftype not in ['d','D','f','F']:
            data = data*1.0
            ftype = 'd'
        return ftype, nnz, data, rowind, col_ptr

               
class spmatrix:
    def __init__(self,s,i=None,j=None,M=None,N=None,nzmax=None,
                 typecode=Float):
        if type(s) in [types.ListType, ArrayType]:
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
        return A.ftype, A.lastel+1, A.data, A.index[0]-1, A.index[1]-1
            
def isspmatrix(x):
    return hasattr(x,'__class__') and x.__class__ is spmatrix


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

def sparse_linear_solve(A,b,permc_spec=0):
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
    gssv = eval('_superlu.' + ftype + 'gssv')
    return gssv(M,N,lastel,data,index0,index1,b,csc,permc_spec)
    

splinsolve = sparse_linear_solve
solve = splinsolve

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
    x = splinsolve(a,b)
    print x
    print "Error: ",a*x-b

    print "Solve: double precision complex."
    a = a.astype('D')
    x = splinsolve(a,b)
    print x
    print "Error: ",a*x-b

    print "Solve: double precision."
    a = a.astype('d')
    x = splinsolve(a,b)
    print x
    print "Error: ",a*x-b

    print "Solve: single precision."
    a = a.astype('f')
    x = splinsolve(a,b.astype('f'))
    print x
    print "Error: ",a*x-b
