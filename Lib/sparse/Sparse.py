from scipy.numeric import *
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

# a = spMatrix(M,N,nzmax,typecode=Float)
#      Create a MxN matrix with room for nzmax non-zero elements of type typecode
# a = spMatrix(data,row,col{,M,N,nzmax}) 
class spMatrix:
    def __init__(self,s,i=None,j=None,M=None,N=None,nzmax=None,typecode=Float):
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
            nzmax = int(j)
            self.ptype = typecode
            self.ftype = _transtabl[self.ptype]
            self.data = zeros((nzmax,),typecode)
            self.index = [zeros((nzmax,)),zeros((M+1,))]
            self.lastel = 0
        elif isspMatrix(s) and s.storage=='CSR':  # make a copy
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
        return "<%dx%d spMatrix of type '%s' with %d elements in %s>" % (self.shape + (self.ptype, self.nzmax, _formats[self.storage][1]))

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
                b = spMatrix(nr[0],nc[0],nels)
                b.lastel = nels-1
                b.data = ao
                b.index = [jao,iao]
                return b                                   
        raise TypeError, "Cannot access sparse matrix that way."


    def astype(self,newtype):
        if newtype == self.ptype:
            return self
        else:
            b = spMatrix(self)
            b.data = b.data.astype(newtype)
            b.ptype = newtype
            b.ftype = _transtabl[newtype]
            return b

    def __add__(self,other):
        if not isspMatrix(other):
            raise TypeError, "Both matrices must be sparse."
        spadd = eval('_sparsekit.'+self.ftype+'aplb')
        assert self.shape == other.shape
        assert self.storage == 'CSR'
        if other.ftype != self.ftype:
            other = other.astype(self.ptype)
        new = spMatrix(self.shape[0],self.shape[1],min((self.nzmax + other.nzmax,product(self.shape))),typecode=self.ptype)
        ierr = array(0)
        iw = zeros((self.shape[1],))
        spadd(array(self.shape[0]),array(self.shape[1]),array(1),self.data,self.index[0],self.index[1],other.data,other.index[0],other.index[1],new.data,new.index[0],new.index[1],array(new.nzmax),iw,ierr)
        nels = max(new.index[1])-1
        new.data = array(new.data[:nels],copy=1)
        new.index[0] = array(new.index[0][:nels],copy=1)
        new.lastel = nels - 1
        new.nzmax = nels
        return new

    def __neg__(self):
        new = spMatrix(self.shape[0],self.shape[1],self.nzmax)
        new.data = -self.data
        new.index = self.index
        new.ptype = self.ptype
        new.ftype = self.ftype
        new.lastel = self.lastel
        return new

    def __sub__(self,other):
        if not isspMatrix(other):
            raise TypeError, "Right operand must also be sparse."
        return self + (-other)

    def __mul__(self,other):
        if isspMatrix(other):
            assert other.shape[0] == self.shape[1]
            assert self.storage == 'CSR'
            new_nz = self.nzmax + other.nzmax
            new = spMatrix(self.shape[0],other.shape[1],new_nz,typecode=self.ptype)
            mult = eval('_sparsekit.'+self.ftype+'amub')
            iw = zeros((self.shape[1],))
            ierr = array(0)
            while 1:  # mult returns error if array wasn't big enough
                mult(array(self.shape[0]),array(other.shape[1]),array(1),
                     self.data, self.index[0], self.index[1], other.data,
                     other.index[0], other.index[1], new.data, new.index[0],
                     new.index[1], array(new.nzmax), iw, ierr)
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
            other = array(other,copy=0,typecode=self.ptype)
            assert len(other.shape)==1 and len(other) == self.shape[1]
            matvec = eval('_sparsekit.'+self.ftype+'amux')
            y = zeros((self.shape[0]),self.ptype)
            matvec(array(self.shape[0]),other,y,self.data,self.index[0],self.index[1])
            return y
        
        elif type(other) in [types.IntType, types.FloatType, types.ComplexType]:
            new = spMatrix(self)           # make a copy
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
            new = spMatrix(self) # make a copy
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
            new = spMatrix(self)
        else:
            new = self
        new.data = conjugate(self.data)
        return new
            
def isspMatrix(x):
    return hasattr(x,'__class__') and x.__class__ is spMatrix


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
    diags = array(diags,copy=0)
    if diags.typecode() not in 'fdFD':
        diags = diags.astype('d')
    offsets = array(offsets,copy=0)
    mtype = diags.typecode()
    assert(len(offsets) == diags.shape[0])
    # set correct diagonal to csr conversion routine for this type
    diagfunc = eval('_sparsekit.'+_transtabl[mtype]+'diacsr')
    # construct empty sparse Matrix and pass it's main parameters to
    #  the diagonal to csr conversion routine.
    nzmax = diags.shape[0]*diags.shape[1]
    s = spMatrix(m,n,nzmax,typecode=mtype)
    diagfunc(array(m), array(n), array(0), array(diags.shape[0]), diags,
             array(diags.shape[1]), offsets, s.data, s.index[0], s.index[1])

    # compute how-many elements were actually filled
    s.lastel = min([m,n])*len(offsets) - 1 - \
               sum(_spdiags_tosub(offsets, a=min([n-m,0]), b=max([n-m,0])))
    return s

def sparse_linear_solve(A,b):
    assert isspMatrix(A)
    assert A.storage=='CSR'
    gssv = eval('_superlu.' + A.ftype + 'gssv')
    return gssv(A.shape[0],A.shape[1],A.lastel+1,A.data,A.index[0]-1,A.index[1]-1,b)
    

splinsolve = sparse_linear_solve

if __name__ == "__main__":
    a = spMatrix(arange(1,9),[0,1,1,2,2,3,3,4],[0,1,3,0,2,3,4,4])
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
