"""UMFSparse:  A minimal sparse-matrix class for Python, useful for interfacing
with the UMFPACK library for solving large linear systems,

This file defines a class: UMFmatrix that defines a sparse-matrix class.
It is designed so that the attribute data points to the non-zero values in
the matrix and index is a list of arrays which detail how to place the
non-zero values in the array.  Currently only a compressed sparse row
storage is implemented. 
"""

import scipy_base
import scipy.io as io
import types
import struct
import umfpack

MAXPRINT=50
MEMALLOCSIZE=100

# A sparse matrix class.  A sparse matrix can be initialized as:
# a = UMFmatrix(N,nzmax,typecode=Float)
#      Create an NxN matrix with room for nzmax non-zero elements of
#      type typecode 
class UMFmatrix:
    """Dolmatrix defines a sparse-matrix class.

    A = UMFmatrix(N,nzmax=None,typecode=Complex64)
       Initialize an N x N sparse matrix with room for nzmax non-zero
       elements of the given typecode (Defaults to 'D')

    A = UMFmatrix(otherSparse)
       Make a copy of the Sparse matrix (see copy() method)
    """

    def __init__(self,N,nzmax=None,typecode=scipy_base.Complex64,
                 zeroval=1e-16, chunks=MEMALLOCSIZE):
        if isUMFmatrix(N):  # make a copy
            self = N.copy()
            return
        if nzmax is None:
            nzmax = chunks
        self.data = scipy_base.zeros((nzmax,),typecode)
        self.index = scipy_base.zeros((2*nzmax,),scipy_base.Int32)
        self.filled=0
        self.shape = (N,N)
        self.maxprint = MAXPRINT
        self._memallocsize = chunks
        return

    def copy(self):
        s = UMFmatrix(1,0)
        for attr in dir(self):
            if attr not in ['data','index']:
                setattr(s,attr,getattr(self,attr))
        s.data = scipy_base.array(self.data,copy=1)
        s.index = scipy_base.array(self.index,copy=1)
        return s

    def write(self,filename):
        fid = io.fopen(filename,'w')
        fid.write(struct.pack('BB',scipy_base.LittleEndian,
                              scipy_base.array(3,'i').itemsize()))
        if scipy_base.LittleEndian:
            swapchar = '<'
        else:
            swapchar = '>'
        fid.write(struct.pack('%sc3siiiii'%swapchar, self.data.typecode(),
                              self.storage, self.shape[0],
                              self.shape[1], len(self.data), self.maxprint,
                              self._memallocsize))
        fid.fwrite(len(self.data), self.data)
        fid.fwrite(2*len(self.data), self.index[0])

    def read(self,filename):
        fid = io.fopen(filename,'r')
        endianness,lsize = struct.unpack('BB',fid.read(2))
        byteswap = (endianness != scipy_base.LittleEndian)
        tsize = scipy_base.array([1],'i').itemsize()
        if lsize != tsize:
            raise ValueError, "Stored size of an int, %d, is not the same as current size of an int, %d." % (lsize, tsize)
        str = fid.read(4+5*lsize)
        if endianness:
            swapchar = '<'
        else:
            swapchar = '>'
        typecode, self.storage, shape0, shape1, nzmax, self.maxprint, self._memallocsize = struct.unpack('%sc3siiiii'%swapchar, str)
        self.shape = (shape0,shape1)
        self.data = fid.fread(nzmax, typecode, typecode, byteswap)
        self.index = [None]*2
        self.index = fid.fread(2*nzmax, 'i', 'i', byteswap)
        fid.close()
        return

    def rowcol(self,key):
        numels = self.filled
        nmax = len(self.data)
        if key >= numels:
            raise ValueError, "There are only %d nonzero entries. " % numels \
                  + "You asked for element %d." % key
        row = self.index[key]
        col = self.index[key+nmax]
        return row,col

    def listprint(self,start,stop):
        val = ''
        for ind in xrange(start,stop):
            val = val + '  %s\t%s\n' % (self.rowcol(ind),self.data[ind])
        return val

    def __repr__(self):
        return "<%dx%d UMFmatrix of type '%s' with %d elements>" % \
               (self.shape + (self.data.typecode(), len(self.data)))

        
    def __str__(self):
        val = ''
        numels = self.filled
        if numels > self.maxprint:
            val = val + self.listprint(0,self.maxprint/2)
            val = val + "  :\t:\n"
            val = val + self.listprint(numels-1-self.maxprint/2,numels)
        else:
            val = val + self.listprint(0,numels)
        return val[:-1]

    def __cmp__(self,other):
        raise TypeError, "Comparison of sparse matrices is not implemented."

    def __len__(self):
        return len(self.data)

    def __setitem__(self,key,value):
        if type(key) == types.TupleType:
            # Check to see if two integers are the keys
            if len(key)==2 and \
               (type(key[0]),type(key[1]))==(types.IntType,)*2:
                nzmax = len(self.data)
                filled = self.filled
                self.data[filled] = value
                self.index[filled] = key[0]
                self.index[filled+nzmax] = key[1]
                filled = filled + 1
                self.filled = filled
            else:
                if len(key)!=2:
                    raise IndexError, "Too many indices: %s" % (key,)
                else:
                    raise IndexError, "Unsupported Slicing operation: %s" % (key,)
        else:
            raise IndexError, "Incorrect selection operation: %s" % (key,)
                            
    def __getitem__(self,key):
        if key is None:
            return 0.0
        if type(key) == types.IntType:
            return (self.data[key], self.rowcol(key))
        elif type(key) == types.TupleType:
            if len(key) == 2 and (type(key[0]),type(key[1])) == (types.IntType,)*2:

                assert(key[0] >= 0 and key[1] >= 0)

                raise ValueError, "Can't do it."
        raise IndexError, "Cannot access sparse matrix that way."

    def astype(self,newtype):
        if newtype == self.data.typecode():
            return self
        else:
            b = self.copy()
            b.data = b.data.astype(newtype)
            return b

    def __neg__(self):
        new = UMFmatrix(self.shape[0],0)
        new.data = -self.data
        new.index = self.index
        return new

    def conjugate(self,inplace=0):
        if inplace == 0:
            new = UMFmatrix(self)
        else:
            new = self
        new.data = conjugate(self.data)
        return new

def isUMFmatrix(x):
    return isinstance(x,UMFmatrix)

def splinsolve(A,b,factored=0,fulloutput=0,lvalcoef=(6,70,1)):

    if not factored:
        assert(isUMFmatrix(A))
        #print A.shape[1], b.shape[0]
        assert(A.shape[1] == b.shape[0])
        assert(A.shape[0] == A.shape[1])
        keep = scipy_base.zeros((20,),'i')
        icntl = scipy_base.zeros((20,),'i')
        cntl = scipy_base.zeros((10,),'d')
        # Call initialization routine.
        umfpack.umz21i(keep,cntl,icntl)

        # Call sparse factorization routine.
        n = A.shape[0]
        ne = A.filled
        nzmax = len(A.data)
        lvalue = lvalcoef[0]*nzmax + lvalcoef[1]*n + lvalcoef[2]
        lindex = 6*nzmax + 70*n + 1
        one = scipy_base.array(1,'i')
        value = scipy_base.zeros((lvalue,),'D')
        index = scipy_base.zeros((lindex,),'i')
        value[:A.filled] = A.data[:A.filled]
        index[:A.filled] = A.index[:A.filled]+one
        index[A.filled:2*A.filled] = A.index[nzmax:nzmax+A.filled]+one
        info = scipy_base.zeros((40,),'i')
        rinfo = scipy_base.zeros((20,),'d')
        umfpack.umz2fa(n,ne,0,0,value, index,
                       keep, cntl, icntl, info, rinfo)
        # Call sparse linear solve routine.
        x = scipy_base.zeros((n,),'D')
        w = scipy_base.zeros((2*n,),'D')
        if info[0] < 0:
            print "Problem with factorization: %d" % info[0]
            return x, info[:24], rinfo[:8]
        umfpack.umz2so(0,0,value,index,keep,b,x,w,cntl,icntl,info,rinfo)

    else:
        assert(len(b.shape)==1)
        n = len(b)
        x = scipy_base.zeros((n,),'D')
        value, index, keep, w, cntl, icntl, info, rinfo = A
        assert(2*n==len(w))
        umfpack.umz2so(0,0,value,index,keep,b,x,w,cntl,icntl,info,rinfo)

    if fulloutput:
        return x, (value, index, keep, w, cntl, icntl, info, rinfo)
    else:
        return x,info[:24],rinfo[:8]


if __name__ == "__main__":
    print "From scratch: 10 x 10 with room for 10 elements initially."
    print ">>> a = UMFSparse.UMFmatrix(10,10)"
    b = UMFmatrix(10,10) 
    print "Representation of a matrix:"
    print repr(b)
    print "How a matrix prints."
    print b
    print "You can set elements of the matrix: "
    b[20,40] = 10

    








