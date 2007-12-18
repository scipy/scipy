"""Sparse DIAgonal format"""

__all__ = ['dia_matrix','isspmatrix_dia']

from numpy import asarray, asmatrix, matrix, zeros, arange, unique, \
        searchsorted, intc, atleast_1d, atleast_2d

from base import spmatrix, isspmatrix, isdense
from sputils import isscalarlike, isshape, upcast, getdtype

from sets import Set

class dia_matrix(spmatrix):
    #TODO add description/motivation
    def __init__(self, arg1, shape=None, dtype=None, copy=False):
        spmatrix.__init__(self)

        if isdense(arg1) or isspmatrix(arg1):
            #convert to COO first, then to DIA
            from coo import coo_matrix
            if isdense(arg1):
                A = coo_matrix(arg1)
            elif arg1.format in ['csr','csc','coo']:
                A = arg1.tocoo(copy=False)
            else:
                A = arg1.tocoo()
            
            ks = A.col - A.row  #the offset for each nonzero          
            diags = unique(ks)

            if len(diags) > 10:
                pass #do something

            #initialize and fill in data array
            data_shape = (len(diags), A.col.max()+1)
            data_dtype = getdtype(dtype, default=A.data.dtype)
            data = zeros(data_shape, dtype=data_dtype)
            data[ searchsorted(diags,ks), A.col ] = A.data
  
            self.data,self.diags = data,diags
            self.shape = A.shape
        elif isinstance(arg1, tuple):
            if isshape(arg1):
                # It's a tuple of matrix dimensions (M, N)
                # create empty matrix
                self.shape = arg1   #spmatrix checks for errors here
                self.data   = zeros( (0,0), getdtype(dtype, default=float))
                self.diags = zeros( (0), dtype=intc)
            else:
                try:
                    # Try interpreting it as (data, diags)
                    data, diags = arg1
                except:
                    raise ValueError, "unrecognized form for dia_matrix constructor"
                else:
                    if shape is None:
                        raise ValueError,'expected a shape argument'
                    self.data   = atleast_2d(asarray(arg1[0],dtype=dtype))
                    self.diags = atleast_1d(asarray(arg1[1],dtype='i'))
                    self.shape   = shape

        #check format
        if self.diags.ndim != 1:
            raise ValueError,'diags array must have rank 1'

        if self.data.ndim != 2:
            raise ValueError,'data array must have rank 2'

        if self.data.shape[0] != len(self.diags):
            raise ValueError,'number of diagonals (%d) ' \
                    'does not match the number of diags (%d)' \
                    % (self.data.shape[0], len(self.diags))
        
        if len(Set(self.diags)) != len(self.diags):
            raise ValueError,'offset array contains duplicate values'



    def __getdtype(self):
        return self.data.dtype

    def __setdtype(self,newtype):
        self.dtype = newtype

    dtype = property(fget=__getdtype,fset=__setdtype)

    def getnnz(self):
        """number of nonzero values

        explicit zero values are included in this number
        """
        M,N = self.shape
        nnz = 0
        for k in self.diags:
            if k > 0:
                nnz += min(M,N-k)
            else:
                nnz += min(M+k,N)
        return nnz

    nnz = property(fget=getnnz)



    def __mul__(self, other): # self * other
        """ Scalar, vector, or matrix multiplication
        """
        if isscalarlike(other):
            return dia_matrix((other * self.data, self.diags),self.shape)
        else:
            return self.dot(other)


    def matmat(self, other):
        if isspmatrix(other):
            M, K1 = self.shape
            K2, N = other.shape
            if (K1 != K2):
                raise ValueError, "shape mismatch error"

            return self.tocsr() * other
            #TODO handle sparse matmat here
        else:
            # matvec handles multiple columns
            return self.matvec(other)


    def matvec(self,other):
        #TODO why is this so slow? it should run at BLAS-speed

        x = asarray(other)

        if x.ndim == 1:
            x = x.reshape(-1,1)
        if self.shape[1] != x.shape[0]:
            raise ValueError, "dimension mismatch"

        y = zeros((self.shape[0],x.shape[1]), dtype=upcast(self.dtype,x.dtype))

        L = self.data.shape[1]
        M,N = self.shape

        for diag,k in zip(self.data,self.diags):
            j_start = max(0,k)
            j_end   = min(M+k,N,L)

            i_start = max(0,-k)
            i_end   = i_start + (j_end - j_start)

            y[i_start:i_end] += diag[j_start:j_end].reshape(-1,1) * x[j_start:j_end]

        
        if isinstance(other, matrix):
            y = asmatrix(y)

        if other.ndim == 1:
            # If 'other' was a 1d array, reshape the result
            y = y.reshape(-1)

        return y

    def todia(self,copy=False):
        if copy:
            return self.copy()
        else:
            return self

    def tocsr(self):
        #TODO optimize COO->CSR
        return self.tocoo().tocsr()

    def tocsc(self):
        #TODO optimize COO->CSC
        return self.tocoo().tocsc()

    def tocoo(self):
        num_data = len(self.data)
        len_data = self.data.shape[1]
        
        row = arange(len_data).reshape(1,-1).repeat(num_data,axis=0)
        col = row.copy()

        for i,k in enumerate(self.diags):
            row[i,:] -= k
        
        mask = (row >= 0) & (row < self.shape[0]) & (col < self.shape[1])
        row,col,data = row[mask],col[mask],self.data[mask]
        row,col,data = row.reshape(-1),col.reshape(-1),data.reshape(-1)
       
        from coo import coo_matrix
        return coo_matrix((data,(row,col)),shape=self.shape)



from sputils import _isinstance

def isspmatrix_dia(x):
    return _isinstance(x, dia_matrix)


