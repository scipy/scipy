__all__ =['approximate_spectral_radius','infinity_norm','diag_sparse',
          'hstack_csr','vstack_csr']

import numpy
import scipy
from scipy import ravel, arange, concatenate, tile, asarray, sqrt, diff, \
                  rand, zeros, ones, empty, asmatrix, dot
from scipy.linalg import norm, eigvals
from scipy.sparse import isspmatrix, isspmatrix_csr, isspmatrix_csc, \
        isspmatrix_bsr, csr_matrix, csc_matrix, bsr_matrix, coo_matrix, \
        extract_diagonal
from scipy.sparse.sputils import upcast


def approximate_spectral_radius(A,tol=0.1,maxiter=10,symmetric=None):
    """approximate the spectral radius of a matrix

    *Parameters*:
        A : dense or sparse matrix 
            E.g. csr_matrix, csc_matrix, ndarray, etc.

        tol : {scalar}
            Tolerance of approximation

        maxiter : {integer}
            Maximum number of iterations to perform

        symmetric : {None,boolean}
            True  - if A is symmetric
                    Lanczos iteration is used (more efficient)
            False - if A is non-symmetric
                    Arnoldi iteration is used (less efficient)
            None  - symmetry of A unknown
                    Method chosen automatically (default)
    *Returns*:
        An approximation to the spectral radius of A (scalar value)

    """
    #from scipy.sandbox.arpack import eigen
    #return norm(eigen(A, k=1, ncv=min(10,A.shape[0]), which='LM', tol=tol, return_eigenvectors=False))
   
    if not isspmatrix(A):
        A = asmatrix(A) #convert dense arrays to matrix type
    
    if A.shape[0] != A.shape[1]:
        raise ValueError,'expected square matrix'

    #TODO make method adaptive

    numpy.random.seed(0)  #make results deterministic

    v0  = rand(A.shape[1],1)
    v0 /= norm(v0)

    H  = zeros((maxiter+1,maxiter))
    V = [v0]

    #save past estimates
    #estimates = []

    for j in range(maxiter):
        w = A * V[-1]
   
        if symmetric:
            if j >= 1:
                H[j-1,j] = beta
                w -= beta * V[-2]

            alpha = dot(ravel(w),ravel(V[-1]))
            H[j,j] = alpha
            w -= alpha * V[-1]
            
            beta = norm(w)
            if (H[j+1,j] < 1e-10): break
            
            w /= beta
            H[j+1,j] = beta

            V.append(w)
            V = V[-2:] #retain only last two vectors
        else:
            #orthogonalize against Vs
            for i,v in enumerate(V):
                H[i,j] = dot(ravel(w),ravel(v))
                w -= H[i,j]*v
            H[j+1,j] = norm(w)
            if (H[j+1,j] < 1e-10): break
            
            w /= H[j+1,j] 
            V.append(w)
   
            # if upper 2x2 block of Hessenberg matrix H is almost symmetric,
            # and the user has not explicitly specified symmetric=False,
            # then switch to symmetric Lanczos algorithm
            if symmetric is not False and j == 1:
                if abs(H[1,0] - H[0,1]) < 1e-12:
                    symmetric = True
                    V = V[1:]
                    H[1,0] = H[0,1]
                    beta = H[2,1]
    
    return max([norm(x) for x in eigvals(H[:j+1,:j+1])])      



def infinity_norm(A):
    """
    Infinity norm of a sparse matrix (maximum absolute row sum).  This serves
    as an upper bound on spectral radius.
    """

    if isspmatrix_csr(A) or isspmatrix_csc(A):
        #avoid copying index and ptr arrays
        abs_A = A.__class__((abs(A.data),A.indices,A.indptr),shape=A.shape)
        return (abs_A * ones(A.shape[1],dtype=A.dtype)).max()
    else:
        return (abs(A) * ones(A.shape[1],dtype=A.dtype)).max()

def diag_sparse(A):
    """
    If A is a sparse matrix (e.g. csr_matrix or csc_matrix)
       - return the diagonal of A as an array

    Otherwise
       - return a csr_matrix with A on the diagonal
    """

    #TODO integrate into SciPy?
    if isspmatrix(A):
        return extract_diagonal(A)
    else:
        return csr_matrix((asarray(A),arange(len(A)),arange(len(A)+1)),(len(A),len(A)))

def scale_rows(A,v,copy=True):
    from scipy.sparse.sparsetools import csr_scale_rows, bsr_scale_rows

    v = ravel(v)

    if isspmatrix_csr(A) or isspmatrix_bsr(A):
        M,N = A.shape
        if M != len(v):
            raise ValueError,'scale vector has incompatible shape'

        if copy:
            A = A.copy()
            A.data = asarray(A.data,dtype=upcast(A.dtype,v.dtype))
        else:
            v = asarray(v,dtype=A.dtype)

        if isspmatrix_csr(A):
            csr_scale_rows(M, N, A.indptr, A.indices, A.data, v)
        else:
            R,C = A.blocksize
            bsr_scale_rows(M/R, N/C, R, C, A.indptr, A.indices, ravel(A.data), v)

        return A
    elif isspmatrix_csc(A):
        return scale_columns(A.T,v)
    else:
        return scale_rows(csr_matrix(A),v)
        
def scale_columns(A,v,copy=True):
    from scipy.sparse.sparsetools import csr_scale_columns, bsr_scale_columns

    v = ravel(v)

    if isspmatrix_csr(A) or isspmatrix_bsr(A):
        M,N = A.shape
        if N != len(v):
            raise ValueError,'scale vector has incompatible shape'

        if copy:
            A = A.copy()
            A.data = asarray(A.data,dtype=upcast(A.dtype,v.dtype))
        else:
            v = asarray(v,dtype=A.dtype)

        if isspmatrix_csr(A):
            csr_scale_columns(M, N, A.indptr, A.indices, A.data, v)
        else:
            R,C = A.blocksize
            bsr_scale_columns(M/R, N/C, R, C, A.indptr, A.indices, ravel(A.data), v)

        return A
    elif isspmatrix_csc(A):
        return scale_rows(A.T,v)
    else:
        return scale_rows(csr_matrix(A),v)

def symmetric_rescaling(A,copy=True):
    if isspmatrix_csr(A) or isspmatrix_csc(A) or isspmatrix_bsr(A):
        if A.shape[0] != A.shape[1]:
            raise ValueError,'expected square matrix'

        D = diag_sparse(A)
        mask = D == 0

        D_sqrt = sqrt(abs(D))
        D_sqrt_inv = 1.0/D_sqrt
        D_sqrt_inv[mask] = 0

        DAD = scale_rows(A,D_sqrt_inv,copy=copy)
        DAD = scale_columns(DAD,D_sqrt_inv,copy=False)

        return D_sqrt,D_sqrt_inv,DAD

    else:
        return symmetric_rescaling(csr_matrix(A))


def hstack_csr(A,B):
    if not isspmatrix(A) or not isspmatrix(B):
        raise TypeError,'expected sparse matrix'

    if A.shape[0] != B.shape[0]:
        raise ValueError,'row dimensions must agree'

    A = A.tocoo()
    B = B.tocoo()
    I = concatenate((A.row,B.row))
    J = concatenate((A.col,B.col+A.shape[1]))
    V = concatenate((A.data,B.data))
    return coo_matrix((V,(I,J)),shape=(A.shape[0],A.shape[1]+B.shape[1])).tocsr()

def vstack_csr(A,B):
    #TODO OPTIMIZE THIS
    if not isspmatrix(A) or not isspmatrix(B):
        raise TypeError,'expected sparse matrix'

    if A.shape[1] != B.shape[1]:
        raise ValueError,'column dimensions must agree'

    A = A.tocoo()
    B = B.tocoo()
    I = concatenate((A.row,B.row+A.shape[0]))
    J = concatenate((A.col,B.col))
    V = concatenate((A.data,B.data))
    return coo_matrix((V,(I,J)),shape=(A.shape[0]+B.shape[0],A.shape[1])).tocsr()








































