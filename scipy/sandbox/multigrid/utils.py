__all__ =['approximate_spectral_radius','infinity_norm','diag_sparse',
          'hstack_csr','vstack_csr','expand_into_blocks']

import numpy
import scipy
from scipy import ravel,arange,concatenate,tile,asarray,sqrt,diff, \
                  rand,zeros,empty,asmatrix,dot
from scipy.linalg import norm,eig
from scipy.sparse import isspmatrix,isspmatrix_csr,isspmatrix_csc, \
                        csr_matrix,csc_matrix,extract_diagonal, \
                        coo_matrix


def approximate_spectral_radius(A,tol=0.1,maxiter=8):
    """
    Approximate the spectral radius of a matrix
    """
    #from scipy.sandbox.arpack import eigen
    #return norm(eigen(A, k=1, ncv=10, which='LM', maxiter=maxiter, tol=tol, return_eigenvectors=False))
   
    numpy.random.seed(0)  #make results deterministic

    #TODO profile vs V -> V.T
    #TODO make algorithm adaptive

    if not isspmatrix(A):
        #convert dense arrays to matrix type
        A = asmatrix(A)

    v0 = rand(A.shape[0])
    H  = zeros((maxiter+1,maxiter))
    V  = zeros((A.shape[0],maxiter+1))
 
    V[:,0]= v0/norm(v0)
    for j in range(maxiter):
        w = A * V[:,j]
        for i in range(j+1):
            H[i,j] = dot(w,V[:,i])
            w -= H[i,j]*V[:,i]
        H[j+1,j] = norm(w)
        if (H[j+1,j] < 1e-12): break
        V[:,j+1] = (1.0/H[j+1,j]) * w
    # end
    m=j
    Heig,tmp = eig(H[0:m,0:m])
    return max( [norm(x) for x in Heig] )



def infinity_norm(A):
    """
    Infinity norm of a sparse matrix (maximum absolute row sum).  This serves
    as an upper bound on spectral radius.
    """

    if isspmatrix_csr(A) or isspmatrix_csc(A):
        #avoid copying index and ptr arrays
        abs_A = A.__class__((abs(A.data),A.indices,A.indptr),dims=A.shape)
        return (abs_A * numpy.ones(A.shape[1],dtype=A.dtype)).max()
    else:
        return (abs(A) * numpy.ones(A.shape[1],dtype=A.dtype)).max()

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


def symmetric_rescaling(A):
    if not (isspmatrix_csr(A) or isspmatrix_csc(A)):
        raise TypeError,'expected csr_matrix or csc_matrix'

    if A.shape[0] != A.shape[1]:
        raise ValueError,'expected square matrix'

    D = diag_sparse(A)
    mask = D == 0

    #D[mask] = 0
    D_sqrt = sqrt(abs(D))
    D_sqrt_inv = 1.0/D_sqrt
    D_sqrt_inv[mask] = 0

    #TODO time this against simple implementation
    data = A.data * D_sqrt_inv[A.indices]
    data *= D_sqrt_inv[arange(A.shape[0]).repeat(diff(A.indptr))]

    DAD = A.__class__((data,A.indices,A.indptr),dims=A.shape)

    return D_sqrt,D_sqrt_inv,DAD


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
    return coo_matrix((V,(I,J)),dims=(A.shape[0],A.shape[1]+B.shape[1])).tocsr()

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
    return coo_matrix((V,(I,J)),dims=(A.shape[0]+B.shape[0],A.shape[1])).tocsr()


def expand_into_blocks(A,m,n):
    """Expand each element in a sparse matrix A into an m-by-n block.

          Example:
          >>> A.todense()
          matrix([[ 1.,  2.],
                  [ 4.,  5.]])

          >>> expand_into_blocks(A,2,2).todense()
          matrix([[ 1.,  1.,  2.,  2.],
                  [ 1.,  1.,  2.,  2.],
                  [ 4.,  4.,  5.,  5.],
                  [ 4.,  4.,  5.,  5.]])

    """
    #TODO EXPLAIN MORE
    #TODO use spkron instead, time for compairson

    if m == 1 and n == 1:
        return A #nothing to do

    A = A.tocoo()

    # expand 1x1 -> mxn
    row  = ( m*A.row ).repeat(m*n).reshape(-1,m,n)
    col  = ( n*A.col ).repeat(m*n).reshape(-1,m,n)

    # increment indices
    row += tile(arange(m).reshape(-1,1),(1,n))
    col += tile(arange(n).reshape(1,-1),(m,1))

    # flatten
    row = row.reshape(-1)
    col = col.reshape(-1)

    data = A.data.repeat(m*n)

    return coo_matrix((data,(row,col)),dims=(m*A.shape[0],n*A.shape[1]))
