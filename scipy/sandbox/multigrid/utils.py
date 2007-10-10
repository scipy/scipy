__all__ =['approximate_spectral_radius','infinity_norm','diag_sparse']

import numpy
import scipy
from numpy import ravel,arange
from scipy.linalg import norm
from scipy.sparse import isspmatrix,isspmatrix_csr,isspmatrix_csc, \
                        csr_matrix,csc_matrix,extract_diagonal


def approximate_spectral_radius(A,tol=0.1,maxiter=20):
    """
    Approximate the spectral radius of a symmetric matrix using ARPACK
    """
    from scipy.sandbox.arpack import eigen
    return norm(eigen(A, k=1, ncv=10, which='LM', maxiter=maxiter, tol=tol, return_eigenvectors=False)[0])



def infinity_norm(A):
    """
    Infinity norm of a sparse matrix (maximum absolute row sum).  This serves 
    as an upper bound on spectral radius.
    """
    
    if isspmatrix_csr(A) or isspmatrix_csc(A):
        #avoid copying index and ptr arrays
        abs_A = A.__class__((abs(A.data),A.indices,A.indptr),dims=A.shape,check=False)
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
    
    if isspmatrix(A):
        return extract_diagonal(A)
    else:
        return csr_matrix((A,arange(len(A)),arange(len(A)+1)),(len(A),len(A)))


