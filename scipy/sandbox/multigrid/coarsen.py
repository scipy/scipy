from scipy import *

import multigridtools
import scipy
import numpy
    
from utils import diag_sparse,inf_norm


def rs_strong_connections(A,theta):
    if not scipy.sparse.isspmatrix_csr(A): raise TypeError('expected sparse.csr_matrix')

    Sp,Sj,Sx = multigridtools.rs_strong_connections(A.shape[0],theta,A.indptr,A.indices,A.data)
    return scipy.sparse.csr_matrix((Sx,Sj,Sp),A.shape)


def rs_interpolation(A,theta=0.25):
    if not scipy.sparse.isspmatrix_csr(A): raise TypeError('expected sparse.csr_matrix')
    
    S = rs_strong_connections(A,theta)

    T = S.T.tocsr()

    Ip,Ij,Ix = multigridtools.rs_interpolation(A.shape[0],\
                                               A.indptr,A.indices,A.data,\
                                               S.indptr,S.indices,S.data,\
                                               T.indptr,T.indices,T.data)

    return scipy.sparse.csr_matrix((Ix,Ij,Ip))


def sa_strong_connections(A,epsilon):
    if not scipy.sparse.isspmatrix_csr(A): raise TypeError('expected sparse.csr_matrix')

    Sp,Sj,Sx = multigridtools.sa_strong_connections(A.shape[0],epsilon,A.indptr,A.indices,A.data)
    return scipy.sparse.csr_matrix((Sx,Sj,Sp),A.shape)

def sa_constant_interpolation(A,epsilon=None):
    if not scipy.sparse.isspmatrix_csr(A): raise TypeError('expected sparse.csr_matrix')
    
    if epsilon is not None:
        S = sa_strong_connections(A,epsilon)
    else:
        S = A
    
    #tentative (non-smooth) interpolation operator I
    Ij = multigridtools.sa_get_aggregates(A.shape[0],S.indptr,S.indices)
    Ip = numpy.arange(len(Ij)+1)
    Ix = numpy.ones(len(Ij))
    
    return scipy.sparse.csr_matrix((Ix,Ij,Ip))


def sa_interpolation(A,epsilon,omega=4.0/3.0):
    if not scipy.sparse.isspmatrix_csr(A): raise TypeError('expected sparse.csr_matrix')
    
    I = sa_constant_interpolation(A,epsilon)

    D_inv = diag_sparse(1.0/diag_sparse(A))       
    
    D_inv_A  = D_inv * A
    D_inv_A *= -omega/inf_norm(D_inv_A)

    P = I + (D_inv_A*I)  #same as P=S*I, (faster?)
        
    return P



