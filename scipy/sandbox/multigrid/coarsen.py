
import multigridtools
import scipy
import numpy
    
from utils import diag_sparse,infinity_norm


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

def sa_constant_interpolation(A,epsilon):
    if not scipy.sparse.isspmatrix_csr(A): raise TypeError('expected sparse.csr_matrix')
    
    S = sa_strong_connections(A,epsilon)

    #S.ensure_sorted_indices()

    #tentative (non-smooth) interpolation operator I
    Pj = multigridtools.sa_get_aggregates(S.shape[0],S.indptr,S.indices)
    Pp = numpy.arange(len(Pj)+1)
    Px = numpy.ones(len(Pj))
    
    return scipy.sparse.csr_matrix((Px,Pj,Pp))

##def sa_smoother(A,S,omega):
##    Bp,Bj,Bx = multigridtools.sa_smoother(A.shape[0],omega,A.indptr,A.indices,A.data,S.indptr,S.indices,S.data)
##
##    return csr_matrix((Bx,Bj,Bp),dims=A.shape)
    
def sa_interpolation(A,epsilon,omega=4.0/3.0):
    if not scipy.sparse.isspmatrix_csr(A): raise TypeError('expected sparse.csr_matrix')
   
    P  = sa_constant_interpolation(A,epsilon)

##    As = sa_strong_connections(A,epsilon)
##    S  = sa_smoother(A,S,omega)
    

    D_inv = diag_sparse(1.0/diag_sparse(A))       
    
    D_inv_A  = D_inv * A
    D_inv_A *= omega/infinity_norm(D_inv_A)

    I = P - (D_inv_A*P)  #same as I=S*P, (faster?)
           
    return I



