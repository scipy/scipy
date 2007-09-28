import multigridtools
import scipy,numpy,scipy.sparse
from scipy.sparse import csr_matrix,isspmatrix_csr

from utils import diag_sparse,approximate_spectral_radius


def rs_strong_connections(A,theta):
    """
    Return a strength of connection matrix using the method of Ruge and Stuben

        An off-diagonal entry A[i.j] is a strong connection iff

                -A[i,j] >= theta * max( -A[i,k] )   where k != i
    """
    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')
    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')

    Sp,Sj,Sx = multigridtools.rs_strong_connections(A.shape[0],theta,A.indptr,A.indices,A.data)
    return scipy.sparse.csr_matrix((Sx,Sj,Sp),A.shape)


def rs_interpolation(A,theta=0.25):
    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')
    
    S = rs_strong_connections(A,theta)

    T = S.T.tocsr()  #transpose S for efficient column access

    Ip,Ij,Ix = multigridtools.rs_interpolation(A.shape[0],\
                                               A.indptr,A.indices,A.data,\
                                               S.indptr,S.indices,S.data,\
                                               T.indptr,T.indices,T.data)

    return scipy.sparse.csr_matrix((Ix,Ij,Ip))


def sa_strong_connections(A,epsilon):
    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')

    Sp,Sj,Sx = multigridtools.sa_strong_connections(A.shape[0],epsilon,A.indptr,A.indices,A.data)

    return scipy.sparse.csr_matrix((Sx,Sj,Sp),A.shape)

def sa_constant_interpolation(A,epsilon,blocks=None):
    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')
    
    if blocks is not None:
        num_dofs   = A.shape[0]
        num_blocks = blocks.max()
        
        if num_dofs != len(blocks):
            raise ValueError,'improper block specification'
        
        # for non-scalar problems, use pre-defined blocks in aggregation
        # the strength of connection matrix is based on the Frobenius norms of the blocks
        
        B  = csr_matrix((ones(num_dofs),blocks,arange(num_dofs + 1)),dims=(num_dofs,num_blocks))
        Block_Frob = B.T.tocsr() * csr_matrix((A.data**2,A.indices,A.indptr),dims=A.shape) * B #Frobenius norms of blocks entries of A

        S = sa_strong_connections(Block_Frob,epsilon)
    
        Pj = multigridtools.sa_get_aggregates(S.shape[0],S.indptr,S.indices)
        Pj = Pj[blocks] #expand block aggregates into constituent dofs
        Pp = B.indptr
        Px = B.data
    else:
        S = sa_strong_connections(A,epsilon)
        
        Pj = multigridtools.sa_get_aggregates(S.shape[0],S.indptr,S.indices)
        Pp = numpy.arange(len(Pj)+1)
        Px = numpy.ones(len(Pj))
    
    return scipy.sparse.csr_matrix((Px,Pj,Pp))

    
##    S = sa_strong_connections(A,epsilon)
##
##    #tentative (non-smooth) interpolation operator I
##    Pj = multigridtools.sa_get_aggregates(S.shape[0],S.indptr,S.indices)
##    Pp = numpy.arange(len(Pj)+1)
##    Px = numpy.ones(len(Pj))
##    
##    return scipy.sparse.csr_matrix((Px,Pj,Pp))

##def sa_smoother(A,S,omega):
##    Bp,Bj,Bx = multigridtools.sa_smoother(A.shape[0],omega,A.indptr,A.indices,A.data,S.indptr,S.indices,S.data)
##
##    return csr_matrix((Bx,Bj,Bp),dims=A.shape)
    
def sa_interpolation(A,epsilon,omega=4.0/3.0,blocks=None):
    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')
   
    P  = sa_constant_interpolation(A,epsilon=epsilon,blocks=blocks)

    D_inv = diag_sparse(1.0/diag_sparse(A))       
    
    D_inv_A  = D_inv * A
    D_inv_A *= omega/approximate_spectral_radius(D_inv_A)

    I = P - (D_inv_A*P)  #same as I=S*P, (faster?)
           
    return I



