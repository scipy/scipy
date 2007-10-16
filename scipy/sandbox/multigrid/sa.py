import scipy
import numpy
from numpy import array,arange,ones,zeros,sqrt,isinf,asarray,empty,diff,\
                  ascontiguousarray
from scipy.sparse import csr_matrix,isspmatrix_csr

from utils import diag_sparse,approximate_spectral_radius
import multigridtools

__all__ = ['sa_filtered_matrix','sa_strong_connections','sa_constant_interpolation',
           'sa_interpolation','sa_fit_candidates']


def sa_filtered_matrix(A,epsilon,blocks=None):
    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')

    if epsilon == 0:
        A_filtered = A

    else:
        if blocks is None:
            Sp,Sj,Sx = multigridtools.sa_strong_connections(A.shape[0],epsilon,A.indptr,A.indices,A.data)
            A_filtered = csr_matrix((Sx,Sj,Sp),A.shape)
        else:
            A_filtered = A  #TODO subtract weak blocks from diagonal blocks?

##            num_dofs   = A.shape[0]
##            num_blocks = blocks.max() + 1
##            
##            if num_dofs != len(blocks):
##                raise ValueError,'improper block specification'
##            
##            # for non-scalar problems, use pre-defined blocks in aggregation
##            # the strength of connection matrix is based on the 1-norms of the blocks
##            
##            B  = csr_matrix((ones(num_dofs),blocks,arange(num_dofs + 1)),dims=(num_dofs,num_blocks))
##            Bt = B.T.tocsr()
##            
##            #1-norms of blocks entries of A
##            Block_A = Bt * csr_matrix((abs(A.data),A.indices,A.indptr),dims=A.shape) * B 
##    
##            S = sa_strong_connections(Block_A,epsilon)       
##            S.data[:] = 1
##
##            Mask = B * S * Bt
##
##            A_strong = A ** Mask
##            #A_weak   = A - A_strong
##            A_filtered = A_strong

    return A_filtered          


def sa_strong_connections(A,epsilon):
    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')
    
    Sp,Sj,Sx = multigridtools.sa_strong_connections(A.shape[0],epsilon,A.indptr,A.indices,A.data)

    return csr_matrix((Sx,Sj,Sp),A.shape)


def sa_constant_interpolation(A,epsilon,blocks=None):
    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')
    
    if blocks is not None:
        num_dofs   = A.shape[0]
        num_blocks = blocks.max() + 1
        
        if num_dofs != len(blocks):
            raise ValueError,'improper block specification'
       
        #print "SA has blocks"
        # for non-scalar problems, use pre-defined blocks in aggregation
        # the strength of connection matrix is based on the Frobenius norms of the blocks
       
        B  = csr_matrix((ones(num_dofs),blocks,arange(num_dofs + 1)),dims=(num_dofs,num_blocks))
        #1-norms of blocks entries of A
        Block_A = B.T.tocsr() * csr_matrix((abs(A.data),A.indices,A.indptr),dims=A.shape) * B 

        S = sa_strong_connections(Block_A,epsilon)
    
        Pj = multigridtools.sa_get_aggregates(S.shape[0],S.indptr,S.indices)
        Pj = Pj[blocks] #expand block aggregates into constituent dofs
        Pp = B.indptr
        Px = B.data
    else:
        S = sa_strong_connections(A,epsilon)
        
        Pj = multigridtools.sa_get_aggregates(S.shape[0],S.indptr,S.indices)
        Pp = numpy.arange(len(Pj)+1)
        Px = numpy.ones(len(Pj))
    
    return csr_matrix((Px,Pj,Pp))


def sa_fit_candidates(AggOp,candidates):
    K = len(candidates)

    N_fine,N_coarse = AggOp.shape

    if K > 1 and len(candidates[0]) == K*N_fine:
        #see if fine space has been expanded (all levels except for first)
        AggOp = csr_matrix((AggOp.data.repeat(K),AggOp.indices.repeat(K),arange(K*N_fine + 1)),dims=(K*N_fine,N_coarse))
        N_fine = K*N_fine

    #TODO convert this to list of coarse candidates
    R = zeros((K*N_coarse,K)) #storage for coarse candidates

    candidate_matrices = []
    for i,c in enumerate(candidates):
        c = c[diff(AggOp.indptr) == 1]  #eliminate DOFs that aggregation misses
        X = csr_matrix((c,AggOp.indices,AggOp.indptr),dims=AggOp.shape)
       
        #orthogonalize X against previous
        for j,A in enumerate(candidate_matrices):
            D_AtX = csr_matrix((A.data*X.data,X.indices,X.indptr),dims=X.shape).sum(axis=0).A.flatten() #same as diagonal of A.T * X            
            R[j::K,i] = D_AtX
            X.data -= D_AtX[X.indices] * A.data
         
        #normalize X
        D_XtX = csr_matrix((X.data**2,X.indices,X.indptr),dims=X.shape).sum(axis=0).A.flatten() #same as diagonal of X.T * X            
        col_norms = sqrt(D_XtX)
        R[i::K,i] = col_norms
        col_norms = 1.0/col_norms
        col_norms[isinf(col_norms)] = 0
        X.data *= col_norms[X.indices]

        candidate_matrices.append(X)

    Q_indptr  = K*AggOp.indptr
    Q_indices = (K*AggOp.indices).repeat(K)
    for i in range(K):
        Q_indices[i::K] += i
    Q_data = empty(AggOp.indptr[-1] * K) #if AggOp includes all nodes, then this is (N_fine * K)
    for i,X in enumerate(candidate_matrices):
        Q_data[i::K] = X.data
    Q = csr_matrix((Q_data,Q_indices,Q_indptr),dims=(N_fine,K*N_coarse))

    coarse_candidates = [ascontiguousarray(R[:,i]) for i in range(K)]

    return Q,coarse_candidates
    
    
def sa_interpolation(A,candidates,epsilon=0.0,omega=4.0/3.0,blocks=None,AggOp=None):
    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')

    if AggOp is None:
        AggOp = sa_constant_interpolation(A,epsilon=epsilon,blocks=blocks)
    else:
        if not isspmatrix_csr(AggOp):
            raise TypeError,'aggregation is specified by a list of csr_matrix objects'
        if A.shape[1] != AggOp.shape[0]:
            raise ValueError,'incompatible aggregation operator'

    T,coarse_candidates = sa_fit_candidates(AggOp,candidates)  

    A_filtered = sa_filtered_matrix(A,epsilon,blocks) #use filtered matrix for anisotropic problems

    D_inv    = diag_sparse(1.0/diag_sparse(A_filtered))       
    D_inv_A  = D_inv * A_filtered
    D_inv_A *= omega/approximate_spectral_radius(D_inv_A)

    # smooth tentative prolongator T
    P = T - (D_inv_A*T)
           
    return P,coarse_candidates



