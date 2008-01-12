import scipy
import numpy
from numpy import array,arange,ones,zeros,sqrt,isinf,asarray,empty,diff,\
                  ascontiguousarray
from scipy.sparse import csr_matrix, isspmatrix_csr, bsr_matrix, isspmatrix_bsr

from utils import diag_sparse, approximate_spectral_radius, \
                  symmetric_rescaling, expand_into_blocks, scale_columns
import multigridtools

__all__ = ['sa_filtered_matrix','sa_strong_connections','sa_constant_interpolation',
           'sa_interpolation','sa_smoothed_prolongator','sa_fit_candidates']



def sa_filtered_matrix(A,epsilon):
    """The filtered matrix is obtained from A by lumping all weak off-diagonal
    entries onto the diagonal.  Weak off-diagonals are determined by
    the standard strength of connection measure using the parameter epsilon.

    In the case epsilon = 0.0, (i.e. no weak connections) A is returned.
    """

    if epsilon == 0:
        return A

    if isspmatrix_csr(A): 
        Sp,Sj,Sx = multigridtools.sa_strong_connections(A.shape[0],epsilon,A.indptr,A.indices,A.data)
        return csr_matrix((Sx,Sj,Sp),shape=A.shape)
    elif ispmatrix_bsr(A):
        raise NotImplementedError,'blocks not handled yet'
    else:
        return sa_filtered_matrix(csr_matrix(A),epsilon)
##            #TODO subtract weak blocks from diagonal blocks?
##            num_dofs   = A.shape[0]
##            num_blocks = blocks.max() + 1
##
##            if num_dofs != len(blocks):
##                raise ValueError,'improper block specification'
##
##            # for non-scalar problems, use pre-defined blocks in aggregation
##            # the strength of connection matrix is based on the 1-norms of the blocks
##
##            B  = csr_matrix((ones(num_dofs),blocks,arange(num_dofs + 1)),shape=(num_dofs,num_blocks))
##            Bt = B.T.tocsr()
##
##            #1-norms of blocks entries of A
##            Block_A = Bt * csr_matrix((abs(A.data),A.indices,A.indptr),shape=A.shape) * B
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


def sa_constant_interpolation(A,epsilon):
    """Compute the sparsity pattern of the tentative prolongator
    """
    if isspmatrix_csr(A): 
        S = sa_strong_connections(A,epsilon)

        Pj = multigridtools.sa_get_aggregates(S.shape[0],S.indptr,S.indices)
        Pp = numpy.arange(len(Pj)+1)
        Px = numpy.ones(len(Pj)) #TODO replace this with something else?
        
        return csr_matrix((Px,Pj,Pp))

    elif isspmatrix_bsr(A):

        # the strength of connection matrix is based on the Frobenius norms of the blocks
        M,N = A.shape
        R,C = A.blocksize

        if R != C:
            raise ValueError,'matrix must have square blocks'

        f_norms = (A.data*A.data).reshape(-1,R*C).sum(axis=1) #Frobenius norm of each block
        
        A = csr_matrix((f_norms,A.indices,A.indptr),shape=(M/R,N/C))

        return sa_constant_interpolation(A,epsilon)

    else:
        sa_constant_interpolation(csr_matrix(A),epsilon)

def sa_fit_candidates(AggOp,candidates,tol=1e-10):
    if candidates.dtype != 'float32':
        candidates = asarray(candidates,dtype='float64')

    K = candidates.shape[1] # number of near-nullspace candidates
    blocksize = candidates.shape[0] / AggOp.shape[0]

    N_fine,N_coarse = AggOp.shape

    #if blocksize > 1:
    #    #see if fine space has been expanded (all levels except for first)
    #    AggOp = expand_into_blocks(AggOp,blocksize,1).tocsr()

    R = zeros((N_coarse,K,K),dtype=candidates.dtype) #storage for coarse candidates

    candidate_matrices = []

    threshold = tol * abs(candidates).max()   # cutoff for small basis functions

    for i in range(K):
        c = candidates[:,i]
        c = c.reshape(-1,blocksize,1)[diff(AggOp.indptr) == 1]     # eliminate DOFs that aggregation misses

        X = bsr_matrix( (c, AggOp.indices, AggOp.indptr), \
                shape=(blocksize*N_fine, N_coarse) )

        #orthogonalize X against previous
        for j,A in enumerate(candidate_matrices):
            #import pdb; pdb.set_trace()
            D_AtX = bsr_matrix((A.data*X.data,X.indices,X.indptr),shape=X.shape).sum(axis=0).A.flatten() #same as diagonal of A.T * X
            R[:,j,i] = D_AtX
            X.data -= scale_columns(A,D_AtX).data

        #normalize X
        D_XtX = bsr_matrix((X.data**2,X.indices,X.indptr),shape=X.shape).sum(axis=0).A.flatten() #same as diagonal of X.T * X
        col_norms = sqrt(D_XtX)
        mask = col_norms < threshold   # set small basis functions to 0
        col_norms[mask] = 0
        R[:,i,i] = col_norms
        col_norms = 1.0/col_norms
        col_norms[mask] = 0
        scale_columns(X,col_norms,copy=False)

        candidate_matrices.append(X)

    # expand AggOp blocks horizontally
    Q_indptr  = AggOp.indptr
    Q_indices = AggOp.indices
    Q_data = empty((AggOp.nnz,blocksize,K)) #if AggOp includes all nodes, then this is (N_fine * K)
    for i,X in enumerate(candidate_matrices):
        Q_data[:,:,i] = X.data.reshape(-1,blocksize)
    Q = bsr_matrix((Q_data,Q_indices,Q_indptr),shape=(blocksize*N_fine,K*N_coarse))

    R = R.reshape(-1,K)

    return Q,R

def sa_smoothed_prolongator(A,T,epsilon,omega):
    """For a given matrix A and tentative prolongator T return the
    smoothed prolongator P

        P = (I - omega/rho(S) S) * T

    where S is a Jacobi smoothing operator defined as follows:

        omega      - damping parameter
        rho(S)     - spectral radius of S (estimated)
        S          - inv(diag(A_filtered)) * A_filtered   (Jacobi smoother)
        A_filtered - sa_filtered_matrix(A,epsilon)
    """


    A_filtered = sa_filtered_matrix(A,epsilon) #use filtered matrix for anisotropic problems

    D_inv    = diag_sparse(1.0/diag_sparse(A_filtered))
    D_inv_A  = D_inv * A_filtered
    D_inv_A *= omega/approximate_spectral_radius(D_inv_A)

    # smooth tentative prolongator T
    P = T - (D_inv_A*T)

    return P

def sa_interpolation(A,candidates,epsilon=0.0,omega=4.0/3.0,AggOp=None):

    if not (isspmatrix_csr(A) or isspmatrix_bsr(A)):
        A = csr_matrix(A)

    if AggOp is None:
        AggOp = sa_constant_interpolation(A,epsilon=epsilon)
    else:
        if not isspmatrix_csr(AggOp):
            AggOp = csr_matrix(AggOp)
        if A.shape[1] != AggOp.shape[0]:
            raise ValueError,'incompatible aggregation operator'


    T,coarse_candidates = sa_fit_candidates(AggOp,candidates)
    P = sa_smoothed_prolongator(A,T,epsilon,omega)
    return P,coarse_candidates

