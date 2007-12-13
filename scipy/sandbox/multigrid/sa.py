import scipy
import numpy
from numpy import array,arange,ones,zeros,sqrt,isinf,asarray,empty,diff,\
                  ascontiguousarray
from scipy.sparse import csr_matrix,isspmatrix_csr

from utils import diag_sparse, approximate_spectral_radius, \
                  symmetric_rescaling, expand_into_blocks
import multigridtools

__all__ = ['sa_filtered_matrix','sa_strong_connections','sa_constant_interpolation',
           'sa_interpolation','sa_smoothed_prolongator','sa_fit_candidates']



def sa_filtered_matrix(A,epsilon,blocks=None):
    """The filtered matrix is obtained from A by lumping all weak off-diagonal
    entries onto the diagonal.  Weak off-diagonals are determined by
    the standard strength of connection measure using the parameter epsilon.

    In the case epsilon = 0.0, (i.e. no weak connections) A is returned.
    """

    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')

    if epsilon == 0:
        A_filtered = A
    else:
        if blocks is None:
            Sp,Sj,Sx = multigridtools.sa_strong_connections(A.shape[0],epsilon,A.indptr,A.indices,A.data)
            A_filtered = csr_matrix((Sx,Sj,Sp),A.shape)
        else:
            raise NotImplementedError,'blocks not handled yet'
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
        #TODO figure out what to do for blocks here
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


def sa_fit_candidates(AggOp,candidates,tol=1e-10):
    #TODO handle non-floating point candidates better
    candidates = candidates.astype('float64')

    K = candidates.shape[1] # num candidates

    N_fine,N_coarse = AggOp.shape

    if K > 1 and candidates.shape[0] == K*N_fine:
        #see if fine space has been expanded (all levels except for first)
        AggOp = expand_into_blocks(AggOp,K,1).tocsr()
        N_fine = K*N_fine

    R = zeros((N_coarse,K,K)) #storage for coarse candidates

    candidate_matrices = []

    for i in range(K):
        c = candidates[:,i]
        c = c[diff(AggOp.indptr) == 1]     # eliminate DOFs that aggregation misses

        threshold = tol * abs(c).max()   # cutoff for small basis functions

        X = csr_matrix((c,AggOp.indices,AggOp.indptr),dims=AggOp.shape)

        #orthogonalize X against previous
        for j,A in enumerate(candidate_matrices):
            D_AtX = csr_matrix((A.data*X.data,X.indices,X.indptr),dims=X.shape).sum(axis=0).A.flatten() #same as diagonal of A.T * X
            R[:,j,i] = D_AtX
            X.data -= D_AtX[X.indices] * A.data

        #normalize X
        D_XtX = csr_matrix((X.data**2,X.indices,X.indptr),dims=X.shape).sum(axis=0).A.flatten() #same as diagonal of X.T * X
        col_norms = sqrt(D_XtX)
        mask = col_norms < threshold   # set small basis functions to 0
        col_norms[mask] = 0
        R[:,i,i] = col_norms
        col_norms = 1.0/col_norms
        col_norms[mask] = 0
        X.data *= col_norms[X.indices]

        candidate_matrices.append(X)

    # expand AggOp blocks horizontally
    Q_indptr  = K*AggOp.indptr
    Q_indices = (K*AggOp.indices).repeat(K)
    for i in range(K):
        Q_indices[i::K] += i
    Q_data = empty(AggOp.indptr[-1] * K) #if AggOp includes all nodes, then this is (N_fine * K)
    for i,X in enumerate(candidate_matrices):
        Q_data[i::K] = X.data
    Q = csr_matrix((Q_data,Q_indices,Q_indptr),dims=(N_fine,K*N_coarse))

    R = R.reshape(-1,K)

    return Q,R

def sa_smoothed_prolongator(A,T,epsilon,omega,blocks=None):
    """For a given matrix A and tentative prolongator T return the
    smoothed prolongator P

        P = (I - omega/rho(S) S) * T

    where S is a Jacobi smoothing operator defined as follows:

        omega      - damping parameter
        rho(S)     - spectral radius of S (estimated)
        S          - inv(diag(A_filtered)) * A_filtered   (Jacobi smoother)
        A_filtered - sa_filtered_matrix(A,epsilon)
    """


    A_filtered = sa_filtered_matrix(A,epsilon,blocks) #use filtered matrix for anisotropic problems

    D_inv    = diag_sparse(1.0/diag_sparse(A_filtered))
    D_inv_A  = D_inv * A_filtered
    D_inv_A *= omega/approximate_spectral_radius(D_inv_A)

    # smooth tentative prolongator T
    P = T - (D_inv_A*T)

    return P

def sa_interpolation(A,candidates,epsilon=0.0,omega=4.0/3.0,blocks=None,AggOp=None):
    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')

    if AggOp is None:
        AggOp = sa_constant_interpolation(A,epsilon=epsilon,blocks=blocks)
    else:
        if not isspmatrix_csr(AggOp):
            raise TypeError,'expected csr_matrix for argument AggOp'
        if A.shape[1] != AggOp.shape[0]:
            raise ValueError,'incompatible aggregation operator'


    T,coarse_candidates = sa_fit_candidates(AggOp,candidates)
    #T = AggOp #TODO test

    P = sa_smoothed_prolongator(A,T,epsilon,omega,blocks)

    if blocks is not None:
        blocks = arange(AggOp.shape[1]).repeat(candidates.shape[1])

    return P,coarse_candidates,blocks
