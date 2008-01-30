from numpy import array, arange, ones, zeros, sqrt, asarray, empty, diff
from scipy.sparse import csr_matrix, isspmatrix_csr, bsr_matrix, isspmatrix_bsr

import multigridtools
from multilevel import multilevel_solver
from utils import diag_sparse, approximate_spectral_radius, \
                  symmetric_rescaling, scale_columns

__all__ = ['smoothed_aggregation_solver',
        'sa_filtered_matrix','sa_strong_connections','sa_standard_aggregation',
        'sa_smoothed_prolongator','sa_fit_candidates']



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

def sa_strong_connections(A,epsilon=0):
    """Compute a strength of connection matrix C

        C[i,j] = 1 if abs(A[i,j]) >= epsilon * abs(A[i,i] * A[j,j])
        C[i,j] = 0 otherwise

    """

    if isspmatrix_csr(A):
        if epsilon == 0:
            return A
        else:
            fn = multigridtools.sa_strong_connections
            Sp,Sj,Sx = fn(A.shape[0],epsilon,A.indptr,A.indices,A.data)
            return csr_matrix((Sx,Sj,Sp),A.shape)

    elif isspmatrix_bsr(A):
        M,N = A.shape
        R,C = A.blocksize

        if R != C:
            raise ValueError,'matrix must have square blocks'

        if epsilon == 0:
            data = ones( len(A.indices), dtype=A.dtype )
            return csr_matrix((data,A.indices,A.indptr),shape=(M/R,N/C))
        else:
            # the strength of connection matrix is based on the 
            # Frobenius norms of the blocks
            data = (A.data*A.data).reshape(-1,R*C).sum(axis=1) 
            A = csr_matrix((data,A.indices,A.indptr),shape=(M/R,N/C))
            return sa_strong_connections(A,epsilon)
    else:
        raise TypeError('expected csr_matrix or bsr_matrix') 

def sa_standard_aggregation(C):
    """Compute the sparsity pattern of the tentative prolongator 
    from a strength of connection matrix C
    """
    if isspmatrix_csr(C): 
        
        Pj = multigridtools.sa_get_aggregates(C.shape[0],C.indptr,C.indices)
        Pp = arange(len(Pj)+1)
        Px = ones(len(Pj)) #TODO replace this with something else?
        
        return csr_matrix((Px,Pj,Pp))
    else:
        raise TypeError('expected csr_matrix') 

def sa_fit_candidates(AggOp,B,tol=1e-10):
    if not isspmatrix_csr(AggOp):
        raise TypeError,'expected csr_matrix for argument AggOp'

    if B.dtype != 'float32':
        B = asarray(B,dtype='float64')

    if len(B.shape) != 2:
        raise ValueError,'expected rank 2 array for argument B'

    if B.shape[0] % AggOp.shape[0] != 0:
        raise ValueError,'dimensions of AggOp %s and B %s are incompatible' % (AggOp.shape, B.shape)
    

    K = B.shape[1] # number of near-nullspace candidates
    blocksize = B.shape[0] / AggOp.shape[0]

    N_fine,N_coarse = AggOp.shape

    R = zeros((N_coarse,K,K), dtype=B.dtype) #storage for coarse candidates

    candidate_matrices = []

    for i in range(K):
        c = B[:,i]
        c = c.reshape(-1,blocksize,1)[diff(AggOp.indptr) == 1]     # eliminate DOFs that aggregation misses

        X = bsr_matrix( (c, AggOp.indices, AggOp.indptr), \
                shape=(blocksize*N_fine, N_coarse) )

        col_thresholds = tol * bsr_matrix((X.data**2,X.indices,X.indptr),shape=X.shape).sum(axis=0).A.flatten() 

        #orthogonalize X against previous
        for j,A in enumerate(candidate_matrices):
            D_AtX = bsr_matrix((A.data*X.data,X.indices,X.indptr),shape=X.shape).sum(axis=0).A.flatten() #same as diagonal of A.T * X
            R[:,j,i] = D_AtX
            X.data -= scale_columns(A,D_AtX).data

        #normalize X
        col_norms = bsr_matrix((X.data**2,X.indices,X.indptr),shape=X.shape).sum(axis=0).A.flatten() #same as diagonal of X.T * X
        mask = col_norms < col_thresholds   # set small basis functions to 0

        col_norms = sqrt(col_norms)
        col_norms[mask] = 0
        R[:,i,i] = col_norms
        col_norms = 1.0/col_norms
        col_norms[mask] = 0
        scale_columns(X,col_norms,copy=False)

        candidate_matrices.append(X)

    Q_indptr  = AggOp.indptr
    Q_indices = AggOp.indices
    Q_data = empty((AggOp.nnz,blocksize,K)) #if AggOp includes all nodes, then this is (N_fine * K)
    for i,X in enumerate(candidate_matrices):
        Q_data[:,:,i] = X.data.reshape(-1,blocksize)
    Q = bsr_matrix((Q_data,Q_indices,Q_indptr),shape=(blocksize*N_fine,K*N_coarse))

    R = R.reshape(-1,K)

    return Q,R

def sa_smoothed_prolongator(A,T,epsilon=0.0,omega=4.0/3.0):
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

    # TODO use scale_rows()
    D = A_filtered.diagonal()
    D_inv = 1.0 / D
    D_inv[D == 0] = 0

    D_inv_A  = diag_sparse(D_inv) * A_filtered
    D_inv_A *= omega/approximate_spectral_radius(D_inv_A)

    # smooth tentative prolongator T
    P = T - (D_inv_A*T)

    return P



def smoothed_aggregation_solver(A, B=None, 
        max_levels = 10, 
        max_coarse = 500,
        strength   = sa_strong_connections, 
        aggregate  = sa_standard_aggregation,
        tentative  = sa_fit_candidates,
        smooth     = sa_smoothed_prolongator):
    """Create a multilevel solver using Smoothed Aggregation (SA)

    *Parameters*:

        A : {csr_matrix, bsr_matrix}
            Square matrix in CSR or BSR format
        B : {None, array_like}
            Near-nullspace candidates stored in the columns of an NxK array.
            The default value B=None is equivalent to B=ones((N,1))
        max_levels: {integer} : default 10
            Maximum number of levels to be used in the multilevel solver.
        max_coarse: {integer} : default 500
            Maximum number of variables permitted on the coarse grid.
        strength :
            Function that computes the strength of connection matrix C
                strength(A) -> C
        aggregate : 
            Function that computes an aggregation operator
                aggregate(C) -> AggOp
        tentative:
            Function that computes a tentative prolongator
                tentative(AggOp,B) -> T,B_coarse
        smooth :
            Function that smooths the tentative prolongator
                smooth(A,C,T) -> P

    Unused Parameters
        epsilon: {float} : default 0.0
            Strength of connection parameter used in aggregation.
        omega: {float} : default 4.0/3.0
            Damping parameter used in prolongator smoothing (0 < omega < 2)
        symmetric: {boolean} : default True
            True if A is symmetric, False otherwise
        rescale: {boolean} : default True
            If True, symmetrically rescale A by the diagonal
            i.e. A -> D * A * D,  where D is diag(A)^-0.5
        aggregation: {None, list of csr_matrix} : optional
            List of csr_matrix objects that describe a user-defined
            multilevel aggregation of the variables.
            TODO ELABORATE

    *Example*:
        TODO

    *References*:
        "Algebraic Multigrid by Smoothed Aggregation for Second and Fourth Order Elliptic Problems",
            Petr Vanek and Jan Mandel and Marian Brezina
            http://citeseer.ist.psu.edu/vanek96algebraic.html

    """

    A = A.asfptype()

    if not (isspmatrix_csr(A) or isspmatrix_bsr(A)):
        raise TypeError('argument A must have type csr_matrix or bsr_matrix')

    if A.shape[0] != A.shape[1]:
        raise ValueError('expected square matrix')

    if B is None:
        B = ones((A.shape[0],1),dtype=A.dtype) # use constant vector
    else:
        B = asarray(B,dtype=A.dtype)

    pre,post = None,None   #preprocess/postprocess

    #if rescale:
    #    D_sqrt,D_sqrt_inv,A = symmetric_rescaling(A)
    #    D_sqrt,D_sqrt_inv = diag_sparse(D_sqrt),diag_sparse(D_sqrt_inv)

    #    B = D_sqrt * B  #scale candidates
    #    def pre(x,b):
    #        return D_sqrt*x,D_sqrt_inv*b
    #    def post(x):
    #        return D_sqrt_inv*x

    As = [A]
    Ps = []
    Rs = []

    while len(As) < max_levels and A.shape[0] > max_coarse:
        C     = strength(A)
        AggOp = aggregate(C)
        T,B   = tentative(AggOp,B)
        P     = smooth(A,T)

        R = P.T.asformat(P.format)

        A = R * A * P     #galerkin operator

        As.append(A)
        Rs.append(R)
        Ps.append(P)


    return multilevel_solver(As,Ps,Rs=Rs,preprocess=pre,postprocess=post)


