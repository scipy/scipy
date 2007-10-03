##import scipy
##import numpy
##
##from numpy import arange,ones,zeros,sqrt,isinf,asarray,empty
##from scipy.sparse import csr_matrix,isspmatrix_csr
##
##from utils import diag_sparse,approximate_spectral_radius
##import multigridtools
##
##
##def rs_strong_connections(A,theta):
##    """
##    Return a strength of connection matrix using the method of Ruge and Stuben
##
##        An off-diagonal entry A[i.j] is a strong connection iff
##
##                -A[i,j] >= theta * max( -A[i,k] )   where k != i
##    """
##    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')
##
##    Sp,Sj,Sx = multigridtools.rs_strong_connections(A.shape[0],theta,A.indptr,A.indices,A.data)
##    return csr_matrix((Sx,Sj,Sp),A.shape)
##
##
##def rs_interpolation(A,theta=0.25):
##    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')
##    
##    S = rs_strong_connections(A,theta)
##
##    T = S.T.tocsr()  #transpose S for efficient column access
##
##    Ip,Ij,Ix = multigridtools.rs_interpolation(A.shape[0],\
##                                               A.indptr,A.indices,A.data,\
##                                               S.indptr,S.indices,S.data,\
##                                               T.indptr,T.indices,T.data)
##
##    return csr_matrix((Ix,Ij,Ip))
##
##
##def sa_strong_connections(A,epsilon):
##    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')
##
##    Sp,Sj,Sx = multigridtools.sa_strong_connections(A.shape[0],epsilon,A.indptr,A.indices,A.data)
##
##    return csr_matrix((Sx,Sj,Sp),A.shape)
##
##def sa_constant_interpolation(A,epsilon,blocks=None):
##    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')
##    
##    #handle epsilon = 0 case without creating strength of connection matrix?
##    
##    if blocks is not None:
##        num_dofs   = A.shape[0]
##        num_blocks = blocks.max()
##        
##        if num_dofs != len(blocks):
##            raise ValueError,'improper block specification'
##        
##        # for non-scalar problems, use pre-defined blocks in aggregation
##        # the strength of connection matrix is based on the Frobenius norms of the blocks
##        
##        B  = csr_matrix((ones(num_dofs),blocks,arange(num_dofs + 1)),dims=(num_dofs,num_blocks))
##        Block_Frob = B.T.tocsr() * csr_matrix((A.data**2,A.indices,A.indptr),dims=A.shape) * B #Frobenius norms of blocks entries of A
##
##        S = sa_strong_connections(Block_Frob,epsilon)
##    
##        Pj = multigridtools.sa_get_aggregates(S.shape[0],S.indptr,S.indices)
##        Pj = Pj[blocks] #expand block aggregates into constituent dofs
##        Pp = B.indptr
##        Px = B.data
##    else:
##        S = sa_strong_connections(A,epsilon)
##        
##        Pj = multigridtools.sa_get_aggregates(S.shape[0],S.indptr,S.indices)
##        Pp = numpy.arange(len(Pj)+1)
##        Px = numpy.ones(len(Pj))
##    
##    return csr_matrix((Px,Pj,Pp))
##
##
##def fit_candidates(AggOp,candidates):
##    K = len(candidates)
##
##    N_fine,N_coarse = AggOp.shape
##
##    if K > 1 and len(candidates[0]) == K*N_fine:
##        #see if fine space has been expanded (all levels except for first)
##        AggOp = csr_matrix((AggOp.data.repeat(K),AggOp.indices.repeat(K),arange(K*N_fine + 1)),dims=(K*N_fine,N_coarse))
##        N_fine = K*N_fine
##
##    R = zeros((K*N_coarse,K))
##
##    candidate_matrices = []
##    for i,c in enumerate(candidates):
##        X = csr_matrix((c.copy(),AggOp.indices,AggOp.indptr),dims=AggOp.shape)
##       
##        #TODO optimize this  
##
##        #orthogonalize X against previous
##        for j,A in enumerate(candidate_matrices):
##            D_AtX = csr_matrix((A.data*X.data,X.indices,X.indptr),dims=X.shape).sum(axis=0).A.flatten() #same as diagonal of A.T * X            
##            R[j::K,i] = D_AtX
##            X.data -= D_AtX[X.indices] * A.data
##
##            #AtX = csr_matrix(A.T.tocsr() * X
##            #R[j::K,i] = AtX.data
##            #X = X - A * AtX 
##    
##        #normalize X
##        XtX = X.T.tocsr() * X
##        col_norms = sqrt(asarray(XtX.sum(axis=0)).flatten())
##        R[i::K,i] = col_norms
##        col_norms = 1.0/col_norms
##        col_norms[isinf(col_norms)] = 0
##        X.data *= col_norms[X.indices]
##
##        candidate_matrices.append(X)
##
##
##    Q_indptr  = K*AggOp.indptr
##    Q_indices = (K*AggOp.indices).repeat(K)
##    for i in range(K):
##        Q_indices[i::K] += i
##    Q_data = empty(N_fine * K)
##    for i,X in enumerate(candidate_matrices):
##        Q_data[i::K] = X.data
##    Q = csr_matrix((Q_data,Q_indices,Q_indptr),dims=(N_fine,K*N_coarse))
##
##    coarse_candidates = [R[:,i] for i in range(K)]
##
##    return Q,coarse_candidates
##    
####    S = sa_strong_connections(A,epsilon)
####
####    #tentative (non-smooth) interpolation operator I
####    Pj = multigridtools.sa_get_aggregates(S.shape[0],S.indptr,S.indices)
####    Pp = numpy.arange(len(Pj)+1)
####    Px = numpy.ones(len(Pj))
####    
####    return scipy.sparse.csr_matrix((Px,Pj,Pp))
##
####def sa_smoother(A,S,omega):
####    Bp,Bj,Bx = multigridtools.sa_smoother(A.shape[0],omega,A.indptr,A.indices,A.data,S.indptr,S.indices,S.data)
####
####    return csr_matrix((Bx,Bj,Bp),dims=A.shape)
##    
##def sa_interpolation(A,candidates,epsilon,omega=4.0/3.0,blocks=None):
##    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')
##   
##    AggOp  = sa_constant_interpolation(A,epsilon=epsilon,blocks=blocks)
##    T,coarse_candidates = fit_candidates(AggOp,candidates)
##
##    D_inv = diag_sparse(1.0/diag_sparse(A))       
##    
##    D_inv_A  = D_inv * A
##    D_inv_A *= omega/approximate_spectral_radius(D_inv_A)
##
##    P = T - (D_inv_A*T)  #same as I=S*P, (faster?)
##           
##    return P,coarse_candidates
##
##
##
