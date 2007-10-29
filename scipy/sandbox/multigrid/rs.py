from scipy.sparse import csr_matrix,isspmatrix_csr

import multigridtools

__all__ = ['rs_strong_connections','rs_interpolation']

def rs_strong_connections(A,theta):
    """Return a strength of connection matrix using the method of Ruge and Stuben

        An off-diagonal entry A[i.j] is a strong connection iff

                -A[i,j] >= theta * max( -A[i,k] )   where k != i
    """
    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')

    Sp,Sj,Sx = multigridtools.rs_strong_connections(A.shape[0],theta,A.indptr,A.indices,A.data)
    return csr_matrix((Sx,Sj,Sp),dims=A.shape)


def rs_interpolation(A,theta=0.25):
    if not isspmatrix_csr(A): raise TypeError('expected csr_matrix')

    S = rs_strong_connections(A,theta)

    T = S.T.tocsr()  #transpose S for efficient column access

    Ip,Ij,Ix = multigridtools.rs_interpolation(A.shape[0],\
                                               A.indptr,A.indices,A.data,\
                                               S.indptr,S.indices,S.data,\
                                               T.indptr,T.indices,T.data)

    return csr_matrix((Ix,Ij,Ip))
