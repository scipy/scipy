#pythran export _Aij(float[:,:], int, int)
#pythran export _Aij(int[:,:], int, int)
def _Aij(A, i, j):
    """Sum of upper-left and lower right blocks of contingency table."""
    # See `somersd` References [2] bottom of page 309
    return A[:i, :j].sum() + A[i+1:, j+1:].sum()

#pythran export _Dij(float[:,:], int, int)
#pythran export _Dij(int[:,:], int, int)
def _Dij(A, i, j):
    """Sum of lower-left and upper-right blocks of contingency table."""
    # See `somersd` References [2] bottom of page 309
    return A[i+1:, :j].sum() + A[:i, j+1:].sum()

#pythran export _P(float[:,:])
#pythran export _P(int[:,:])
def _P(A):
    """Twice the number of concordant pairs, excluding ties."""
    # See `somersd` References [2] bottom of page 309
    m, n = A.shape
    count = 0
    for i in range(m):
        for j in range(n):
            count += A[i, j]*_Aij(A, i, j)
    return count

#pythran export _Q(float[:,:])
#pythran export _Q(int[:,:])
def _Q(A):
    """Twice the number of discordant pairs, excluding ties."""
    # See `somersd` References [2] bottom of page 309
    m, n = A.shape
    count = 0
    for i in range(m):
        for j in range(n):
            count += A[i, j]*_Dij(A, i, j)
    return count

#pythran export _a_ij_Aij_Dij2(float[:,:])
#pythran export _a_ij_Aij_Dij2(int[:,:])
def _a_ij_Aij_Dij2(A):
    """A term that appears in the ASE of Kendall's tau and Somers' D."""
    # See `somersd` References [2] section 4: Modified ASEs to test the null hypothesis...
    m, n = A.shape
    count = 0
    for i in range(m):
        for j in range(n):
            count += A[i, j]*(_Aij(A, i, j) - _Dij(A, i, j))**2
    return count
