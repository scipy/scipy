#pythran export _Aij_pythran(float[:,:], int, int)
#pythran export _Aij_pythran(int[:,:], int, int)
def _Aij_pythran(A, i, j):
    """Sum of upper-left and lower right blocks of contingency table."""
    # See [2] bottom of page 309
    return A[:i, :j].sum() + A[i+1:, j+1:].sum()

#pythran export _Dij_pythran(float[:,:], int, int)
#pythran export _Dij_pythran(int[:,:], int, int)
def _Dij_pythran(A, i, j):
    """Sum of lower-left and upper-right blocks of contingency table."""
    # See [2] bottom of page 309
    return A[i+1:, :j].sum() + A[:i, j+1:].sum()

#pythran export _P_pythran(float[:,:])
#pythran export _P_pythran(int[:,:])
def _P_pythran(A):
    """Twice the number of concordant pairs, excluding ties."""
    # See [2] bottom of page 309
    m, n = A.shape
    count = 0
    for i in range(m):
        for j in range(n):
            count += A[i, j]*_Aij_pythran(A, i, j)
    return count

#pythran export _Q_pythran(float[:,:])
#pythran export _Q_pythran(int[:,:])
def _Q_pythran(A):
    """Twice the number of discordant pairs, excluding ties."""
    # See [2] bottom of page 309
    m, n = A.shape
    count = 0
    for i in range(m):
        for j in range(n):
            count += A[i, j]*_Dij_pythran(A, i, j)
    return count

#pythran export _a_ij_Aij_Dij2_pythran(float[:,:])
#pythran export _a_ij_Aij_Dij2_pythran(int[:,:])
def _a_ij_Aij_Dij2_pythran(A):
    """A term that appears in the ASE of Kendall's tau and Somers' D."""
    # See [2] section 4: Modified ASEs to test the null hypothesis...
    m, n = A.shape
    count = 0
    for i in range(m):
        for j in range(n):
            count += A[i, j]*(_Aij_pythran(A, i, j) - _Dij_pythran(A, i, j))**2
    return count