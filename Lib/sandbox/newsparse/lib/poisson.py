import spmatrix

def poisson1d(n):
    L = spmatrix.ll_mat(n, n)
    for i in range(n):
        L[i,i] = 2
        if i > 0:
            L[i,i-1] = -1
        if i < n-1:
            L[i,i+1] = -1
    return L
 
def poisson1d_sym(n):
    L = spmatrix.ll_mat_sym(n)
    for i in range(n):
        L[i,i] = 2
        if i > 0:
            L[i,i-1] = -1
    return L
        
def poisson2d(n):
    L = spmatrix.ll_mat(n*n, n*n)
    for i in range(n):
        for j in range(n):
            k = i + n*j
            L[k,k] = 4
            if i > 0:
                L[k,k-1] = -1
            if i < n-1:
                L[k,k+1] = -1
            if j > 0:
                L[k,k-n] = -1
            if j < n-1:
                L[k,k+n] = -1
    return L

def poisson2d_sym(n):
    L = spmatrix.ll_mat_sym(n*n)
    for i in range(n):
        for j in range(n):
            k = i + n*j
            L[k,k] = 4
            if i > 0:
                L[k,k-1] = -1
            if j > 0:
                L[k,k-n] = -1
    return L

def poisson2d_sym_blk(n):
    L = spmatrix.ll_mat_sym(n*n)
    I = spmatrix.ll_mat_sym(n)
    for i in range(n):
        I[i,i] = -1
    P = spmatrix.ll_mat_sym(n)
    for i in range(n):
        P[i,i] = 4
        if i > 0:
            P[i,i-1] = -1
    for i in range(0, n*n, n):
        L[i:i+n,i:i+n] = P
        if i > 0: L[i:i+n,i-n:i] = I
    return L
