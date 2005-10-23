import Numeric
import math
import spmatrix
import itsolvers
import precon
import jdsym
import time

def poisson2d_sym_blk(n):
    L = spmatrix.ll_mat_sym(n*n)
    I = spmatrix.ll_mat_sym(n)
    P = spmatrix.ll_mat_sym(n)
    for i in range(n):
        I[i,i] = -1
    for i in range(n):
        P[i,i] = 4
        if i > 0: P[i,i-1] = -1
    for i in range(0, n*n, n):
        L[i:i+n,i:i+n] = P
        if i > 0: L[i:i+n,i-n:i] = I
    return L

n = 200

t1 = time.clock()
L = poisson2d_sym_blk(n)
print 'Time for constructing the matrix: %8.2f sec' % (time.clock() - t1, )

print L.nnz

# ---------------------------------------------------------------------------------------
t1 = time.clock()
jdsym.jdsym(L.to_sss(), None, None, 5, 0.0, 1e-8, 100, itsolvers.qmrs, clvl=1)
print 'Time spend in jdsym: %8.2f sec' % (time.clock() - t1, )
