import time
import spmatrix


n = 1000000

# create 2 nxn tridiag matrices
A = spmatrix.ll_mat(n, n)
B = spmatrix.ll_mat(n, n)

for i in xrange(n):
    A[i,i] = i
    B[i,i] = i
    if i > 0:
        A[i,i-1] = 1
        B[i,i-1] = 1
    if i < n-1:
        A[i,i+1] = 1
        B[i,i-1] = 1


t1 = time.clock()
C = spmatrix.matrixmultiply(A, B)
t_mult = time.clock() - t1
print 'time for multiplying %dx%d matrices: %.2f sec' % (n, n, t_mult)
print C[:10,:10]
