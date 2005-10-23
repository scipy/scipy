from spmatrix import *
import RandomArray
import time

n = 1000
nnz = 50000
A = ll_mat(n, n, nnz)

R = RandomArray.randint(0, n, (nnz,2))

t1 = time.clock()

for k in xrange(nnz):
    A[R[k,0],R[k,1]] = k
    
print 'Time for populating matrix: %8.2f sec' % (time.clock() - t1, )

print A.nnz

B = A[:,:]
A.shift(-1.0, B)
print A

