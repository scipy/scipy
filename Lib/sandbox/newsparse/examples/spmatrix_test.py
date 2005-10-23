import spmatrix
import Numeric
import traceback
import spmatrix_util

def printMatrix(M):
    n, m = M.shape
    Z = Numeric.zeros((n,m), 'd')
    for i in range(n):
        for j in range(m):
            Z[i,j] = M[i,j]
    print str(Z) + '\n'
    
n = 10
A = spmatrix.ll_mat(n,n)
As = spmatrix.ll_mat_sym(n)
Is = spmatrix.ll_mat_sym(n)
I = spmatrix.ll_mat(n,n)
Os = spmatrix.ll_mat_sym(n)
O = spmatrix.ll_mat(n,n)

for i in range(n):
    for j in range(n):
        if i >= j:
            A[i,j] = 10*i + j
        else:
            A[i,j] = 10*j + i
        O[i,j] = 1
            
for i in range(n):
    for j in range(n):
        if i >= j:
            As[i,j] = 10*i + j
            Os[i,j] = 1

for i in range(n):
    I[i,i] = 1
    Is[i,i] = 1

print 'Setting matrix elements'
printMatrix(A)
printMatrix(As)

print 'Extracting submatrices'
printMatrix(A[4:8,1:3])
printMatrix(As[4:8,1:3])
printMatrix(A[1:3,4:8])
printMatrix(As[1:3,4:8])

printMatrix(A[6:9,6:9])
printMatrix(As[6:9,6:9])
print As[5:9,5:9]
print

print 'this should raise an execption...\n'
try:
    As[5:9, 4:10]
except:
    traceback.print_exc()
   
print 'Setting submatrices'
T = spmatrix.ll_mat_sym(n)
T[6:9,6:9] = As[6:9,6:9]
T[4:8,1:3] = As[4:8,1:3]
printMatrix(T)

print 'this should raise execptions...\n'
try:
    T[6:9,6:9] = A[6:9,6:9]
except:
    traceback.print_exc()
try:
    T[5:9, 4:10] = A[5:9, 4:10]
except:
    traceback.print_exc()

print 'Matrix multiplications'

printMatrix(spmatrix.matrixmultiply(I, A))
printMatrix(spmatrix.matrixmultiply(Is, A))

printMatrix(spmatrix.matrixmultiply(O, O))
printMatrix(spmatrix.matrixmultiply(Os, O))

print 'Dot product'
printMatrix(spmatrix.dot(I, A))

print 'Matrix export'
A[:4,:4].export_mtx('A.mtx', 3)
As[:4,:4].export_mtx('As.mtx', 3)

print open('A.mtx').read()
print open('As.mtx').read()

print 'Matrix import'
printMatrix(spmatrix.ll_mat_from_mtx('A.mtx'))
printMatrix(spmatrix.ll_mat_from_mtx('As.mtx'))

print 'Conversion to CSR'
print A[:4,:4]
print A[:4,:4].to_csr()
print As[:4,:4].to_csr()

print 'update_add_mask operations'
ind = Numeric.array([3, 4, 5, 6], 'i')
mask = Numeric.array([1, 1, 1, 1], 'i')
B = Numeric.ones((4,4), 'd')
Ac = A.copy()
Ac.update_add_mask(B, ind, ind, mask, mask)
A.update_add_mask_sym(B, ind, mask)
As.update_add_mask_sym(B, ind, mask)
printMatrix(Ac[2:8,2:8])
printMatrix(A[2:8,2:8])
printMatrix(As[2:8,2:8])

print 'deleting rows'

Atemp = A.copy()
print 'original matrix:'
printMatrix(Atemp)

print 'Matrix with rows 7 and 8 and deleted:'
mask = Numeric.ones(n, 'l')
mask[7:9] = 0
Atemp.delete_rows(mask)
printMatrix(Atemp)
print Atemp.delete_rows.__doc__

print 'deleting rows and column'

Atemp = As.copy()
print 'original matrix:'
printMatrix(Atemp)

print 'Matrix with rows/cols 7 and 8 and deleted:'
mask = Numeric.ones(n, 'l')
mask[7:9] = 0
Atemp.delete_rowcols(mask)
printMatrix(Atemp)

nn = 100
R = spmatrix_util.ll_mat_rand(nn, nn, 0.3)
##print R.nnz
for i in range(nn-5):
    mask = Numeric.ones(nn, 'l')
    mask[0] = 0
    R.delete_rowcols(mask)
    nn -= 1
##print R
