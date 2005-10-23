import Numeric
import math
import spmatrix
import itsolvers
import precon
import time

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

tol = 1e-8
n = 100

t1 = time.clock()
L = poisson2d_sym_blk(n)
print 'Time for constructing the matrix using poisson2d_sym_blk: %8.2f sec' % (time.clock() - t1, )

t1 = time.clock()
L = poisson2d_sym(n)
print 'Time for constructing the matrix using poisson2d_sym    : %8.2f sec' % (time.clock() - t1, )

t1 = time.clock()
L = poisson2d(n)
print 'Time for constructing the matrix using poisson2d        : %8.2f sec' % (time.clock() - t1, )


A = L.to_csr()
S = L.to_sss()
print L.nnz
print S.nnz
print A.nnz
b = Numeric.ones(n*n, 'd')

# ---------------------------------------------------------------------------------------

t1 = time.clock()

x = Numeric.zeros(n*n, 'd')
info, iter, relres = itsolvers.pcg(S, b, x, tol, 2000)
print 'info=%d, iter=%d, relres=%e' % (info, iter, relres)

print 'Time for solving the system using SSS matrix: %8.2f sec' % (time.clock() - t1, )

print 'norm(x) = %g' % math.sqrt(Numeric.dot(x, x))

r = Numeric.zeros(n*n, 'd')
S.matvec(x, r)
r = b - r
print 'norm(b - A*x) = %g' % math.sqrt(Numeric.dot(r, r))

print x[0:10]

# ---------------------------------------------------------------------------------------

t1 = time.clock()

x = Numeric.zeros(n*n, 'd')
info, iter, relres = itsolvers.pcg(A, b, x, tol, 2000)
print 'info=%d, iter=%d, relres=%e' % (info, iter, relres)

print 'Time for solving the system using CSR matrix: %8.2f sec' % (time.clock() - t1, )

print 'norm(x) = %g' % math.sqrt(Numeric.dot(x, x))

r = Numeric.zeros(n*n, 'd')
A.matvec(x, r)
r = b - r
print 'norm(b - A*x) = %g' % math.sqrt(Numeric.dot(r, r))

# ---------------------------------------------------------------------------------------

t1 = time.clock()

x = Numeric.zeros(n*n, 'd')
info, iter, relres = itsolvers.pcg(L, b, x, tol, 2000)
print 'info=%d, iter=%d, relres=%e' % (info, iter, relres)

print 'Time for solving the system using LL matrix: %8.2f sec' % (time.clock() - t1, )

print 'norm(x) = %g' % math.sqrt(Numeric.dot(x, x))

r = Numeric.zeros(n*n, 'd')
A.matvec(x, r)
r = b - r
print 'norm(b - A*x) = %g' % math.sqrt(Numeric.dot(r, r))

# ---------------------------------------------------------------------------------------

K_ssor = precon.ssor(S, 1.9)
t1 = time.clock()

x = Numeric.zeros(n*n, 'd')
info, iter, relres = itsolvers.pcg(S, b, x, tol, 2000, K_ssor)
print 'info=%d, iter=%d, relres=%e' % (info, iter, relres)

print 'Time for solving the system using SSS matrix and SSOR preconditioner: %8.2f sec' % (time.clock() - t1, )

print 'norm(x) = %g' % math.sqrt(Numeric.dot(x, x))

r = Numeric.zeros(n*n, 'd')
S.matvec(x, r)
r = b - r
print 'norm(b - A*x) = %g' % math.sqrt(Numeric.dot(r, r))

# ---------------------------------------------------------------------------------------

import jdsym
jdsym.jdsym(S, None, None, 5, 0.0, 1e-8, 100, itsolvers.qmrs, clvl=1)
