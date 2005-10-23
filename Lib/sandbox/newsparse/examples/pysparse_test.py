import math, os, sys, time
import Numeric
import spmatrix
import itsolvers
import precon

ll = spmatrix.ll_mat(5,5)
print ll
print ll[1,1]
print ll

ll[2,1] = 1.0
ll[1,3] = 2.0
print ll
print ll.to_csr()

print ll[1,3]
print ll[1,-1]
print ll.nnz

ll.export_mtx('test.mtx')

L = spmatrix.ll_mat(10, 10)
for i in range(0, 10):
    L[i,i] = float(i+1)

A = L.to_csr()
x = Numeric.ones([10], 'd')
y = Numeric.zeros([10], 'd')
print A, x, y
A.matvec(x, y)
print y

ll = spmatrix.ll_mat(100, 100)
for i in range(0, 100, 5):
    for j in range(0, 100, 4):
        ll[i,j] = 1.0/float(i+j+1)
A = ll.to_csr()

x = Numeric.arange(100).astype(Numeric.Float)
y = Numeric.zeros(100, 'd')
z = Numeric.zeros(100, 'd')

A.matvec(x, y)
print y
print 'norm(y) = ', math.sqrt(Numeric.add.reduce(y))

##A.matvec_transp(x, z)
##print z
##print 'norm(z) = ', math.sqrt(Numeric.add.reduce(z))

L = spmatrix.ll_mat(10,10)
for i in range(10):
    L[i,i] = float(i+1)
A = L.to_csr()
print A
x = Numeric.zeros(10, 'd')
b = Numeric.ones(10, 'd')
info, iter, relres = itsolvers.pcg(A, b, x, 1e-8, 100)
print info, iter, relres
print x
if (info != 0):
    print >> sys.stderr, 'cg not converged'

L2 = L.copy()
x = Numeric.zeros(10, 'd')
info, iter, relres = itsolvers.pcg(A, b, x, 1e-8, 100)
print info, iter, relres

# -----------------------------------------------------------
print 'remove test'
n = 100
L = spmatrix.ll_mat(n, n)

for run in range(5):

    print 'adding elements...'
    for i in range(0,n,2):
        for j in range (n):
            L[i,j] = i+j+1
            # print L

    print L.nnz

    print 'removing elements...'
    for j in range(0,n,2):
        for i in range (n):
            L[i,j] = 0.0
            # print L

    print L.nnz

# -----------------------------------------------------------
print 'submatrix test'
n = 100
L = spmatrix.ll_mat(n, n)

for i in range (0, n, 2):
    for j in range (1, n, 2):
        L[i,j] = float(n*i + j);
print L[10:18,75:80]
print L[10:15,35:10]
print L[19:15,35:10]

# -----------------------------------------------------------
print 'submatrix assign test'
n = 10
L = spmatrix.ll_mat(n, n);

for i in range (0, n, 1):
    for j in range (0, n, 1):
        L[i,j] = 1.0;

print L
Z = spmatrix.ll_mat(n-2, n-2)
L[1:n-1,1:n-1] = Z
print L
print L.nnz

#------------------------------------------------------------

if 0:
    f = open(os.environ['HOME']+'/matrices/poi2d_300.mtx')
    t1 = time.clock()
    L = ll_mat_from_mtx(f)
    t_read = time.clock() - t1
    f.close()
    print 'time for reading matrix data from file: %.2f sec' % t_read

if 1:
    t1 = time.clock()
    L = spmatrix.ll_mat_from_mtx(os.environ['HOME']+'/matrices/poi2d_300.mtx')
    t_read = time.clock() - t1
    print 'time for reading matrix data from file: %.2f sec' % t_read

#------------------------------------------------------------

L = spmatrix.ll_mat_from_mtx(os.environ['HOME']+'/matrices/node4x3x1_A.mtx')
print L.shape, L.nnz

A = L.to_sss()

class diag_prec:
    def __init__(self, A):
        self.shape = A.shape
        n = self.shape[0]
        self.dinv = Numeric.zeros(n, 'd')
        for i in xrange(n):
            self.dinv[i] = 1.0 / A[i,i]
    def precon(self, x, y):
        Numeric.multiply(x, self.dinv, y)

def resid(A, b, x):
    r = x.copy()
    A.matvec(x, r)
    r = b - r
    return math.sqrt(Numeric.dot(r, r))

K_diag = diag_prec(A)
K_jac = precon.jacobi(A, 1.0, 1)
K_ssor = precon.ssor(A, 1.0, 1)
# K_ilu = precon.ilutp(L)

n = L.shape[0];
b = Numeric.arange(n).astype(Numeric.Float)
x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.pcg(A, b, x, 1e-6, 1000)
print 'pcg, K_none: ', info, iter, relres, resid(A, b, x)
x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.pcg(A, b, x, 1e-6, 1000, K_diag)
print 'pcg, K_diag: ', info, iter, relres, resid(A, b, x)
x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.pcg(A, b, x, 1e-6, 1000, K_jac)
print 'pcg, K_jac: ', info, iter, relres, resid(A, b, x)
x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.pcg(A, b, x, 1e-6, 1000, K_ssor)
print 'pcg, K_ssor: ', info, iter, relres, resid(A, b, x)

x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.minres(A, b, x, 1e-6, 1000)
print 'minres, K_none: ', info, iter, relres, resid(A, b, x)
x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.minres(A, b, x, 1e-6, 1000, K_diag)
print 'minres, K_diag: ', info, iter, relres, resid(A, b, x)
x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.minres(A, b, x, 1e-6, 1000, K_jac)
print 'minres, K_jac: ', info, iter, relres, resid(A, b, x)
x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.minres(A, b, x, 1e-6, 1000, K_ssor)
print 'minres, K_ssor: ', info, iter, relres, resid(A, b, x)

x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.qmrs(A, b, x, 1e-6, 1000)
print 'qmrs, K_none: ', info, iter, relres, resid(A, b, x)
x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.qmrs(A, b, x, 1e-6, 1000, K_diag)
print 'qmrs, K_diag: ', info, iter, relres, resid(A, b, x)
x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.qmrs(A, b, x, 1e-6, 1000, K_jac)
print 'qmrs, K_jac: ', info, iter, relres, resid(A, b, x)
x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.qmrs(A, b, x, 1e-6, 1000, K_ssor)
print 'qmrs, K_ssor: ', info, iter, relres, resid(A, b, x)

x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.cgs(A, b, x, 1e-6, 1000)
print 'cgs, K_none: ', info, iter, relres, resid(A, b, x)
x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.cgs(A, b, x, 1e-6, 1000, K_diag)
print 'cgs, K_diag: ', info, iter, relres, resid(A, b, x)
x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.cgs(A, b, x, 1e-6, 1000, K_jac)
print 'cgs, K_jac: ', info, iter, relres, resid(A, b, x)
x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.cgs(A, b, x, 1e-6, 1000, K_ssor)
print 'cgs, K_ssor: ', info, iter, relres, resid(A, b, x)

x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.bicgstab(A, b, x, 1e-6, 1000)
print 'bicgstab, K_none: ', info, iter, relres, resid(A, b, x)
x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.bicgstab(A, b, x, 1e-6, 1000, K_diag)
print 'bicgstab, K_diag: ', info, iter, relres, resid(A, b, x)
x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.bicgstab(A, b, x, 1e-6, 1000, K_jac)
print 'bicgstab, K_jac: ', info, iter, relres, resid(A, b, x)
x = Numeric.zeros(n, 'd')
info, iter, relres = itsolvers.bicgstab(A, b, x, 1e-6, 1000, K_ssor)
print 'bicgstab, K_ssor: ', info, iter, relres, resid(A, b, x)

#------------------------------------------------------------

import superlu

L = spmatrix.ll_mat_from_mtx(os.environ['HOME']+'/matrices/cop18_el3_A.mtx')
##f = open('cop18_el5_A.mtx')
##L = ll_mat_from_mtx(f)
##f.close()
n11 = 4688
L = L[0:n11, 0:n11]                    # extract (1,1)-block

# make matrix regular
for i in xrange(n11):
    L[i,i] = 1

print L.shape, L.nnz
n = L.shape[0]

B = L.to_csr()
su = superlu.factorize(B, diag_pivot_thresh=0.0)
print su.nnz
b = Numeric.arange(n).astype(Numeric.Float) / n
x = Numeric.zeros(n, 'd')
su.solve(b, x)
print 'norm(b) = %g' % math.sqrt(Numeric.dot(b, b))
print 'norm(x) = %g' % math.sqrt(Numeric.dot(x, x))

r = Numeric.zeros(n, 'd')
B.matvec(x, r)
r = b - r
print 'norm(b - A*x) = %g' % math.sqrt(Numeric.dot(r, r))

if 1:
    for panel_size in [5, 10, 15]:
        for relax in [1, 3, 5]:
            for permc_spec in [0, 1, 2]:
                for diag_pivot_thresh in [0.0, 0.5, 1.0]:

                    t1 = time.clock()
                    su = superlu.factorize(B,
                                           panel_size=panel_size,
                                           relax=relax,
                                           permc_spec=permc_spec,
                                           diag_pivot_thresh=diag_pivot_thresh)
                    t_fact = time.clock() - t1

                    t1 = time.clock()
                    su.solve(b, x)
                    t_solve = time.clock() - t1

                    print 'panel_size=%2d, relax=%d, permc_spec=%d, diag_pivot_thresh=%.1f   nnz=%d, t_fact=%.2f, t_solve=%.2f' % \
                          (panel_size, relax, permc_spec, diag_pivot_thresh, su.nnz, t_fact, t_solve)
