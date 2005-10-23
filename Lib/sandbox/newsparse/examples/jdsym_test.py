import spmatrix, jdsym, itsolvers
from Numeric import zeros, dot, allclose, multiply
from math import sqrt
import RandomArray

class diagPrecShifted:
    def __init__(self, A, M, sigma):
        self.shape = A.shape
        n = self.shape[0]
        self.dinv = zeros(n, 'd')
        for i in xrange(n):
            self.dinv[i] = 1.0 / (A[i,i] - sigma*M[i,i])
    def precon(self, x, y):
        multiply(x, self.dinv, y)
    
def computeResiduals(A, M, lmbd, Q):
    kconv = lmbd.shape[0]
    residuals = zeros((kconv, ), 'd')
    r = zeros((n, ), 'd')
    u = zeros((n, ), 'd')
    t = zeros((n, ), 'd')
    for k in xrange(kconv):
        u = Q[:,k].copy()
        A.matvec(u, r)
        if M <> None:
            M.matvec(u, t)
        else:
            t = u
        r = r - lmbd[k]*t
        residuals[k] = sqrt(dot(r,r))
    return residuals
    
n = 1000; ncv = 5; tol = 1e-6

A = spmatrix.ll_mat_sym(n)
for i in xrange(n):
    A[i,i] = i+1.0
As = A.to_sss()

M = spmatrix.ll_mat_sym(n)
for i in xrange(n):
    M[i,i] = float(n/2) + i
Ms = M.to_sss()
normM = M[n-1,n-1]

K = diagPrecShifted(A, M, 0.006)

#-------------------------------------------------------------------------------
# Test 1: M = K = None

print 'Test 1'

lmbd_exact = zeros(ncv, 'd')
for k in xrange(ncv):
    lmbd_exact[k] =  A[k,k]

kconv, lmbd, Q, it, it_inner = jdsym.jdsym(As, None, None, ncv, 0.0, tol, 150, itsolvers.qmrs,
                                           jmin=5, jmax=10, eps_tr=1e-4, clvl=1)
    
assert ncv == kconv
assert allclose(computeResiduals(As, None, lmbd, Q), zeros(kconv), 0.0, tol)
assert allclose(lmbd, lmbd_exact, tol*tol, 0.0)

print 'OK'

#-------------------------------------------------------------------------------
# Test 2: K = None

print 'Test 2',

lmbd_exact = zeros(ncv, 'd')
for k in xrange(ncv):
    lmbd_exact[k] =  A[k,k]/M[k,k]


X0 = RandomArray.random((n,ncv))

kconv, lmbd, Q, it, it_inner = jdsym.jdsym(As, Ms, None, ncv, 0.0, tol, 150, itsolvers.qmrs,
                                           jmin=5, jmax=10, eps_tr=1e-4, clvl=1)
    
assert ncv == kconv
assert allclose(computeResiduals(As, Ms, lmbd, Q), zeros(kconv), 0.0, normM*tol)
assert allclose(lmbd, lmbd_exact, normM*tol*tol, 0.0)

print 'OK'

#-------------------------------------------------------------------------------
# Test 3: general case

print 'Test 3',

lmbd_exact = zeros(ncv, 'd')
for k in xrange(ncv):
    lmbd_exact[k] =  A[k,k]/M[k,k]

kconv, lmbd, Q, it, it_inner = jdsym.jdsym(As, Ms, K, ncv, 0.0, tol, 150, itsolvers.qmrs,
                                           jmin=5, jmax=10, eps_tr=1e-4, clvl=1)
assert ncv == kconv
assert allclose(computeResiduals(As, Ms, lmbd, Q), zeros(kconv), 0.0, normM*tol)
assert allclose(lmbd, lmbd_exact, normM*tol*tol, 0.0)

print 'OK'

#-------------------------------------------------------------------------------
# Test 4: K = None, with X0

print 'Test 4',

lmbd_exact = zeros(ncv, 'd')
for k in xrange(ncv):
    lmbd_exact[k] =  A[k,k]/M[k,k]

# Fixme: RandomArray.random is broken AMD64
# X0 = RandomArray.random((n,ncv))
X0 = zeros((n,ncv), 'd')
for k in xrange(ncv):
    X0[k,k] = 10000

kconv, lmbd, Q, it, it_inner = jdsym.jdsym(As, Ms, None, ncv, 0.0, tol, 150, itsolvers.qmrs,
                                           jmin=5, jmax=10, eps_tr=1e-4, clvl=1, V0=X0)

assert ncv == kconv
assert allclose(computeResiduals(As, Ms, lmbd, Q), zeros(kconv), 0.0, normM*tol)
assert allclose(lmbd, lmbd_exact, normM*tol*tol, 0.0)

print 'OK'
