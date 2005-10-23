import math, unittest
import Numeric
import spmatrix, superlu
import poisson

def residual(A, x, b):
    n = A.shape[0]
    r = Numeric.zeros(n, 'd')
    A.matvec(x, r)
    r -= b
    return math.sqrt(Numeric.dot(r, r))

def error(x, y):
    n = len(x)
    t = x.copy(); t -= y
    return math.sqrt(Numeric.dot(t, t))

def speye(n):
    A = spmatrix.ll_mat_sym(n, n)
    for i in xrange(n):
        A[i,i] = 1.0
    return A

def macheps():
    "compute machine epsilon"
    eps = 1.0
    while (1.0 + eps > 1.0):
        eps /= 2.0
    return 2.0 * eps
    
class EyeTestCase(unittest.TestCase):
    def setUp(self):
        self.n = 10000
        self.A = speye(self.n).to_csr()
        self.b = Numeric.ones(self.n, 'd')
        self.x = Numeric.zeros(self.n, 'd')

    def testTrivial(self):
        luA = superlu.factorize(self.A)
        luA.solve(self.b, self.x)
        self.failUnless(residual(self.A, self.x, self.b) == 0.0)

        luA = superlu.factorize(self.A, diag_pivot_thresh=0.0)
        luA.solve(self.b, self.x)
        self.failUnless(residual(self.A, self.x, self.b) == 0.0)

        luA = superlu.factorize(self.A, relax=20)
        luA.solve(self.b, self.x)
        self.failUnless(residual(self.A, self.x, self.b) == 0.0)

        luA = superlu.factorize(self.A, panel_size=1)
        luA.solve(self.b, self.x)
        self.failUnless(residual(self.A, self.x, self.b) == 0.0)

        for permc_spec in [0, 1, 2, 3]:
            luA = superlu.factorize(self.A, permc_spec=permc_spec)
            luA.solve(self.b, self.x)
            self.failUnless(residual(self.A, self.x, self.b) == 0.0)
            
class Poisson1dTestCase(unittest.TestCase):
    def setUp(self):
        self.n = 50000
        self.B = poisson.poisson1d(self.n).to_csr()
        
        self.b = Numeric.zeros(self.n, 'd')
        self.x = Numeric.zeros(self.n, 'd')
        self.x_exact = Numeric.ones(self.n, 'd')
        self.x_exact /= math.sqrt(self.n)
        self.B.matvec(self.x_exact, self.b)
        
        lmbd_min = 4.0 * math.sin(math.pi/2.0/self.n) ** 2
        lmbd_max = 4.0 * math.sin((self.n - 1)*math.pi/2.0/self.n) ** 2
        cond = lmbd_max/lmbd_min
        self.tol = cond * macheps()
        
    def testPoisson1dDefault(self):
        luA = superlu.factorize(self.B)
        luA.solve(self.b, self.x)
        print error(self.x, self.x_exact), self.tol, luA.nnz
        self.failUnless(error(self.x, self.x_exact) < self.tol)

    def testPoisson1dNoPivot(self):
        luA = superlu.factorize(self.B, diag_pivot_thresh=0.0)
        luA.solve(self.b, self.x)
        print error(self.x, self.x_exact), self.tol, luA.nnz
        self.failUnless(error(self.x, self.x_exact) < self.tol)

    def testPoisson1dRelax(self):
        luA = superlu.factorize(self.B, relax=20)
        luA.solve(self.b, self.x)
        print error(self.x, self.x_exact), self.tol, luA.nnz
        self.failUnless(error(self.x, self.x_exact) < self.tol)

    def testPoisson1dSmallPanel(self):
        luA = superlu.factorize(self.B, panel_size=1)
        luA.solve(self.b, self.x)
        print error(self.x, self.x_exact), self.tol, luA.nnz
        self.failUnless(error(self.x, self.x_exact) < self.tol)
        
    def testPoisson1dOrigOrdering(self):
        luA = superlu.factorize(self.B, permc_spec=0)
        luA.solve(self.b, self.x)
        print error(self.x, self.x_exact), self.tol, luA.nnz
        self.failUnless(error(self.x, self.x_exact) < self.tol)

    def testPoisson1dMMD_AtimesA(self):
        luA = superlu.factorize(self.B, permc_spec=1)
        luA.solve(self.b, self.x)
        print error(self.x, self.x_exact), self.tol, luA.nnz
        self.failUnless(error(self.x, self.x_exact) < self.tol)

    def testPoisson1dMMD_AplusA(self):
        luA = superlu.factorize(self.B, permc_spec=2)
        luA.solve(self.b, self.x)
        print error(self.x, self.x_exact), self.tol, luA.nnz
        self.failUnless(error(self.x, self.x_exact) < self.tol)

    def testPoisson1dCOLAMD(self):
        luA = superlu.factorize(self.B, permc_spec=3)
        luA.solve(self.b, self.x)
        print error(self.x, self.x_exact), self.tol, luA.nnz
        self.failUnless(error(self.x, self.x_exact) < self.tol)
        
class Poisson2dTestCase(unittest.TestCase):
    def setUp(self):
        self.n = 200
        self.B = poisson.poisson2d(self.n).to_csr()
        
        self.b = Numeric.zeros(self.n*self.n, 'd')
        self.x = Numeric.zeros(self.n*self.n, 'd')
        self.x_exact = Numeric.ones(self.n*self.n, 'd')
        self.x_exact /= math.sqrt(self.n*self.n)
        self.B.matvec(self.x_exact, self.b)

        h = 1.0 / self.n
        lmbd_min = 4.0/h/h * (math.sin(math.pi*h/2.0) ** 2 +
                              math.sin(math.pi*h/2.0) ** 2)
        lmbd_max = 4.0/h/h * (math.sin((self.n - 1)*math.pi*h/2.0) ** 2 +
                              math.sin((self.n - 1)*math.pi*h/2.0) ** 2)
        cond = lmbd_max/lmbd_min
        self.tol = cond * macheps()
        
    def testPoisson2dDefault(self):
        luA = superlu.factorize(self.B)
        luA.solve(self.b, self.x)
        print error(self.x, self.x_exact), self.tol, luA.nnz
        self.failUnless(error(self.x, self.x_exact) < self.tol)

    def testPoisson2dNoPivot(self):
        luA = superlu.factorize(self.B, diag_pivot_thresh=0.0)
        luA.solve(self.b, self.x)
        print error(self.x, self.x_exact), self.tol, luA.nnz
        self.failUnless(error(self.x, self.x_exact) < self.tol)

    def testPoisson2dRelax(self):
        luA = superlu.factorize(self.B, relax=20)
        luA.solve(self.b, self.x)
        print error(self.x, self.x_exact), self.tol, luA.nnz
        self.failUnless(error(self.x, self.x_exact) < self.tol)

    def testPoisson2dSmallPanel(self):
        luA = superlu.factorize(self.B, panel_size=1)
        luA.solve(self.b, self.x)
        print error(self.x, self.x_exact), self.tol, luA.nnz
        self.failUnless(error(self.x, self.x_exact) < self.tol)
        
    def testPoisson2dOrigOrdering(self):
        luA = superlu.factorize(self.B, permc_spec=0)
        luA.solve(self.b, self.x)
        print error(self.x, self.x_exact), self.tol, luA.nnz
        self.failUnless(error(self.x, self.x_exact) < self.tol)

    def testPoisson2dMMD_AtimesA(self):
        luA = superlu.factorize(self.B, permc_spec=1)
        luA.solve(self.b, self.x)
        print error(self.x, self.x_exact), self.tol, luA.nnz
        self.failUnless(error(self.x, self.x_exact) < self.tol)

    def testPoisson2dMMD_AplusA(self):
        luA = superlu.factorize(self.B, permc_spec=2)
        luA.solve(self.b, self.x)
        print error(self.x, self.x_exact), self.tol, luA.nnz
        self.failUnless(error(self.x, self.x_exact) < self.tol)

    def testPoisson2dCOLAMD(self):
        luA = superlu.factorize(self.B, permc_spec=3)
        luA.solve(self.b, self.x)
        print error(self.x, self.x_exact), self.tol, luA.nnz
        self.failUnless(error(self.x, self.x_exact) < self.tol)

if __name__ == '__main__':
    unittest.main()
