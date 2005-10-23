import unittest
import math, random
import spmatrix
import spmatrix_util
import poisson
import Numeric, RandomArray

def llmat_isEqual(aMat, bMat):
    if aMat.issym and not bMat.issym:
        temp = aMat; aMat = bMat; bMat = temp
    zMat = aMat.copy()
    zMat.shift(-1.0, bMat)
    return zMat.nnz == 0

class LLMatSimpleTestCase(unittest.TestCase):
    def setUp(self):
        self.n = 10
        self.A = spmatrix.ll_mat(self.n, self.n)
        self.S = spmatrix.ll_mat_sym(self.n)

    def testCreate(self):
        self.failUnless(self.A.shape == (self.n, self.n))
        self.failUnless(self.A.nnz == 0)
        self.failUnless(not self.A.issym)
        self.failUnless(self.S.shape == (self.n, self.n))
        self.failUnless(self.S.nnz == 0)
        self.failUnless(self.S.issym)
        
    def testEntry(self):
        def assignUP(): self.S[0,1] = 1.0
        def assignLeft(): self.S[-11,0] = 1.0
        def assignRight(): self.S[10,0] = 1.0
        def assignTop(): self.S[0,-11] = 1.0
        def assignBottom(): self.S[0,10] = 1.0
        
        self.A[0,0] = 1.0
        self.S[0,0] = 1.0
        self.failUnless(self.A[0,0] == 1.0)
        self.failUnless(self.A.nnz == 1)
        self.failUnless(self.S[0,0] == 1.0)
        self.failUnless(self.S.nnz == 1)
        self.failUnlessRaises(spmatrix.error, assignUP)
        self.A[0,0] += 1.0
        self.failUnless(self.A[0,0] == 2.0)
        self.failUnless(self.A.nnz == 1)
        self.A[0,0] -= 2.0
        self.failUnless(self.A[0,0] == 0.0)
        self.failUnless(self.A.nnz == 0)
        # indices out of bounds
        for f in [assignLeft, assignRight, assignTop, assignBottom]:
            self.failUnlessRaises(IndexError, f)
        # negative indices
        I = spmatrix.ll_mat(10, 10, 100)
        for i in range(10):
            for j in range(10):
                I[i,j] = 10*i + j
        for i in range(-10, 0):
            for j in range(-10, 0):
                self.failUnless(I[i,j] == I[10+i,10+j])
        

class LLMatPoissonTestCase(unittest.TestCase):
    
    def setUp(self):
        self.n = 20
        self.A = poisson.poisson2d(self.n)
        self.S = poisson.poisson2d_sym(self.n)
        self.B = poisson.poisson2d_sym_blk(self.n)

    def testBasic(self):
        self.failUnless(self.S.nnz == self.n*(3*self.n - 2))
        self.failUnless(self.A.nnz == self.n*(5*self.n - 4))
        self.failUnless(llmat_isEqual(self.A, self.A))
        self.failUnless(llmat_isEqual(self.S, self.S))
        self.failUnless(llmat_isEqual(self.A, self.S))
        self.failUnless(llmat_isEqual(self.A, self.B))

    def testSubmatrix(self):
        n = self.n
        Psym = poisson.poisson1d_sym(n)
        P = poisson.poisson1d(n)
        for i in range(n):
            P[i,i] = 4.0
            Psym[i,i] = 4.0
        # read and test diagonal blocks
        for i in range(n):
            self.failUnless(llmat_isEqual(self.A[n*i:n*(i+1),n*i:n*(i+1)], P))
            self.failUnless(llmat_isEqual(self.S[n*i:n*(i+1),n*i:n*(i+1)], P))
            self.failUnless(llmat_isEqual(self.A[n*i:n*(i+1),n*i:n*(i+1)], Psym))
            self.failUnless(llmat_isEqual(self.S[n*i:n*(i+1),n*i:n*(i+1)], Psym))
        # store and get diagonal blocks
        R = spmatrix_util.ll_mat_rand(n*n, n*n, 0.01) # random matrix
        for i in range(n):
            R[n*i:n*(i+1),n*i:n*(i+1)] = P
            self.failUnless(llmat_isEqual(R[n*i:n*(i+1),n*i:n*(i+1)], P))
            R[n*i:n*(i+1),n*i:n*(i+1)] = Psym
            self.failUnless(llmat_isEqual(R[n*i:n*(i+1),n*i:n*(i+1)], Psym))
        # store and get off-diagonal blocks
        for i in range(n-1):
            R[n*i:n*(i+1),n*(i+1):n*(i+2)] = P
            self.failUnless(llmat_isEqual(R[n*i:n*(i+1),n*(i+1):n*(i+2)], P))
            R[n*i:n*(i+1),n*(i+1):n*(i+2)] = Psym
            self.failUnless(llmat_isEqual(R[n*i:n*(i+1),n*(i+1):n*(i+2)], Psym))
        # store and get diagonal blocks in symmetric matrix
        R = spmatrix.ll_mat_sym(n*n)
        for i in range(n):
            R[n*i:n*(i+1),n*i:n*(i+1)] = Psym
            self.failUnless(llmat_isEqual(R[n*i:n*(i+1),n*i:n*(i+1)], Psym))
        # store and get off-diagonal blocks in symmetric matrix
        for i in range(n-1):
            R[n*(i+1):n*(i+2),n*i:n*(i+1)] = P
            self.failUnless(llmat_isEqual(R[n*(i+1):n*(i+2),n*i:n*(i+1)], P))
            R[n*(i+1):n*(i+2),n*i:n*(i+1)] = Psym
            self.failUnless(llmat_isEqual(R[n*(i+1):n*(i+2),n*i:n*(i+1)], Psym))

class LLMatDeleteRowColsTestCase(unittest.TestCase):
    
    def setUp(self):
        import Numeric
        
        self.n = 30
        self.P = poisson.poisson1d(self.n)
        for i in range(self.n):
            self.P[i,i] = 4.0
        self.A = poisson.poisson2d(self.n)
        self.S = poisson.poisson2d_sym(self.n)
        self.I = spmatrix.ll_mat_sym(self.n)
        for i in range(self.n):
            self.I[i,i] = -1.0
        self.mask = Numeric.zeros(self.n**2, 'l')
        self.mask[self.n/2*self.n:(self.n/2 + 1)*self.n] = 1
        self.mask1 = Numeric.zeros(self.n**2, 'l')
        self.mask1[(self.n/2 + 1)*self.n:(self.n/2 + 2)*self.n] = 1

    def testDeleteRowColsSym(self):
        self.S.delete_rowcols(self.mask)
        self.failUnless(llmat_isEqual(self.S, self.P))
        
    def testDeleteRowColsGen(self):
        self.A.delete_rowcols(self.mask)
        self.failUnless(llmat_isEqual(self.A, self.P))
        
    def testDeleteRowColsGen2Step(self):
        self.A.delete_rows(self.mask)
        self.A.delete_cols(self.mask)
        self.failUnless(llmat_isEqual(self.A, self.P))
        
    def testDeleteRowColsGen2StepOff(self):
        self.A.delete_rows(self.mask)
        self.A.delete_cols(self.mask1)
        self.failUnless(llmat_isEqual(self.A, self.I))

    def testCompress(self):
        self.A.delete_rows(self.mask)
        self.A.delete_cols(self.mask1)
        norm1 = self.A.norm('fro')
        self.A.compress()
        norm2 = self.A.norm('fro')
        self.failUnless(norm1 == norm2)

    def testCompressStress(self):
        n = 20
        A = spmatrix.ll_mat(n, n)
        for k in range(20):
            for i in range(n*n/2):
                i = random.randrange(n)
                j = random.randrange(n)
                A[i, j] = 1.0
            for i in range(n*n/2):
                i = random.randrange(n)
                j = random.randrange(n)
                A[i, j] = 0.0
        
class LLMatNorm(unittest.TestCase):
    def setUp(self):
        self.n = 30

    def testNormGeneral(self):
        A = poisson.poisson2d(self.n)
        self.failUnless(A.norm('1') == 8)
        self.failUnless(A.norm('inf') == 8)
        self.failUnless(poisson.poisson1d(3).norm('fro') == 4)
            
    def testNormSymmetric(self):
        A = spmatrix.ll_mat_sym(4)
        A[0,0] = 1; A[1,1] = 2; A[2,2] = 3; A[3,3] = 4;
        A[1,0] = 3; A[2,0] = 2; A[3,0] = 2; 
        self.failUnless(A.norm('fro') == 8)
        
    def testNormSymmetricNotImplemented(self):
        def f(): return A.norm('1')
        def g(): return A.norm('inf')
        
        A = poisson.poisson2d_sym(self.n)
        self.failUnlessRaises(NotImplementedError, f)
        self.failUnlessRaises(NotImplementedError, g)

class LLMatMatMul(unittest.TestCase):
      def testRandomMat(self):
          eps = 2.2204460492503131E-16
          n = 30; m = 60; k = 30

          for i in range(100):
              A = spmatrix_util.ll_mat_rand(n, k, 0.9)
              B = spmatrix_util.ll_mat_rand(k, m, 0.4)
              C = spmatrix.matrixmultiply(A, B)
              t = Numeric.zeros(k, 'd')
              y1 = Numeric.zeros(n, 'd')
              y2 = Numeric.zeros(n, 'd')
              for s in range(10):
                  x = RandomArray.random((m, ))
                  C.matvec(x, y1)
                  B.matvec(x, t)
                  A.matvec(t, y2)
                  self.failUnless(math.sqrt(Numeric.dot(y1 - y2, y1 - y2)) < eps * n*m*k)
              
if __name__ == '__main__':
    unittest.main()
