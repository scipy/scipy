from scipy.testing import *

import numpy
import scipy
from scipy.sparse import spdiags, csr_matrix
from scipy import arange, ones, zeros, array, allclose, zeros_like, \
        tril, diag, triu, rand, asmatrix
from scipy.linalg import solve



import scipy.sandbox.multigrid
from scipy.sandbox.multigrid.gallery    import poisson
from scipy.sandbox.multigrid.relaxation import polynomial_smoother,gauss_seidel,jacobi




class TestRelaxation(TestCase):
    def test_polynomial(self):
        N  = 3
        A  = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N,format='csr')
        x0 = arange(N).astype(numpy.float64)
        x  = x0.copy()
        b  = zeros(N)

        r = (b - A*x0)
        polynomial_smoother(A,x,b,[-1.0/3.0])

        assert_almost_equal(x,x0-1.0/3.0*r)

        x  = x0.copy()
        polynomial_smoother(A,x,b,[0.2,-1])
        assert_almost_equal(x,x0 + 0.2*A*r - r)

        x  = x0.copy()
        polynomial_smoother(A,x,b,[0.2,-1])
        assert_almost_equal(x,x0 + 0.2*A*r - r)

        x  = x0.copy()
        polynomial_smoother(A,x,b,[-0.14285714,  1., -2.])
        assert_almost_equal(x,x0 - 0.14285714*A*A*r + A*r - 2*r)

    def test_jacobi(self):
        N = 1
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N,format='csr')
        x = arange(N).astype(numpy.float64)
        b = zeros(N)
        jacobi(A,x,b)
        assert_almost_equal(x,array([0]))

        N = 3
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N,format='csr')
        x = zeros(N)
        b = arange(N).astype(numpy.float64)
        jacobi(A,x,b)
        assert_almost_equal(x,array([0.0,0.5,1.0]))

        N = 3
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N,format='csr')
        x = arange(N).astype(numpy.float64)
        b = zeros(N)
        jacobi(A,x,b)
        assert_almost_equal(x,array([0.5,1.0,0.5]))

        N = 1
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N,format='csr')
        x = arange(N).astype(numpy.float64)
        b = array([10])
        jacobi(A,x,b)
        assert_almost_equal(x,array([5]))

        N = 3
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N,format='csr')
        x = arange(N).astype(numpy.float64)
        b = array([10,20,30])
        jacobi(A,x,b)
        assert_almost_equal(x,array([5.5,11.0,15.5]))

        N = 3
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N,format='csr')
        x = arange(N).astype(numpy.float64)
        x_copy = x.copy()
        b = array([10,20,30])
        jacobi(A,x,b,omega=1.0/3.0)
        assert_almost_equal(x,2.0/3.0*x_copy + 1.0/3.0*array([5.5,11.0,15.5]))

    def test_gauss_seidel_bsr(self):
        cases = []

        for N in [1,2,3,4,5,6,10]:
            A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).tocsr()
            
            divisors = [ n for n in range(1,N+1) if N % n == 0 ]

            x_csr = arange(N).astype(numpy.float64)
            b = x_csr**2
            gauss_seidel(A,x_csr,b)

            for D in divisors:
                B = A.tobsr(blocksize=(D,D))
                x_bsr = arange(N).astype(numpy.float64)
                gauss_seidel(B,x_bsr,b)
                assert_almost_equal(x_bsr,x_csr)
               
    def test_gauss_seidel_new(self):
        scipy.random.seed(0)

        cases = []
        cases.append( poisson( (4,), format='csr' ) )
        cases.append( poisson( (4,4), format='csr' ) )

        temp = asmatrix( rand(4,4) )
        cases.append( csr_matrix( temp.T * temp) )

        # reference implementation
        def gold(A,x,b,iterations,sweep):
            A = A.todense()

            L = tril(A,k=-1)
            D = diag(diag(A))
            U = triu(A,k=1)

            for i in range(iterations):
                if sweep == 'forward':
                    x = solve(L + D, (b - U*x) )
                elif sweep == 'backward':
                    x = solve(U + D, (b - L*x) )
                else:
                    x = solve(L + D, (b - U*x) )
                    x = solve(U + D, (b - L*x) )
            return x            


        for A in cases:

            b = asmatrix(rand(A.shape[0],1))
            x = asmatrix(rand(A.shape[0],1))

            x_copy = x.copy()
            gauss_seidel(A, x, b, iterations=1, sweep='forward')
            assert_almost_equal( x, gold(A,x_copy,b,iterations=1,sweep='forward') )
            
            x_copy = x.copy()
            gauss_seidel(A, x, b, iterations=1, sweep='backward')
            assert_almost_equal( x, gold(A,x_copy,b,iterations=1,sweep='backward') )
            
            x_copy = x.copy()
            gauss_seidel(A, x, b, iterations=1, sweep='symmetric')
            assert_almost_equal( x, gold(A,x_copy,b,iterations=1,sweep='symmetric') )



    def test_gauss_seidel_csr(self):
        N = 1
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N,format='csr')
        x = arange(N).astype(numpy.float64)
        b = zeros(N)
        gauss_seidel(A,x,b)
        assert_almost_equal(x,array([0]))

        N = 3
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N,format='csr')
        x = arange(N).astype(numpy.float64)
        b = zeros(N)
        gauss_seidel(A,x,b)
        assert_almost_equal(x,array([1.0/2.0,5.0/4.0,5.0/8.0]))

        N = 1
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N,format='csr')
        x = arange(N).astype(numpy.float64)
        b = zeros(N)
        gauss_seidel(A,x,b,sweep='backward')
        assert_almost_equal(x,array([0]))

        N = 3
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N,format='csr')
        x = arange(N).astype(numpy.float64)
        b = zeros(N)
        gauss_seidel(A,x,b,sweep='backward')
        assert_almost_equal(x,array([1.0/8.0,1.0/4.0,1.0/2.0]))

        N = 1
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N,format='csr')
        x = arange(N).astype(numpy.float64)
        b = array([10])
        gauss_seidel(A,x,b)
        assert_almost_equal(x,array([5]))

        N = 3
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N,format='csr')
        x = arange(N).astype(numpy.float64)
        b = array([10,20,30])
        gauss_seidel(A,x,b)
        assert_almost_equal(x,array([11.0/2.0,55.0/4,175.0/8.0]))


        #forward and backward passes should give same result with x=ones(N),b=zeros(N)
        N = 100
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N,format='csr')
        x = ones(N)
        b = zeros(N)
        gauss_seidel(A,x,b,iterations=200,sweep='forward')
        resid1 = numpy.linalg.norm(A*x,2)
        x = ones(N)
        gauss_seidel(A,x,b,iterations=200,sweep='backward')
        resid2 = numpy.linalg.norm(A*x,2)
        self.assert_(resid1 < 0.01 and resid2 < 0.01)
        self.assert_(allclose(resid1,resid2))

if __name__ == '__main__':
    nose.run(argv=['', __file__])
