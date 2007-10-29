from numpy.testing import *

import numpy
import scipy
from scipy import arange,ones,zeros,array,allclose
from scipy.sparse import spdiags


set_package_path()
import scipy.sandbox.multigrid
from scipy.sandbox.multigrid.relaxation import polynomial_smoother,gauss_seidel,jacobi
restore_path()


class TestRelaxation(NumpyTestCase):
    def check_polynomial(self):
        N  = 3
        A  = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).T
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

    def check_jacobi(self):
        N = 1
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).T
        x = arange(N).astype(numpy.float64)
        b = zeros(N)
        jacobi(A,x,b)
        assert_almost_equal(x,array([0]))

        N = 3
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).T
        x = zeros(N)
        b = arange(N).astype(numpy.float64)
        jacobi(A,x,b)
        assert_almost_equal(x,array([0.0,0.5,1.0]))

        N = 3
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).T
        x = arange(N).astype(numpy.float64)
        b = zeros(N)
        jacobi(A,x,b)
        assert_almost_equal(x,array([0.5,1.0,0.5]))

        N = 1
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).T
        x = arange(N).astype(numpy.float64)
        b = array([10])
        jacobi(A,x,b)
        assert_almost_equal(x,array([5]))

        N = 3
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).T
        x = arange(N).astype(numpy.float64)
        b = array([10,20,30])
        jacobi(A,x,b)
        assert_almost_equal(x,array([5.5,11.0,15.5]))

        N = 3
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).T
        x = arange(N).astype(numpy.float64)
        x_copy = x.copy()
        b = array([10,20,30])
        jacobi(A,x,b,omega=1.0/3.0)
        assert_almost_equal(x,2.0/3.0*x_copy + 1.0/3.0*array([5.5,11.0,15.5]))


    def check_gauss_seidel(self):
        N = 1
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).T
        x = arange(N).astype(numpy.float64)
        b = zeros(N)
        gauss_seidel(A,x,b)
        assert_almost_equal(x,array([0]))

        N = 3
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).T
        x = arange(N).astype(numpy.float64)
        b = zeros(N)
        gauss_seidel(A,x,b)
        assert_almost_equal(x,array([1.0/2.0,5.0/4.0,5.0/8.0]))

        N = 1
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).T
        x = arange(N).astype(numpy.float64)
        b = zeros(N)
        gauss_seidel(A,x,b,sweep='backward')
        assert_almost_equal(x,array([0]))

        N = 3
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).T
        x = arange(N).astype(numpy.float64)
        b = zeros(N)
        gauss_seidel(A,x,b,sweep='backward')
        assert_almost_equal(x,array([1.0/8.0,1.0/4.0,1.0/2.0]))

        N = 1
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).T
        x = arange(N).astype(numpy.float64)
        b = array([10])
        gauss_seidel(A,x,b)
        assert_almost_equal(x,array([5]))

        N = 3
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).T
        x = arange(N).astype(numpy.float64)
        b = array([10,20,30])
        gauss_seidel(A,x,b)
        assert_almost_equal(x,array([11.0/2.0,55.0/4,175.0/8.0]))


        #forward and backward passes should give same result with x=ones(N),b=zeros(N)
        N = 100
        A = spdiags([2*ones(N),-ones(N),-ones(N)],[0,-1,1],N,N).T
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
    NumpyTest().run()
