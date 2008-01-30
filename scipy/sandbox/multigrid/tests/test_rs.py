from scipy.testing import *
from scipy import rand
from scipy.sparse import spdiags,csr_matrix,lil_matrix, \
                         isspmatrix_csr,isspmatrix_csc,isspmatrix_coo, \
                         isspmatrix_lil
import numpy

from scipy.sandbox.multigrid.rs import ruge_stuben_solver
from scipy.sandbox.multigrid.gallery import poisson

class TestRugeStubenSolver(TestCase):
    def test_poisson(self):
        cases = []
        
        cases.append( (1000,) )
        cases.append( (500,500) )
        cases.append( (50,50,50) )

        for case in cases:
            A = poisson( case, format='csr' )

            numpy.random.seed(0) #make tests repeatable

            x = rand(A.shape[0])
            b = A*rand(A.shape[0]) #zeros_like(x)

            ml = ruge_stuben_solver(A)

            x_sol,residuals = ml.solve(b,x0=x,maxiter=20,tol=1e-12,return_residuals=True)

            avg_convergence_ratio = (residuals[-1]/residuals[0])**(1.0/len(residuals))
            
            print "avg",avg_convergence_ratio

            assert(avg_convergence_ratio < 0.15)

if __name__ == '__main__':
    nose.run(argv=['', __file__])
