__all__ = ['poisson_problem1D','poisson_problem2D',
           'ruge_stuben_solver','smoothed_aggregation_solver',
           'multilevel_solver']

from numpy.linalg import norm
from numpy import zeros,zeros_like,array
import scipy
import numpy

from coarsen import sa_interpolation,rs_interpolation
from relaxation import gauss_seidel,jacobi,sor
from utils import infinity_norm


def poisson_problem1D(N):
    """
    Return a sparse CSR matrix for the 1d poisson problem
    with standard 3-point finite difference stencil on a
    grid with N points.
    """
    D = 2*numpy.ones(N)
    O =  -numpy.ones(N)
    return scipy.sparse.spdiags([D,O,O],[0,-1,1],N,N).tocoo().tocsr() #eliminate zeros

def poisson_problem2D(N):
    """
    Return a sparse CSR matrix for the 2d poisson problem
    with standard 5-point finite difference stencil on a
    square N-by-N grid.
    """
    D = 4*numpy.ones(N*N)
    T =  -numpy.ones(N*N)
    O =  -numpy.ones(N*N)
    T[N-1::N] = 0
    return scipy.sparse.spdiags([D,O,T,T,O],[0,-N,-1,1,N],N*N,N*N).tocoo().tocsr() #eliminate zeros
    

def ruge_stuben_solver(A,max_levels=10,max_coarse=500):
    """
    Create a multilevel solver using Ruge-Stuben coarsening (Classical AMG)
    
        References:
            "Multigrid"
                Trottenberg, U., C. W. Oosterlee, and Anton Schuller. San Diego: Academic Press, 2001.
                See Appendix A
    
    """
    As = [A]
    Ps = []
    
    while len(As) < max_levels  and A.shape[0] > max_coarse:
        P = rs_interpolation(A)
        
        A = (P.T.tocsr() * A) * P     #galerkin operator

        As.append(A)
        Ps.append(P)
        
    return multilevel_solver(As,Ps)

def smoothed_aggregation_solver(A,candidates=None,blocks=None,max_levels=10,max_coarse=500,epsilon=0.08,omega=4.0/3.0):
    """
    Create a multilevel solver using Smoothed Aggregation (SA)

        References:
            "Algebraic Multigrid by Smoothed Aggregation for Second and Fourth Order Elliptic Problems",
                Petr Vanek and Jan Mandel and Marian Brezina
                http://citeseer.ist.psu.edu/vanek96algebraic.html
    
    """
    As = [A]
    Ps = []
    
    if candidates is None:
        candidates = [ ones(A.shape[0]) ] # use constant vector
        
    while len(As) < max_levels  and A.shape[0] > max_coarse:
        P,candidates = sa_interpolation(A,candidates,epsilon*0.5**(len(As)-1),omega=omega,blocks=blocks)
        #P = sa_interpolation(A,epsilon=0.0)

        A = (P.T.tocsr() * A) * P     #galerkin operator

        As.append(A)
        Ps.append(P)
        
    return multilevel_solver(As,Ps)


class multilevel_solver:
    def __init__(self,As,Ps):
        self.As = As
        self.Ps = Ps
            
    def __repr__(self):
        output = 'multilevel_solver\n'
        output += 'Number of Levels:     %d\n' % len(self.As)
        output += 'Operator Complexity: %6.3f\n' % self.operator_complexity()
        output += 'Grid Complexity:     %6.3f\n' % self.grid_complexity()

        total_nnz =  sum([A.nnz for A in self.As])
        
        for n,A in enumerate(self.As):
            output += '   [level %2d]  unknowns: %10d  nnz: %5.2f%%\n' % (n,A.shape[1],(100*float(A.nnz)/float(total_nnz)))

        return output

    def operator_complexity(self):
        """number of nonzeros on all levels / number of nonzeros on the finest level"""
        return sum([A.nnz for A in self.As])/float(self.As[0].nnz)
    
    def grid_complexity(self):
        """number of unknowns on all levels / number of unknowns on the finest level"""
        return sum([A.shape[0] for A in self.As])/float(self.As[0].shape[0])
    

    def psolve(self, b):
        return self.solve(b,maxiter=1)

    def solve(self, b, x0=None, tol=1e-5, maxiter=100, callback=None, return_residuals=False):
        """
        TODO
        """

        if x0 is None:
            x = zeros_like(b)
        else:
            x = array(x0)


        #TODO change use of tol (relative tolerance) to agree with other iterative solvers
        A = self.As[0]
        residuals = [scipy.linalg.norm(b-A*x)]

        while len(residuals) <= maxiter and residuals[-1]/residuals[0] > tol:
            self.__solve(0,x,b)

            residuals.append(scipy.linalg.norm(b-A*x))

            if callback is not None:
                callback(x)

        if return_residuals:
            return x,residuals
        else:
            return x
        
        
    def __solve(self,lvl,x,b):
        A = self.As[lvl]
        
        if len(self.As) == 1:
            x[:] = scipy.linsolve.spsolve(A,b)
            return 

        self.presmoother(A,x,b)

        residual = b - A*x                                

        coarse_x = zeros((self.As[lvl+1].shape[0]))
        coarse_b = self.Ps[lvl].T * residual
        
        if lvl == len(self.As) - 2:
            #use direct solver on coarsest level
            coarse_x[:] = scipy.linsolve.spsolve(self.As[-1],coarse_b)
            #coarse_x[:] = scipy.linalg.cg(self.As[-1],coarse_b,tol=1e-12)[0]
            #print "coarse residual norm",scipy.linalg.norm(coarse_b - self.As[-1]*coarse_x)
        else:   
            self.__solve(lvl+1,coarse_x,coarse_b)
                
        x += self.Ps[lvl] * coarse_x   #coarse grid correction

        self.postsmoother(A,x,b)


    def presmoother(self,A,x,b):
        gauss_seidel(A,x,b,iterations=1,sweep="symmetric")
    
    def postsmoother(self,A,x,b):
        gauss_seidel(A,x,b,iterations=1,sweep="symmetric")



if __name__ == '__main__':
    from scipy import *
    candidates = None
    A = poisson_problem2D(100)
    #A = io.mmread("rocker_arm_surface.mtx").tocsr()
    #A = io.mmread("9pt-100x100.mtx").tocsr()
    #A = io.mmread("/home/nathan/Desktop/9pt/9pt-100x100.mtx").tocsr()
    #A = io.mmread("/home/nathan/Desktop/BasisShift_W_EnergyMin_Luke/9pt-5x5.mtx").tocsr()
    #A = io.mmread('tests/sample_data/elas30_A.mtx').tocsr()
    #candidates = io.mmread('tests/sample_data/elas30_nullspace.mtx')
    #candidates = [ array(candidates[:,x]) for x in range(candidates.shape[1]) ]
    
    ml = smoothed_aggregation_solver(A,candidates,max_coarse=100,max_levels=3)
    #ml = ruge_stuben_solver(A)

    x = rand(A.shape[0])
    b = zeros_like(x)
    #b = rand(A.shape[0])
    
    if True:
        x_sol,residuals = ml.solve(b,x0=x,maxiter=30,tol=1e-12,return_residuals=True)
    else:
        residuals = []
        def add_resid(x):
            residuals.append(linalg.norm(b - A*x))
        A.psolve = ml.psolve
        x_sol = linalg.cg(A,b,x0=x,maxiter=12,tol=1e-100,callback=add_resid)[0]

    residuals = array(residuals)/residuals[0]

    print residuals



