from numpy.linalg import norm
from numpy import zeros_like
import scipy
import numpy

from multigrid import sa_interpolation,rs_interpolation
from relaxation import gauss_seidel,jacobi 



def poisson_problem1D(N):
    """
    Return a sparse CSC matrix for the 1d poisson problem
    with standard 3-point finite difference stencil on a
    grid with N points.
    """
    D = 2*numpy.ones(N)
    O =  -numpy.ones(N)
    return scipy.sparse.spdiags([D,O,O],[0,-1,1],N,N)

def poisson_problem2D(N):
    """
    Return a sparse CSC matrix for the 2d poisson problem
    with standard 5-point finite difference stencil on a
    square N-by-N grid.
    """
    D = 4*numpy.ones(N*N)
    T =  -numpy.ones(N*N)
    O =  -numpy.ones(N*N)
    T[N-1::N] = 0
    return scipy.sparse.spdiags([D,O,T,T,O],[0,-N,-1,1,N],N*N,N*N)

def ruge_stuben_solver(A,max_levels=10,max_coarse=500):
    As = [A]
    Ps = []
    
    while len(As) < max_levels  and A.shape[0] > max_coarse:
        P = rs_interpolation(A)
        
        A = (P.T.tocsr() * A) * P     #galerkin operator

        As.append(A)
        Ps.append(P)
        
    return multilevel_solver(As,Ps)

def smoothed_aggregation_solver(A,max_levels=10,max_coarse=500):
    As = [A]
    Ps = []
    
    while len(As) < max_levels  and A.shape[0] > max_coarse:
        P = sa_interpolation(A)
        
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
        residuals = [norm(b-A*x,2)]

        while len(residuals) <= maxiter and residuals[-1]/residuals[0] > tol:
            self.__solve(0,x,b)

            residuals.append(scipy.linalg.norm(b-A*x,2))

            if callback is not None:
                callback(x)

        if return_residuals:
            return x,residuals
        else:
            return x
        
        
    def __solve(self,lvl,x,b):
        A = self.As[lvl]
        
        if len(self.As) == 1:
            x[:] = scipy.linalg.solve(A.todense(),b)
            return x

        self.presmoother(A,x,b)
            
        residual = b - A*x                                

        coarse_x = zeros((self.As[lvl+1].shape[0]))
        coarse_b = self.Ps[lvl].T * residual
        
        if lvl == len(self.As) - 2:
            coarse_x[:] = scipy.linalg.solve(self.As[-1].todense(),coarse_b)
        else:   
            self.__solve(lvl+1,coarse_x,coarse_b)
                
        x += self.Ps[lvl] * coarse_x   #coarse grid correction

        self.postsmoother(A,x,b)


    def presmoother(self,A,x,b):
        gauss_seidel(A,x,b,iterations=1,sweep="forward")
    
    def postsmoother(self,A,x,b):
        gauss_seidel(A,x,b,iterations=1,sweep="backward")



if __name__ == '__main__':
    from scipy import *
    A = poisson_problem2D(100).T
    asa = smoothed_aggregation_solver(A)
    #asa = ruge_stuben_solver(A)
    x = rand(A.shape[0])
    b = zeros_like(x)
    
    resid = []
    
    for n in range(10):
        x = asa.solve(b,x,maxiter=1)
        resid.append(linalg.norm(A*x))




