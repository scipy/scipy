__all__ = ['poisson_problem1D','poisson_problem2D',
           'ruge_stuben_solver','smoothed_aggregation_solver',
           'multilevel_solver']

import scipy
import numpy
from numpy import ones,zeros,zeros_like,array,asarray,empty
from numpy.linalg import norm
from scipy.linsolve import spsolve
from scipy.sparse import dia_matrix

from sa import sa_interpolation
from rs import rs_interpolation
from relaxation import gauss_seidel,jacobi,sor
from utils import symmetric_rescaling, diag_sparse


def poisson_problem1D(N):
    """
    Return a sparse CSR matrix for the 1d poisson problem
    with standard 3-point finite difference stencil on a
    grid with N points.
    """
    D = 2*ones(N)
    O =  -ones(N)
    return scipy.sparse.spdiags([D,O,O],[0,-1,1],N,N).tocoo().tocsr() #eliminate explicit zeros

def poisson_problem2D(N, epsilon=1.0, dtype='d', format=None):
    """
    Return a sparse matrix for the 2d poisson problem
    with standard 5-point finite difference stencil on a
    square N-by-N grid.
    """
    if N == 1:
        diags   = asarray( [[4]],dtype=dtype)
        return dia_matrix((diags,[0]), shape=(1,1)).asformat(format)

    offsets = array([0,-N,N,-1,1])

    diags = empty((5,N**2),dtype=dtype)

    diags[0]  =  (2 + 2*epsilon) #main diagonal
    diags[1,:] = -1
    diags[2,:] = -1
    
    diags[3,:] = -epsilon #first lower diagonal 
    diags[4,:] = -epsilon #first upper diagonal 
    diags[3,N-1::N] = 0  
    diags[4,N::N]   = 0    
    
    return dia_matrix((diags,offsets),shape=(N**2,N**2)).tocoo().tocsr() #eliminate explicit zeros
    #return dia_matrix((diags,offsets),shape=(N**2,N**2)).asformat(format)


def ruge_stuben_solver(A,max_levels=10,max_coarse=500):
    """
    Create a multilevel solver using Ruge-Stuben coarsening (Classical AMG)

        References:
            "Multigrid"
                Trottenberg, U., C. W. Oosterlee, and Anton Schuller.
                San Diego: Academic Press, 2001.
                Appendix A

    """
    As = [A]
    Ps = []

    while len(As) < max_levels  and A.shape[0] > max_coarse:
        P = rs_interpolation(A)

        A = (P.T.tocsr() * A) * P     #galerkin operator

        As.append(A)
        Ps.append(P)

    return multilevel_solver(As,Ps)



def smoothed_aggregation_solver(A, B=None,  \
                                aggregation=None, max_levels=10, \
                                max_coarse=500, epsilon=0.0, \
                                omega=4.0/3.0, symmetric=True, \
                                rescale = True):
    """Create a multilevel solver using Smoothed Aggregation (SA)

    *Parameters*:

        A : {csr_matrix}
            NxN matrix in CSR or BSR format
        B : {None, array_like} : optional
            Near-nullspace candidates stored in the columns of an NxK array.
            The default value B=None is equivalent to B=ones((N,1))
        blocks : {None, array_like} : optional
            Array of length N that groups the variables into 'superblocks'.
            For example, in a 2d vector-valued problem where the even
            variables [0,2,4,...N-2] correspond to the x-components and the
            odd variables [1,3,5,...,N-1] correspond to the y-components then
            blocks=[0,0,1,1,2,2,...,N/2,N/2] is expected.  The default
            value blocks=None is equivalent to blocks=[0,1,2,..,N] which
            implies that each variable should be aggregated seperately.
            The default is appropriate for scalar valued problems.
        aggregation: {None, list of csr_matrix} : optional
            List of csr_matrix objects that describe a user-defined
            multilevel aggregation of the variables.
            TODO ELABORATE
        max_levels: {integer} : default 10
            Maximum number of levels to be used in the multilevel solver.
        max_coarse: {integer} : default 500
            Maximum number of variables permitted on the coarse grid.
        epsilon: {float} : default 0.0
            Strength of connection parameter used in aggregation.
        omega: {float} : default 4.0/3.0
            Damping parameter used in prolongator smoothing (0 < omega < 2)
        symmetric: {boolean} : default True
            True if A is symmetric, False otherwise
        rescale: {boolean} : default True
            If True, symmetrically rescale A by the diagonal
            i.e. A -> D * A * D,  where D is diag(A)^-0.5

    *Example*:
        TODO

    *References*:
        "Algebraic Multigrid by Smoothed Aggregation for Second and Fourth Order Elliptic Problems",
            Petr Vanek and Jan Mandel and Marian Brezina
            http://citeseer.ist.psu.edu/vanek96algebraic.html

    """

    A = A.asfptype()

    if B is None:
        B = ones((A.shape[0],1),dtype=A.dtype) # use constant vector
    else:
        B = asarray(B,dtype=A.dtype)

    pre,post = None,None   #preprocess/postprocess

    if rescale:
        D_sqrt,D_sqrt_inv,A = symmetric_rescaling(A)
        D_sqrt,D_sqrt_inv = diag_sparse(D_sqrt),diag_sparse(D_sqrt_inv)

        B = D_sqrt * B  #scale candidates
        def pre(x,b):
            return D_sqrt*x,D_sqrt_inv*b
        def post(x):
            return D_sqrt_inv*x

    As = [A]
    Ps = []

    if aggregation is None:
        while len(As) < max_levels and A.shape[0] > max_coarse:
            P,B = sa_interpolation(A,B,epsilon*0.5**(len(As)-1),omega=omega)

            A = (P.T.asformat(P.format) * A) * P     #galerkin operator

            As.append(A)
            Ps.append(P)
    else:
        #use user-defined aggregation
        for AggOp in aggregation:
            P,B = sa_interpolation(A,B,omega=omega,AggOp=AggOp)

            A = (P.T.tocsr() * A) * P     #galerkin operator

            As.append(A)
            Ps.append(P)

    return multilevel_solver(As,Ps,preprocess=pre,postprocess=post)


class multilevel_solver:
    def __init__(self,As,Ps,preprocess=None,postprocess=None):
        self.As = As
        self.Ps = Ps
        self.preprocess = preprocess
        self.postprocess = postprocess

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
            x = array(x0) #copy

        if self.preprocess is not None:
            x,b = self.preprocess(x,b)

        #TODO change use of tol (relative tolerance) to agree with other iterative solvers
        A = self.As[0]
        residuals = [ norm(b-A*x) ]

        while len(residuals) <= maxiter and residuals[-1]/residuals[0] > tol:
            self.__solve(0,x,b)

            residuals.append( norm(b-A*x) )

            if callback is not None:
                callback(x)

        if self.postprocess is not None:
            x = self.postprocess(x)

        if return_residuals:
            return x,residuals
        else:
            return x


    def __solve(self,lvl,x,b):
        A = self.As[lvl]

        if len(self.As) == 1:
            #TODO make spsolve preserve dimensions
            x[:] = spsolve(A,b).reshape(x.shape)
            return

        self.presmoother(A,x,b)

        residual = b - A*x

        coarse_b = self.Ps[lvl].T * residual
        coarse_x = zeros_like(coarse_b)

        if lvl == len(self.As) - 2:
            #use direct solver on coarsest level
            #TODO reuse factors for efficiency?
            coarse_x[:] = spsolve(self.As[-1],coarse_b).reshape(coarse_x.shape)
            #coarse_x[:] = scipy.linalg.cg(self.As[-1],coarse_b,tol=1e-12)[0].reshape(coarse_x.shape)
            #A_inv = asarray(scipy.linalg.pinv2(self.As[-1].todense()))
            #coarse_x[:] = scipy.dot(A_inv,coarse_b)
            #print "coarse residual norm",scipy.linalg.norm(coarse_b - self.As[-1]*coarse_x)
        else:
            self.__solve(lvl+1,coarse_x,coarse_b)

        x += self.Ps[lvl] * coarse_x   #coarse grid correction

        self.postsmoother(A,x,b)


    def presmoother(self,A,x,b):
        gauss_seidel(A,x,b,iterations=1,sweep="symmetric")

    def postsmoother(self,A,x,b):
        gauss_seidel(A,x,b,iterations=1,sweep="symmetric")


