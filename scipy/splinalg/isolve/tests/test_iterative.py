#!/usr/bin/env python
""" Test functions for the splinalg.isolve module
"""

from scipy.testing import *

from numpy import zeros, dot, diag, ones, arange, array
from numpy.random import rand
from scipy.linalg import norm
from scipy.sparse import spdiags

from scipy.splinalg.isolve import cg, cgs, bicg, bicgstab, gmres, qmr, minres

#def callback(x):
#    global A, b
#    res = b-dot(A,x)
#    #print "||A.x - b|| = " + str(norm(dot(A,x)-b))


#TODO check that method preserve shape and type
#TODO test complex matrices
#TODO test both preconditioner methods

data = ones((3,10))
data[0,:] =  2
data[1,:] = -1
data[2,:] = -1
Poisson1D = spdiags( data, [0,-1,1], 10, 10, format='csr')

data = array([[6, -5, 2, 7, -1, 10, 4, -3, -8, 9]],dtype='d')
RandDiag = spdiags( data, [0], 10, 10, format='csr' )

class TestIterative(TestCase):
    def setUp(self):
        # list of tuples (solver, symmetric, positive_definite )
        self.solvers = []
        self.solvers.append( (cg,       True,  True) )
        self.solvers.append( (cgs,      False, False) )
        self.solvers.append( (bicg,     False, False) )
        self.solvers.append( (bicgstab, False, False) )
        self.solvers.append( (gmres,    False, False) )
        self.solvers.append( (qmr,      False, False) )
        self.solvers.append( (minres,   True,  False) )
        
        # list of tuples (A, symmetric, positive_definite )
        self.cases = []

        # Symmetric and Positive Definite
        self.cases.append( (Poisson1D,True,True) )

        # Symmetric and Negative Definite
        self.cases.append( (-Poisson1D,True,False) )

        # Symmetric and Indefinite 
        self.cases.append( (RandDiag,True,False) )

        # Non-symmetric and Positive Definite
        # bicg and cgs fail to converge on this one
        #data = ones((2,10))
        #data[0,:] =  2
        #data[1,:] = -1
        #A = spdiags( data, [0,-1], 10, 10, format='csr')
        #self.cases.append( (A,False,True) )



    def test_convergence(self):
        """test whether all methods converge"""

        tol = 1e-8

        for solver,req_sym,req_pos in self.solvers:
            for A,sym,pos in self.cases:
                if req_sym and not sym: continue
                if req_pos and not pos: continue

                b  = arange(A.shape[0], dtype=float)
                x0 = 0*b

                x, info = solver(A, b, x0=x0, tol=tol)
                
                assert_array_equal(x0, 0*b) #ensure that x0 is not overwritten
                assert_equal(info,0)

                assert( norm(b - A*x) < tol*norm(b) )
    
    def test_precond(self):
        """test whether all methods accept a preconditioner"""

        tol = 1e-8

        for solver,req_sym,req_pos in self.solvers:
            for A,sym,pos in self.cases:
                if req_sym and not sym: continue
                if req_pos and not pos: continue

                M,N = A.shape
                D = spdiags( [abs(1.0/A.diagonal())], [0], M, N)
                def precond(b,which=None):
                    return D*b

                A = A.copy()
                A.psolve  = precond
                A.rpsolve = precond

                b  = arange(A.shape[0], dtype=float)
                x0 = 0*b

                x, info = solver(A, b, x0=x0, tol=tol)
                
                assert_equal(info,0)
                assert( norm(b - A*x) < tol*norm(b) )


class TestQMR(TestCase):
    def test_leftright_precond(self):
        """Check that QMR works with left and right preconditioners"""

        from scipy.splinalg.dsolve import splu
        from scipy.splinalg.interface import LinearOperator

        n = 100

        dat = ones(n)
        A = spdiags([-2*dat, 4*dat, -dat], [-1,0,1] ,n,n)
        b = arange(n,dtype='d')

        L = spdiags([-dat/2, dat], [-1,0], n, n)
        U = spdiags([4*dat, -dat], [ 0,1], n, n)
        
        L_solver = splu(L)
        U_solver = splu(U)

        def L_solve(b):
            return L_solver.solve(b)
        def U_solve(b):
            return U_solver.solve(b)
        def LT_solve(b):
            return L_solver.solve(b,'T')
        def UT_solve(b):
            return U_solver.solve(b,'T')

        M1 = LinearOperator( (n,n), matvec=L_solve, rmatvec=LT_solve )
        M2 = LinearOperator( (n,n), matvec=U_solve, rmatvec=UT_solve )

        x,info = qmr(A, b, tol=1e-8, maxiter=15, M1=M1, M2=M2)
       
        assert_equal(info,0)
        assert( norm(b - A*x) < 1e-8*norm(b) )


if __name__ == "__main__":
    nose.run(argv=['', __file__])
