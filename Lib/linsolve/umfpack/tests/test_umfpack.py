#!/usr/bin/env python
#

""" Test functions for UMFPACK wrappers

"""

from numpy import transpose, array, arange

import random
from numpy.testing import *
set_package_path()
from scipy import linsolve, rand, matrix, diag, eye
from scipy.sparse import csc_matrix, dok_matrix, spdiags

import numpy as nm
import scipy.linsolve.umfpack as um

restore_path()

class test_solvers(ScipyTestCase):
    """Tests inverting a sparse linear system"""
    
    def check_solve_complex_without_umfpack(self):
        """Solve: single precision complex"""
        linsolve.use_solver( {'useUmfpack' :  False} )
        a = self.a.astype('F')
        b = self.b
        x = linsolve.spsolve(a, b)
        #print x
        #print "Error: ", a*x-b
        assert_array_almost_equal(a*x, b)
        
        
    def check_solve_without_umfpack(self): 
        """Solve: single precision"""
        linsolve.use_solver( {'useUmfpack' :  False} )
        a = self.a.astype('f')
        b = self.b
        x = linsolve.spsolve(a, b.astype('f'))
        #print x
        #print "Error: ", a*x-b
        assert_array_almost_equal(a*x, b)


    def check_solve_complex_umfpack(self):
        """Solve with UMFPACK: double precision complex"""
        linsolve.use_solver( {'useUmfpack' :  True} )
        a = self.a.astype('D')
        b = self.b
        x = linsolve.spsolve(a, b)
        #print x
        #print "Error: ", a*x-b
        assert_array_almost_equal(a*x, b)

    def check_solve_umfpack(self):
        """Solve with UMFPACK: double precision"""
        linsolve.use_solver( {'useUmfpack' :  True} )
        a = self.a.astype('d')
        b = self.b
        x = linsolve.spsolve(a, b)
        #print x
        #print "Error: ", a*x-b
        assert_array_almost_equal(a*x, b)


    def setUp(self):
        self.a = spdiags([[1, 2, 3, 4, 5], [6, 5, 8, 9, 10]], [0, 1], 5, 5)
        #print "The sparse matrix (constructed from diagonals):"
        #print self.a
        self.b = array([1, 2, 3, 4, 5])


        
class test_factorization(ScipyTestCase):
    """Tests factorizing a sparse linear system"""
    
    def check_complex_lu(self):
        """Getting factors of complex matrix"""
        umfpack = um.UmfpackContext("zi")

        for A in self.complex_matrices:                             
            umfpack.numeric(A)
            
            (L,U,P,Q,R,do_recip) = umfpack.lu(A)
            
            L = L.todense()
            U = U.todense()
            A = A.todense()
            if not do_recip: R = 1.0/R
            R = matrix(diag(R))
            P = eye(A.shape[0])[P,:]
            Q = eye(A.shape[1])[:,Q]
            
            assert_array_almost_equal(P*R*A*Q,L*U)

    def check_real_lu(self):
        """Getting factors of real matrix"""
        umfpack = um.UmfpackContext("di")

        for A in self.real_matrices:           
            umfpack.numeric(A)
            
            (L,U,P,Q,R,do_recip) = umfpack.lu(A)
            
            L = L.todense()
            U = U.todense()
            A = A.todense()
            if not do_recip: R = 1.0/R
            R = matrix(diag(R))
            P = eye(A.shape[0])[P,:]
            Q = eye(A.shape[1])[:,Q]
            
            assert_array_almost_equal(P*R*A*Q,L*U)
            

    def setUp(self):     
        random.seed(0) #make tests repeatable     
        self.real_matrices = []        
        self.real_matrices.append(spdiags([[1, 2, 3, 4, 5], [6, 5, 8, 9, 10]],
                                          [0, 1], 5, 5))
        self.real_matrices.append(spdiags([[1, 2, 3, 4, 5], [6, 5, 8, 9, 10]],
                                          [0, 1], 4, 5))
        self.real_matrices.append(spdiags([[1, 2, 3, 4, 5], [6, 5, 8, 9, 10]],
                                          [0, 2], 5, 5))
        self.real_matrices.append(csc_matrix(rand(3,3)))
        self.real_matrices.append(csc_matrix(rand(5,4)))
        self.real_matrices.append(csc_matrix(rand(4,5)))
        
        self.complex_matrices = [x.astype(nm.complex128)
                                 for x in self.real_matrices]
        
        
if __name__ == "__main__":
    ScipyTest().run()
