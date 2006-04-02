#!/usr/bin/env python
#

""" Test functions for UMFPACK wrappers

"""

#import numpy
#from numpy import arange, zeros, array, dot, ones, matrix, asmatrix, asarray, \
#        float32, float64, complex64, complex128
from numpy import transpose, array, arange

import random
from numpy.testing import *
set_package_path()
from scipy import linsolve
from scipy.sparse import csc_matrix, dok_matrix, spdiags
restore_path()

class test_solvers(ScipyTestCase):
    """Tests inverting a sparse linear system"""
    
    # def check_solve_complex_without_umfpack(self):
    #     """Solve: single precision complex"""
    #     linsolve.useUmfpack = False   # What should this be? This doesn't help!
    #     a = self.a.astype('F')
    #     b = self.b
    #     x = linsolve.spsolve(a, b)
    #     #print x
    #     #print "Error: ", a*x-b
    #     assert_array_almost_equal(a*x, b)
    #     
    #     
    # def check_solve_without_umfpack(self): 
    #     """Solve: single precision"""
    #     linsolve.useUmfpack = False   # What should this be? This doesn't help!
    #     a = self.a.astype('f')
    #     b = self.b
    #     x = linsolve.spsolve(a, b.astype('f'))
    #     #print x
    #     #print "Error: ", a*x-b
    #     assert_array_almost_equal(a*x, b)


    def check_solve_complex_umfpack(self):
        """Solve with UMFPACK: double precision complex"""
        # globals()['useUmfpack'] = True
        linsolve.useUmfpack = True
        a = self.a.astype('D')
        b = self.b
        x = linsolve.spsolve(a, b)
        print x
        #print "Error: ", a*x-b
        assert_array_almost_equal(a*x, b)

    def check_solve_umfpack(self):
        """Solve with UMFPACK: double precision"""
        # globals()['useUmfpack'] = True
        linsolve.useUmfpack = True
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


if __name__ == "__main__":
    ScipyTest().run()
