#!/usr/bin/env python
#
# Created by: Travis Oliphant, April 2004
#
""" Test functions for Sparse matrices

"""
__usage__ = """
Build linalg:
  python setup_linalg.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.sparse.test(<level>)'
Run tests if linalg is not installed:
  python tests/test_Sparse.py [<level>]
"""

import scipy_base
from scipy_base import arange, zeros, array, dot

import sys
from scipy_test.testing import *
from scipy.sparse import csc_matrix

import unittest


class test_csc(ScipyTestCase):

    dat = array([[1,0,0,2],[3,0,1,0],[0,2,0,0]],'d')
    datsp = csc_matrix(dat)

    def check_constructor1(self):
        b = array([[1,0,0],[3,0,1],[0,2,0]],'d')
        bsp = csc_matrix(b)
        assert_array_almost_equal(bsp.data,[1,3,2,1])
        assert_array_equal(bsp.rowind,[0,1,2,1])
        assert_array_equal(bsp.indptr,[0,2,3,4])
        assert_equal(bsp.getnnz(),4)
        assert_equal(bsp.getformat(),'csc')
        
    def check_constructor2(self):
        b = zeros((6,6),'d')
        b[2,4] = 5
        bsp = csc_matrix(b)
        assert_array_almost_equal(bsp.data,[5])
        assert_array_equal(bsp.rowind,[2])
        assert_array_equal(bsp.indptr,[0,0,0,0,0,1,1])

    def check_constructor3(self):
        b = array([[1,0],[0,2],[3,0]],'d')
        bsp = csc_matrix(b)
        assert_array_almost_equal(bsp.data,[1,3,2])
        assert_array_equal(bsp.rowind,[0,2,1])
        assert_array_equal(bsp.indptr,[0,2,3])

    def check_getelement(self):
        assert_equal(self.datsp[0,0],1)
        assert_equal(self.datsp[0,1],0)
        assert_equal(self.datsp[1,0],3)
        assert_equal(self.datsp[2,1],2)

    def check_todense(self):
        chk = self.datsp.todense()
        assert_array_equal(chk,self.dat)

    def check_setelement(self):
        a = self.datsp - self.datsp
        a[1,2] = 4.0
        a[0,1] = 3
        a[2,0] = 2.0
        assert_array_equal(a.todense(),[[0,3,0,0],[0,0,4,0],[2,0,0,0]])

    def check_add(self):
        a = self.datsp
        b = self.datsp.copy()
        b[0,2] = 2.0
        c = a + b
        assert_array_equal(c.todense(),[[2,0,2,4],[6,0,2,0],[0,4,0,0]])

    def check_elmul(self):
        a = self.datsp
        b = self.datsp.copy()
        b[0,2] = 2.0
        c = a ** b
        assert_array_equal(c.todense(),[[1,0,0,4],[9,0,1,0],[0,4,0,0]])

    def check_matvec(self):
        b = csc_matrix(array([[3,0,0],[0,1,0],[2,0,3.0],[2,3,0]]))
        assert_array_almost_equal(b*[1,2,3],dot(b.todense(),[1,2,3]))
        assert_array_almost_equal([1,2,3,4]*b,dot([1,2,3,4],b.todense()))

    def check_matmat(self):
        asp = csc_matrix(array([[3,0,0],[0,1,0],[2,0,3.0],[2,3,0]]))
        bsp = csc_matrix(array([[0,1],[1,0],[0,2]],'d'))
        assert_array_almost_equal((asp*bsp).todense(),dot(asp.todense(),bsp.todense()))
    
if __name__ == "__main__":
    ScipyTest('sparse.Sparse').run()
