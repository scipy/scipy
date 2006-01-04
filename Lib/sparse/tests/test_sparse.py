#!/usr/bin/env python
#
# Created by: Travis Oliphant, April 2004
#
""" Test functions for sparse matrices

"""
__usage__ = """
Build sparse:
  python setup.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.sparse.test(<level>)'
Run tests if sparse is not installed:
  python tests/test_sparse.py [<level>]
"""

import numpy
from numpy import arange, zeros, array, dot

import random
from numpy.testing import *
set_package_path()
from scipy.sparse import csc_matrix, csr_matrix, dok_matrix
restore_path()

class _test_cs(ScipyTestCase):

    def setUp(self):
        self.dat = array([[1,0,0,2],[3,0,1,0],[0,2,0,0]],'d')
        self.datsp = self.spmatrix(self.dat)

    def check_getelement(self):
        assert_equal(self.datsp[0,0],1)
        assert_equal(self.datsp[0,1],0)
        assert_equal(self.datsp[1,0],3)
        assert_equal(self.datsp[2,1],2)

    def check_todense(self):
        chk = self.datsp.todense()
        assert_array_equal(chk,self.dat)
        a = array([1.,2.,3.])
        dense_dot_dense = dot(a, self.dat)
        check = dot(a, self.datsp.todense())
        assert_array_equal(dense_dot_dense, check)
        b = array([1.,2.,3.,4.])
        dense_dot_dense = dot(self.dat, b)
        check2 = dot(self.datsp.todense(), b)
        assert_array_equal(dense_dot_dense, check2)

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
        b = self.spmatrix(array([[3,0,0],[0,1,0],[2,0,3.0],[2,3,0]]))
        assert_array_almost_equal(b*[1,2,3],dot(b.todense(),[1,2,3]))
        assert_array_almost_equal([1,2,3,4]*b,dot([1,2,3,4],b.todense()))

    def check_matmat(self):
        a = array([[3,0,0],[0,1,0],[2,0,3.0],[2,3,0]])
        b = array([[0,1],[1,0],[0,2]],'d')
        asp = self.spmatrix(a)
        bsp = self.spmatrix(b)
        assert_array_almost_equal((asp*bsp).todense(), dot(a,b))
        assert_array_almost_equal((asp*b).todense(), dot(a,b))
        assert_array_almost_equal((a*bsp).todense(), dot(a,b))
        # Now try performing cross-type multplication:
        csp = bsp.tocsc()
        c = b
        assert_array_almost_equal((asp*csp).todense(), dot(a,c))
        assert_array_almost_equal((asp.matmat(csp)).todense(), dot(a,c))
        assert_array_almost_equal((asp*c).todense(), dot(a,c))
        assert_array_almost_equal((a*csp).todense(), dot(a,c))
        csp = bsp.tocsr()
        assert_array_almost_equal((asp*csp).todense(), dot(a,c))
        assert_array_almost_equal((asp.matmat(csp)).todense(), dot(a,c))
        assert_array_almost_equal((asp*c).todense(), dot(a,c))
        assert_array_almost_equal((a*csp).todense(), dot(a,c))
        csp = bsp.tocoo()
        assert_array_almost_equal((asp*csp).todense(), dot(a,c))
        assert_array_almost_equal((asp.matmat(csp)).todense(), dot(a,c))
        assert_array_almost_equal((asp*c).todense(), dot(a,c))
        assert_array_almost_equal((a*csp).todense(), dot(a,c))

    def check_tocoo(self):
        a = self.datsp.tocoo()
        assert_array_almost_equal(a.todense(),self.dat)

    def check_tocsc(self):
        a = self.datsp.tocsc()
        assert_array_almost_equal(a.todense(),self.dat)

    def check_tocsr(self):
        a = self.datsp.tocsr()
        assert_array_almost_equal(a.todense(),self.dat)

    def check_transpose(self):
        a = self.datsp.transpose()
        b = self.dat.transpose()
        assert_array_equal(a.todense(), b)
        assert_array_equal(a.transpose().todense(), self.dat)
        assert_array_equal(a.transpose().todense(), self.datsp.todense())

    def check_large(self):
        # Create a 100x100 matrix with 100 non-zero elements
        # and play around with it
        A = dok_matrix()
        for k in range(100):
            i = random.randrange(100)
            j = random.randrange(100)
            A[i,j] = 1.
        csr = A.tocsr()
        csc = A.tocsc()
        csc2 = csr.tocsc()
        coo = A.tocoo()
        csr2 = coo.tocsr()
        assert_array_equal(A.transpose().todense(), csr.transpose().todense())
        assert_array_equal(csc.todense(), csr.todense())
        assert_array_equal(csr.todense(), csr2.todense())
        assert_array_equal(csr2.todense().transpose(), coo.transpose().todense())
        assert_array_equal(csr2.todense(), csc2.todense())
        csr_plus_csc = csr + csc
        csc_plus_csr = csc + csr
        assert_array_equal(csr_plus_csc.todense(), (2*A).todense())
        assert_array_equal(csr_plus_csc.todense(), csc_plus_csr.todense())

    def check_add_dense(self):
        """ Check whether adding a dense matrix to a sparse matrix works
        """
        sum1 = self.dat + self.datsp
        assert_array_equal(sum1, 2*self.dat)
        sum2 = self.datsp + self.dat
        assert_array_equal(sum2, 2*self.dat)

    def check_copy(self):
        """ Check whether the copy=True and copy=False keywords work
        """
        pass

    # Eventually we'd like to allow matrix products between dense
    # and sparse matrices using the normal dot() function:
    #def check_dense_dot_sparse(self):
    #    a = array([1.,2.,3.])
    #    dense_dot_dense = dot(a, self.dat)
    #    dense_dot_sparse = dot(a, self.datsp)
    #    assert_array_equal(dense_dot_dense, dense_dot_sparse)

    #def check_sparse_dot_dense(self):
    #    b = array([1.,2.,3.,4.])
    #    dense_dot_dense = dot(self.dat, b)
    #    dense_dot_sparse = dot(self.datsp, b)
    #    assert_array_equal(dense_dot_dense, dense_dot_sparse)

class test_csr(_test_cs):

    spmatrix = csr_matrix

    def check_constructor1(self):
        b = array([[0,4,0],
                   [3,0,1],
                   [0,2,0]],'d')
        bsp = csr_matrix(b)
        assert_array_almost_equal(bsp.data,[4,3,1,2])
        assert_array_equal(bsp.colind,[1,0,2,1])
        assert_array_equal(bsp.indptr,[0,1,3,4])
        assert_equal(bsp.getnnz(),4)
        assert_equal(bsp.getformat(),'csr')
        assert_array_almost_equal(bsp.todense(),b)

    def check_constructor2(self):
        b = zeros((6,6),'d')
        b[3,4] = 5
        bsp = csr_matrix(b)
        assert_array_almost_equal(bsp.data,[5])
        assert_array_equal(bsp.colind,[4])
        assert_array_equal(bsp.indptr,[0,0,0,0,1,1,1])
        assert_array_almost_equal(bsp.todense(),b)
        
    def check_constructor3(self):
        b = array([[1,0],
                   [0,2],
                   [3,0]],'d')
        bsp = csr_matrix(b)
        assert_array_almost_equal(bsp.data,[1,2,3])
        assert_array_equal(bsp.colind,[0,1,0])
        assert_array_equal(bsp.indptr,[0,1,2,3])
        assert_array_almost_equal(bsp.todense(),b)

class test_csc(_test_cs):

    spmatrix = csc_matrix

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

class test_dok(_test_cs):
    spmatrix = csr_matrix

    def check_mult(self):
        A = dok_matrix()
        A[0,3] = 10
        A[5,6] = 20
        D = A*A.T
        E = A*A.H
        assert_array_equal(D.A, E.A)
    
    def check_add(self):
        A = dok_matrix()
        A[0,1] = -10
        A[2,0] = 20
        A += 10
        B = array([[10, 0], [10, 10], [30, 10]])
        assert_array_equal(A.todense(), B)

    def check_convert(self):
        """Test provided by Andrew Straw.  Fails in SciPy <= r1477.
        """
        (m, n) = (6, 7)
        a=dok_matrix((m, n))
        
        # set a few elements, but none in the last column
        a[2,1]=1
        a[0,2]=2
        a[3,1]=3
        a[1,5]=4
        a[4,3]=5
        a[4,2]=6
        
        # assert that the last column is all zeros
        assert_array_equal( a.todense()[:,n-1], zeros(m,) )
        
        # make sure it still works for CSC format
        csc=a.tocsc()
        assert_array_equal( csc.todense()[:,n-1], zeros(m,) )

        # now test CSR
        (m, n) = (n, m)
        b = a.transpose()
        assert b.shape == (m, n)
        # assert that the last row is all zeros
        assert_array_equal( b.todense()[m-1,:], zeros(n,) )

        # make sure it still works for CSR format
        csr=b.tocsr()
        assert_array_equal( csr.todense()[m-1,:], zeros(n,))


if __name__ == "__main__":
    ScipyTest().run()
