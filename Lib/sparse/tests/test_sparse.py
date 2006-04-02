#!/usr/bin/env python
#
# Authors: Travis Oliphant, Ed Schofield, Robert Cimrman, and others

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
from numpy import arange, zeros, array, dot, ones, matrix, asmatrix, asarray, \
        float32, float64, complex64, complex128

import random
from numpy.testing import *
set_package_path()
from scipy.sparse import csc_matrix, csr_matrix, dok_matrix, spidentity, \
        speye, lil_matrix
from scipy.linsolve import splu
restore_path()

class _test_cs(ScipyTestCase):

    def setUp(self):
        self.dat = matrix([[1,0,0,2],[3,0,1,0],[0,2,0,0]],'d')
        self.datsp = self.spmatrix(self.dat)

    def check_getelement(self):
        assert_equal(self.datsp[0,0],1)
        assert_equal(self.datsp[0,1],0)
        assert_equal(self.datsp[1,0],3)
        assert_equal(self.datsp[2,1],2)

    def check_sum(self):
        """Does the matrix's sum() method work?
        """
        assert_array_equal(self.dat.sum(), self.datsp.sum())
        assert_array_equal(self.dat.sum(axis=None), self.datsp.sum(axis=None))
        assert_array_equal(self.dat.sum(axis=0), self.datsp.sum(axis=0))
        assert_array_equal(self.dat.sum(axis=1), self.datsp.sum(axis=1))

    def check_todense(self):
        chk = self.datsp.todense()
        assert_array_equal(chk,self.dat)
        a = matrix([1.,2.,3.])
        dense_dot_dense = a * self.dat
        check = a * self.datsp.todense()
        assert_array_equal(dense_dot_dense, check)
        b = matrix([1.,2.,3.,4.]).T
        dense_dot_dense = self.dat * b
        check2 = self.datsp.todense() * b
        assert_array_equal(dense_dot_dense, check2)

    def check_toarray(self):
        dat = asarray(self.dat)
        chk = self.datsp.toarray()
        assert_array_equal(chk, dat)
        a = array([1.,2.,3.])
        dense_dot_dense = dot(a, dat)
        check = dot(a, self.datsp.toarray())
        assert_array_equal(dense_dot_dense, check)
        b = array([1.,2.,3.,4.])
        dense_dot_dense = dot(dat, b)
        check2 = dot(self.datsp.toarray(), b)
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

    def check_rmatvec(self):
        M = self.spmatrix(matrix([[3,0,0],[0,1,0],[2,0,3.0],[2,3,0]]))
        assert_array_almost_equal([1,2,3,4]*M, dot([1,2,3,4], M.toarray()))
        row = matrix([[1,2,3,4]])
        assert_array_almost_equal(row*M, row*M.todense())

    def check_matvec(self):
        M = self.spmatrix(matrix([[3,0,0],[0,1,0],[2,0,3.0],[2,3,0]]))
        col = matrix([1,2,3]).T
        assert_array_almost_equal(M * col, M.todense() * col)

        # Should this be supported or not?!
        #flat = array([1,2,3])
        #assert_array_almost_equal(M*flat, M.todense()*flat)
        # Currently numpy dense matrices promote the result to a 1x3 matrix,
        # whereas sparse matrices leave the result as a rank-1 array.  Which
        # is preferable?

        # Note: the following command does not work.  Both NumPy matrices
        # and spmatrices should raise exceptions!
        # assert_array_almost_equal(M*[1,2,3], M.todense()*[1,2,3])

        # The current relationship between sparse matrix products and array
        # products is as follows:
        assert_array_almost_equal(M*array([1,2,3]), dot(M.A,[1,2,3]))
        assert_array_almost_equal(M*[[1],[2],[3]], asmatrix(dot(M.A,[1,2,3])).T)
        # Note that the result of M * x is dense if x has a singleton dimension.

        # Currently M.matvec(asarray(col)) is rank-1, whereas M.matvec(col)
        # is rank-2.  Is this desirable?

    def check_matmat(self):
        a = matrix([[3,0,0],[0,1,0],[2,0,3.0],[2,3,0]])
        a2 = array([[3,0,0],[0,1,0],[2,0,3.0],[2,3,0]])
        b = matrix([[0,1],[1,0],[0,2]],'d')
        asp = self.spmatrix(a)
        bsp = self.spmatrix(b)
        assert_array_almost_equal((asp*bsp).todense(), a*b)
        assert_array_almost_equal((asp*b).todense(), a*b)
        assert_array_almost_equal((a*bsp).todense(), a*b)
        assert_array_almost_equal((a2*bsp).todense(), a*b)

        # Now try performing cross-type multplication:
        csp = bsp.tocsc()
        c = b
        assert_array_almost_equal((asp*csp).todense(), a*c)
        assert_array_almost_equal((asp.matmat(csp)).todense(), a*c)
        assert_array_almost_equal((asp*c).todense(), a*c)
        
        assert_array_almost_equal((a*csp).todense(), a*c)
        assert_array_almost_equal((a2*csp).todense(), a*c)
        csp = bsp.tocsr()
        assert_array_almost_equal((asp*csp).todense(), a*c)
        assert_array_almost_equal((asp.matmat(csp)).todense(), a*c)
        assert_array_almost_equal((asp*c).todense(), a*c)

        assert_array_almost_equal((a*csp).todense(), a*c)
        assert_array_almost_equal((a2*csp).todense(), a*c)
        csp = bsp.tocoo()
        assert_array_almost_equal((asp*csp).todense(), a*c)
        assert_array_almost_equal((asp.matmat(csp)).todense(), a*c)
        assert_array_almost_equal((asp*c).todense(), a*c)

        assert_array_almost_equal((a*csp).todense(), a*c)
        assert_array_almost_equal((a2*csp).todense(), a*c)

        # Test provided by Andy Fraser, 2006-03-26
        L = 30
        frac = .3
        random.seed(0) # make runs repeatable
        A = self.spmatrix((L,2))
        for i in xrange(L):
            for j in xrange(2):
                r = random.random()
                if r < frac:
                    A[i,j] = r/frac
        B = A*A.T
        assert_array_almost_equal(B.todense(), A.todense() * A.T.todense())
        assert_array_almost_equal(B.todense(), A.todense() * A.todense().T)
    
    
    def check_tocoo(self):
        a = self.datsp.tocoo()
        assert_array_almost_equal(a.todense(), self.dat)

    def check_tocsc(self):
        a = self.datsp.tocsc()
        assert_array_almost_equal(a.todense(), self.dat)
        b = complexsp = self.spmatrix(self.dat+3j)
        c = b.tocsc()
        assert_array_almost_equal(c.todense(), self.dat+3j)

    def check_tocsr(self):
        a = self.datsp.tocsr()
        assert_array_almost_equal(a.todense(), self.dat)

    def check_transpose(self):
        a = self.datsp.transpose()
        b = self.dat.transpose()
        assert_array_equal(a.todense(), b)
        assert_array_equal(a.transpose().todense(), self.dat)
        assert_array_equal(a.transpose().todense(), self.datsp.todense())

    def check_large(self):
        # Create a 100x100 matrix with 100 non-zero elements
        # and play around with it
        A = dok_matrix((100,100))
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

    def check_solve(self):
        """ Test whether the lu_solve command segfaults, as reported by Nils
        Wagner for a 64-bit machine, 02 March 2005 (EJS)
        """
        n = 20
        A = self.spmatrix((n,n), dtype=complex)
        x = numpy.rand(n)
        y = numpy.rand(n-1)+1j*numpy.rand(n-1)
        r = numpy.rand(n)
        for i in range(len(x)):
            A[i,i] = x[i]
        for i in range(len(y)):
            A[i,i+1] = y[i]
            A[i+1,i] = numpy.conjugate(y[i])
        B = A.tocsc()
        xx = splu(B).solve(r)
        # Don't actually test the output until we know what it should be ...


class _test_horiz_slicing(ScipyTestCase):
    """Tests vertical slicing (e.g. [:, 0]).  Tests for individual sparse
    matrix types that implement this should derive from this class.
    """
    def check_get_horiz_slice(self):
        """Test for new slice functionality (EJS)"""
        B = asmatrix(arange(50.).reshape(5,10))
        A = self.spmatrix(B)
        assert_array_equal(B[1,:], A[1,:].todense())
        assert_array_equal(B[1,2:5], A[1,2:5].todense())

        C = matrix([[1, 2, 1], [4, 0, 6], [0, 0, 0], [0, 0, 1]])
        D = self.spmatrix(C)
        assert_array_equal(C[1, 1:3], D[1, 1:3].todense())

        # Now test slicing when a row contains only zeros
        E = matrix([[1, 2, 1], [4, 0, 0], [0, 0, 0], [0, 0, 1]])
        F = self.spmatrix(E)
        assert_array_equal(E[1, 1:3], F[1, 1:3].todense())
        assert_array_equal(E[2, -2:], F[2, -2:].A)
        
        # The following should raise exceptions:
        caught = 0
        try:
            a = A[:,11]
        except IndexError:
            caught += 1
        try:
            a = A[6,3:7]
        except IndexError:
            caught += 1
        assert caught == 2


class _test_vert_slicing(ScipyTestCase):
    """Tests vertical slicing (e.g. [:, 0]).  Tests for individual sparse
    matrix types that implement this should derive from this class.
    """
    def check_get_vert_slice(self):
        """Test for new slice functionality (EJS)"""
        B = asmatrix(arange(50.).reshape(5,10))
        A = self.spmatrix(B)
        assert_array_equal(B[2:5,0], A[2:5,0].todense())
        assert_array_equal(B[:,1], A[:,1].todense())

        C = matrix([[1, 2, 1], [4, 0, 6], [0, 0, 0], [0, 0, 1]])
        D = self.spmatrix(C)
        assert_array_equal(C[1:3, 1], D[1:3, 1].todense())
        assert_array_equal(C[:, 2], D[:, 2].todense())

        # Now test slicing when a column contains only zeros
        E = matrix([[1, 0, 1], [4, 0, 0], [0, 0, 0], [0, 0, 1]])
        F = self.spmatrix(E)
        assert_array_equal(E[:, 1], F[:, 1].todense())
        assert_array_equal(E[-2:, 2], F[-2:, 2].todense())
        
        # The following should raise exceptions:
        caught = 0
        try:
            a = A[:,11]
        except IndexError:
            caught += 1
        try:
            a = A[6,3:7]
        except IndexError:
            caught += 1
        assert caught == 2


class _test_fancy_indexing(ScipyTestCase):
    """Tests fancy indexing features.  The tests for any matrix formats
    that implement these features should derive from this class.
    """
    # This isn't supported by any matrix objects yet:
    def check_sequence_indexing(self):
        B = asmatrix(arange(50.).reshape(5,10))
        A = self.spmatrix(B)
        assert_array_equal(B[(1,2),(3,4)], A[(1,2),(3,4)].todense())
        assert_array_equal(B[(1,2,3),(3,4,5)], A[(1,2,3),(3,4,5)].todense())

    def check_fancy_indexing(self):
        """Test for new indexing functionality"""
        B = ones((5,10), float)
        A = dok_matrix(B)
        # Write me!
        
        # Both slicing and fancy indexing: not yet supported
        # assert_array_equal(B[(1,2),:], A[(1,2),:].todense())
        # assert_array_equal(B[(1,2,3),:], A[(1,2,3),:].todense())



class test_csr(_test_cs, _test_horiz_slicing):
    spmatrix = csr_matrix

    def check_constructor1(self):
        b = matrix([[0,4,0],
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
        b = matrix([[1,0],
                   [0,2],
                   [3,0]],'d')
        bsp = csr_matrix(b)
        assert_array_almost_equal(bsp.data,[1,2,3])
        assert_array_equal(bsp.colind,[0,1,0])
        assert_array_equal(bsp.indptr,[0,1,2,3])
        assert_array_almost_equal(bsp.todense(),b)
    
    def check_empty(self):
        """Test manipulating empty matrices. Fails in SciPy SVN <= r1768
        """
        # This test should be made global (in _test_cs), but first we
        # need a uniform argument order / syntax for constructing an
        # empty sparse matrix. (coo_matrix is currently different).
        shape = (5, 5)
        for mytype in [float32, float64, complex64, complex128]:
            a = self.spmatrix(shape, dtype=mytype)
            b = a + a
            c = 2 * a
            d = a + a.tocsc()
            e = a * a.tocoo()
            assert_equal(e.A, a.A*a.A)
            # These fail in all revisions <= r1768:
            assert(e.dtype.type == mytype)
            assert(e.A.dtype.type == mytype)


class test_csc(_test_cs, _test_vert_slicing):
    spmatrix = csc_matrix

    def check_constructor1(self):
        b = matrix([[1,0,0],[3,0,1],[0,2,0]],'d')
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
        b = matrix([[1,0],[0,2],[3,0]],'d')
        bsp = csc_matrix(b)
        assert_array_almost_equal(bsp.data,[1,3,2])
        assert_array_equal(bsp.rowind,[0,2,1])
        assert_array_equal(bsp.indptr,[0,2,3])

    def check_empty(self):
        """Test manipulating empty matrices. Fails in SciPy SVN <= r1768
        """
        # This test should be made global (in _test_cs), but first we
        # need a uniform argument order / syntax for constructing an
        # empty sparse matrix. (coo_matrix is currently different).
        shape = (5, 5)
        for mytype in [float32, float64, complex64, complex128]:
            a = self.spmatrix(shape, dtype=mytype)
            b = a + a
            c = 2 * a
            d = a + a.tocsc()
            e = a * a.tocoo()
            assert_equal(e.A, a.A*a.A)
            assert(e.dtype.type == mytype)
            assert(e.A.dtype.type == mytype)


class test_dok(_test_cs):
    spmatrix = dok_matrix

    def check_mult(self):
        A = dok_matrix((10,10))
        A[0,3] = 10
        A[5,6] = 20
        D = A*A.T
        E = A*A.H
        assert_array_equal(D.A, E.A)

    def check_add(self):
        A = dok_matrix((3,2))
        A[0,1] = -10
        A[2,0] = 20
        A += 10
        B = matrix([[10, 0], [10, 10], [30, 10]])
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
        assert_array_equal( a.toarray()[:,n-1], zeros(m,) )

        # make sure it still works for CSC format
        csc=a.tocsc()
        assert_array_equal( csc.toarray()[:,n-1], zeros(m,) )

        # now test CSR
        (m, n) = (n, m)
        b = a.transpose()
        assert b.shape == (m, n)
        # assert that the last row is all zeros
        assert_array_equal( b.toarray()[m-1,:], zeros(n,) )

        # make sure it still works for CSR format
        csr=b.tocsr()
        assert_array_equal( csr.toarray()[m-1,:], zeros(n,))

    def check_set_slice(self):
        """Test for slice functionality (EJS)"""
        A = dok_matrix((5,10))
        B = zeros((5,10), float)
        A[:,0] = 1
        B[:,0] = 1
        assert_array_equal(A.todense(), B)
        A[1,:] = 2
        B[1,:] = 2
        assert_array_equal(A.todense(), B)
        A[:,:] = 3
        B[:,:] = 3
        assert_array_equal(A.todense(), B)
        A[1:5, 3] = 4
        B[1:5, 3] = 4
        assert_array_equal(A.todense(), B)
        A[1, 3:6] = 5
        B[1, 3:6] = 5
        assert_array_equal(A.todense(), B)
        A[1:4, 3:6] = 6
        B[1:4, 3:6] = 6
        assert_array_equal(A.todense(), B)
        A[1, 3:10:3] = 7
        B[1, 3:10:3] = 7
        assert_array_equal(A.todense(), B)
        A[1:5, 0] = range(1,5)
        B[1:5, 0] = range(1,5)
        assert_array_equal(A.todense(), B)
        A[0, 1:10:2] = xrange(1,10,2)
        B[0, 1:10:2] = xrange(1,10,2)
        assert_array_equal(A.todense(), B)
        caught = 0
        # The next 6 commands should raise exceptions
        try:
            A[0,0] = range(100)
        except TypeError:
            caught += 1
        try:
            A[0,0] = arange(100)
        except TypeError:
            caught += 1
        try:
            A[0,:] = range(100)
        except ValueError:
            caught += 1
        try:
            A[:,1] = range(100)
        except ValueError:
            caught += 1
        try:
            A[:,1] = A.copy()
        except:
            caught += 1
        try:
            A[:,-1] = range(5)
        except IndexError:
            caught += 1
        assert caught == 6


class test_lil(_test_cs, _test_horiz_slicing):
    spmatrix = lil_matrix
    def check_mult(self):
        A = matrix(zeros((10,10)))
        A[0,3] = 10
        A[5,6] = 20

        B = lil_matrix((10,10))
        B[0,3] = 10
        B[5,6] = 20
        assert_array_equal(A * A.T, (B * B.T).todense())
        assert_array_equal(A * A.H, (B * B.H).todense())

class test_construct_utils(ScipyTestCase):
    def check_identity(self):
        a = spidentity(3)
        b = array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype='d')
        assert_array_equal(a.toarray(), b)

    def check_eye(self):
        a = speye(2, 3 )
#        print a, a.__repr__
        b = array([[1, 0, 0], [0, 1, 0]], dtype='d')
        assert_array_equal(a.toarray(), b)

        a = speye(3, 2)
#        print a, a.__repr__
        b = array([[1, 0], [0, 1], [0, 0]], dtype='d')
        assert_array_equal( a.toarray(), b)

        a = speye(3, 3)
#        print a, a.__repr__
        b = array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype='d')
        assert_array_equal(a.toarray(), b)


if __name__ == "__main__":
    ScipyTest().run()
