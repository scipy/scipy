#!/usr/bin/env python
#
# Authors: Travis Oliphant, Ed Schofield, Robert Cimrman, Nathan Bell, and others

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
from numpy import arange, zeros, array, dot, ones, matrix, asmatrix, \
        asarray, vstack

import random
from numpy.testing import *
set_package_path()
from scipy.sparse import csc_matrix, csr_matrix, dok_matrix, coo_matrix, \
     spidentity, speye, spkron, extract_diagonal, lil_matrix, lil_eye, lil_diags
from scipy.linsolve import splu
restore_path()

class _test_cs:
    def setUp(self):
        self.dat = matrix([[1,0,0,2],[3,0,1,0],[0,2,0,0]],'d')
        self.datsp = self.spmatrix(self.dat)

    def check_getelement(self):
        assert_equal(self.datsp[0,0],1)
        assert_equal(self.datsp[0,1],0)
        assert_equal(self.datsp[1,0],3)
        assert_equal(self.datsp[2,1],2)

    def check_abs(self):
        A = matrix([[-1, 0, 17],[0, -5, 0],[1, -4, 0],[0,0,0]],'d')
        assert_equal(abs(A),abs(self.spmatrix(A)).todense())

    def check_neg(self):
        A = matrix([[-1, 0, 17],[0, -5, 0],[1, -4, 0],[0,0,0]],'d')
        assert_equal(-A,(-self.spmatrix(A)).todense())

    def check_sum(self):
        """Does the matrix's sum(,axis=0) method work?
        """
        assert_array_equal(self.dat.sum(), self.datsp.sum())
        assert_array_equal(self.dat.sum(axis=None), self.datsp.sum(axis=None))
        assert_array_equal(self.dat.sum(axis=0), self.datsp.sum(axis=0))
        assert_array_equal(self.dat.sum(axis=1), self.datsp.sum(axis=1))

    def check_mean(self):
        """Does the matrix's mean(,axis=0) method work?
        """
        assert_array_equal(self.dat.mean(), self.datsp.mean())
        assert_array_equal(self.dat.mean(axis=None), self.datsp.mean(axis=None))
        assert_array_equal(self.dat.mean(axis=0), self.datsp.mean(axis=0))
        assert_array_equal(self.dat.mean(axis=1), self.datsp.mean(axis=1))

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
        a[0,-1] = 8
        a[-1,-2] = 7
        assert_array_equal(a.todense(),[[0,3,0,8],[0,0,4,0],[2,0,7,0]])

    def check_mul_scalar(self):
        assert_array_equal(self.dat*2,(self.datsp*2).todense())
        assert_array_equal(self.dat*17.3,(self.datsp*17.3).todense())

    def check_imul_scalar(self):
        a = self.datsp.copy()
        a *= 2
        assert_array_equal(self.dat*2,a.todense())

        a = self.datsp.copy()
        a *= 17.3
        assert_array_equal(self.dat*17.3,a.todense())

    def check_idiv_scalar(self):
        a = self.datsp.copy()
        a /= 2
        assert_array_equal(self.dat/2,a.todense())

        a = self.datsp.copy()
        a /= 17.3
        assert_array_equal(self.dat/17.3,a.todense())

    def check_rmul_scalar(self):
        assert_array_equal(2*self.dat,(2*self.datsp).todense())
        assert_array_equal(17.3*self.dat,(17.3*self.datsp).todense())
        
    def check_add(self):
        a = self.datsp
        b = self.datsp.copy()
        b[0,2] = 2.0
        c = a + b
        assert_array_equal(c.todense(),[[2,0,2,4],[6,0,2,0],[0,4,0,0]])
    
    def check_radd(self):
        a = self.datsp
        b = self.datsp.copy()
        b[0,2] = 2.0
        c = a.todense() + b
        assert_array_equal(c,[[2,0,2,4],[6,0,2,0],[0,4,0,0]])
    
    def check_sub(self):
        assert_array_equal((self.datsp - self.datsp).todense(),[[0,0,0,0],[0,0,0,0],[0,0,0,0]])

        A = self.spmatrix(matrix([[1,0,0,4],[-1,0,0,0],[0,8,0,-5]],'d'))
        assert_array_equal((self.datsp - A).todense(),self.dat - A.todense())
        assert_array_equal((A - self.datsp).todense(),A.todense() - self.dat)
    
    def check_rsub(self):
        assert_array_equal((self.dat - self.datsp),[[0,0,0,0],[0,0,0,0],[0,0,0,0]])
        assert_array_equal((self.datsp - self.dat),[[0,0,0,0],[0,0,0,0],[0,0,0,0]])

        A = self.spmatrix(matrix([[1,0,0,4],[-1,0,0,0],[0,8,0,-5]],'d'))
        assert_array_equal((self.dat - A),self.dat - A.todense())
        assert_array_equal((A - self.dat),A.todense() - self.dat)
        assert_array_equal(A.todense() - self.datsp,A.todense() - self.dat)
        assert_array_equal(self.datsp - A.todense(),self.dat - A.todense())

    def check_elmul(self):
        a = self.datsp
        b = self.datsp.copy()
        b[0,2] = 2.0
        c = a ** b
        assert_array_equal(c.todense(),[[1,0,0,4],[9,0,1,0],[0,4,0,0]])
    
    def check_eldiv(self):
        assert_array_equal((self.datsp / self.datsp).todense(),[[1,0,0,1],[1,0,1,0],[0,1,0,0]])

        denom = self.spmatrix(matrix([[1,0,0,4],[-1,0,0,0],[0,8,0,-5]],'d'))
        res = matrix([[1,0,0,0.5],[-3,0,numpy.inf,0],[0,0.25,0,0]],'d')
        assert_array_equal((self.datsp / denom).todense(),res)

    def check_rmatvec(self):
        M = self.spmatrix(matrix([[3,0,0],[0,1,0],[2,0,3.0],[2,3,0]]))
        assert_array_almost_equal([1,2,3,4]*M, dot([1,2,3,4], M.toarray()))
        row = matrix([[1,2,3,4]])
        assert_array_almost_equal(row*M, row*M.todense())

    def check_matvec(self):
        M = self.spmatrix(matrix([[3,0,0],[0,1,0],[2,0,3.0],[2,3,0]]))
        col = matrix([1,2,3]).T
        assert_array_almost_equal(M * col, M.todense() * col)
        
        #check result dimensions (ticket #514) 
        assert_equal((M * array([1,2,3])).shape,(4,)) 
        assert_equal((M * array([[1],[2],[3]])).shape,(4,1)) 
        assert_equal((M * matrix([[1],[2],[3]])).shape,(4,1)) 

        #ensure exception is raised for improper dimensions
        bad_vecs = [array([1,2]), array([1,2,3,4]), array([[1],[2]]), 
                    matrix([1,2,3]), matrix([[1],[2]])]
        caught = 0
        for x in bad_vecs:
            try:
                y = M * x
            except ValueError:
                caught += 1
        assert_equal(caught,len(bad_vecs))

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

    def check_matmat_sparse(self):
        a = matrix([[3,0,0],[0,1,0],[2,0,3.0],[2,3,0]])
        a2 = array([[3,0,0],[0,1,0],[2,0,3.0],[2,3,0]])
        b = matrix([[0,1],[1,0],[0,2]],'d')
        asp = self.spmatrix(a)
        bsp = self.spmatrix(b)
        assert_array_almost_equal((asp*bsp).todense(), a*b)
        assert_array_almost_equal( asp*b, a*b)
        assert_array_almost_equal( a*bsp, a*b)
        assert_array_almost_equal( a2*bsp, a*b)

        # Now try performing cross-type multplication:
        csp = bsp.tocsc()
        c = b
        assert_array_almost_equal((asp*csp).todense(), a*c)
        assert_array_almost_equal((asp.matmat(csp)).todense(), a*c)
        assert_array_almost_equal( asp*c, a*c)
        
        assert_array_almost_equal( a*csp, a*c)
        assert_array_almost_equal( a2*csp, a*c)
        csp = bsp.tocsr()
        assert_array_almost_equal((asp*csp).todense(), a*c)
        assert_array_almost_equal((asp.matmat(csp)).todense(), a*c)
        assert_array_almost_equal( asp*c, a*c)

        assert_array_almost_equal( a*csp, a*c)
        assert_array_almost_equal( a2*csp, a*c)
        csp = bsp.tocoo()
        assert_array_almost_equal((asp*csp).todense(), a*c)
        assert_array_almost_equal((asp.matmat(csp)).todense(), a*c)
        assert_array_almost_equal( asp*c, a*c)

        assert_array_almost_equal( a*csp, a*c)
        assert_array_almost_equal( a2*csp, a*c)

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
    
    def check_matmat_dense(self):
        a = matrix([[3,0,0],[0,1,0],[2,0,3.0],[2,3,0]])
        asp = self.spmatrix(a)
        
        # check both array and matrix types
        bs = [ array([[1,2],[3,4],[5,6]]), matrix([[1,2],[3,4],[5,6]]) ]

        for b in bs:
            result = asp*b
            assert( isinstance(result, type(b)) )
            assert_equal( result.shape, (4,2) )
            assert_equal( result, dot(a,b) )
    
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

    def check_sub_dense(self):
        """ Check whether adding a dense matrix to a sparse matrix works
        """
        sum1 = 3*self.dat - self.datsp
        assert_array_equal(sum1, 2*self.dat)
        sum2 = 3*self.datsp - self.dat
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
        x = numpy.random.rand(n)
        y = numpy.random.rand(n-1)+1j*numpy.random.rand(n-1)
        r = numpy.random.rand(n)
        for i in range(len(x)):
            A[i,i] = x[i]
        for i in range(len(y)):
            A[i,i+1] = y[i]
            A[i+1,i] = numpy.conjugate(y[i])
        B = A.tocsc()
        xx = splu(B).solve(r)
        # Don't actually test the output until we know what it should be ...

    def check_extract_diagonal(self):
        """
        Test extraction of main diagonal from sparse matrices
        """
        L = []
        L.append(array([[0,0,3],[1,6,4],[5,2,0]]))
        L.append(array([[1,2,3]]))
        L.append(array([[7],[6],[5]]))
        L.append(array([[2]]))

        for A in L:
            assert_array_equal(numpy.diag(A),extract_diagonal(self.spmatrix(A)))
            

class _test_horiz_slicing:
    """Tests horizontal slicing (e.g. [:, 0]).  Tests for individual sparse
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


class _test_vert_slicing:
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

class _test_slicing:
    """Tests vertical and horizontal slicing (e.g. [:,0:2]). Tests for
    individual sparse matrix types that implement this should derive from this
    class.
    """
    def check_get_slices(self):
        B = asmatrix(arange(50.).reshape(5,10))
        A = self.spmatrix(B)
        assert_array_equal(B[2:5,0:3], A[2:5,0:3].todense())
        assert_array_equal(B[1:,:-1], A[1:,:-1].todense())
        assert_array_equal(B[:-1,1:], A[:-1,1:].todense())
        
        # Now test slicing when a column contains only zeros
        E = matrix([[1, 0, 1], [4, 0, 0], [0, 0, 0], [0, 0, 1]])
        F = self.spmatrix(E)
        assert_array_equal(E[1:2, 1:2], F[1:2, 1:2].todense())
        assert_array_equal(E[:, 1:], F[:, 1:].todense())

class _test_fancy_indexing:
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


class _test_arith:
    """
    Test real/complex arithmetic
    """
    def arith_init(self):
        #these can be represented exactly in FP (so arithmetic should be exact)
        self.A = matrix([[   -1.5,      0,       0,    2.25],
                         [  3.125,      0,  -0.125,       0],
                         [      0, -5.375,       0,       0]],'float64')
        self.B = matrix([[      0,  3.375,       0,  -7.875],
                         [  6.625,   4.75,       0,       0],
                         [    3.5, 6.0625,       0,       1]],'float64')
        
        self.C = matrix([[  0.375,       0,      -5,     2.5],
                         [      0,    7.25,       0,  -4.875],
                         [      0, -0.0625,       0,       0]],'complex128')
        self.C.imag = matrix([[    1.25,     0,  0, -3.875],
                              [       0, 4.125,  0,   2.75],
                              [ -0.0625,     0,  0,      1]],'float64')

        #fractions are all x/16ths
        assert_array_equal((self.A*16).astype('int32'),16*self.A)
        assert_array_equal((self.B*16).astype('int32'),16*self.B)
        assert_array_equal((self.C.real*16).astype('int32'),16*self.C.real)
        assert_array_equal((self.C.imag*16).astype('int32'),16*self.C.imag)

        self.Asp = self.spmatrix(self.A)
        self.Bsp = self.spmatrix(self.B)
        self.Csp = self.spmatrix(self.C)
    
        #supported types
        self.dtypes =  ['int8','uint8','int16','int32','int64',
                        'float32','float64','complex64','complex128']

    def check_conversion(self):
        self.arith_init()
        
        #check whether dtype and value is preserved in conversion
        for x in self.dtypes:
            A = self.A.astype(x)
            B = self.B.astype(x)
            C = self.C.astype(x)
                    
            Asp = self.spmatrix(A)
            Bsp = self.spmatrix(B)
            Csp = self.spmatrix(C)
            assert_equal(A.dtype,Asp.dtype)
            assert_equal(B.dtype,Bsp.dtype)
            assert_equal(C.dtype,Csp.dtype)
            assert_array_equal(A,Asp.todense())
            assert_array_equal(B,Bsp.todense())
            assert_array_equal(C,Csp.todense())
    
    def check_add_sub(self):
        self.arith_init()
        
        #basic tests
        assert_array_equal(self.A+self.B,(self.Asp+self.Bsp).todense())
        assert_array_equal(self.A+self.C,(self.Asp+self.Csp).todense())

        #check conversions
        for x in self.dtypes:
            for y in self.dtypes:
                A = self.A.astype(x)
                B = self.B.astype(y)
                C = self.C.astype(y)
                
                Asp = self.spmatrix(A)
                Bsp = self.spmatrix(B)
                Csp = self.spmatrix(C)

                #addition
                D1 = A + B
                D2 = A + C
                D3 = B + C
                S1 = Asp + Bsp
                S2 = Asp + Csp
                S3 = Bsp + Csp
                
                assert_equal(D1.dtype,S1.dtype)
                assert_equal(D2.dtype,S2.dtype)
                assert_equal(D3.dtype,S3.dtype)
                assert_array_equal(D1,S1.todense())
                assert_array_equal(D2,S2.todense())
                assert_array_equal(D3,S3.todense())
                assert_array_equal(D1,Asp + B)          #check sparse + dense
                assert_array_equal(D2,Asp + C)
                assert_array_equal(D3,Bsp + C)
                assert_array_equal(D1,A + Bsp)          #check dense + sparse
                assert_array_equal(D2,A + Csp)
                assert_array_equal(D3,B + Csp)

                #subtraction
                D1 = A - B
                D2 = A - C
                D3 = B - C
                S1 = Asp - Bsp
                S2 = Asp - Csp
                S3 = Bsp - Csp
                
                assert_equal(D1.dtype,S1.dtype)
                assert_equal(D2.dtype,S2.dtype)
                assert_equal(D3.dtype,S3.dtype)                    
                assert_array_equal(D1,S1.todense())
                assert_array_equal(D2,S2.todense())
                assert_array_equal(D3,S3.todense())
                assert_array_equal(D1,Asp - B)          #check sparse - dense
                assert_array_equal(D2,Asp - C)
                assert_array_equal(D3,Bsp - C)
                assert_array_equal(D1,A - Bsp)          #check dense - sparse
                try:
                    assert_array_equal(D2,A - Csp)
                except:
                    import pdb; pdb.set_trace()
                assert_array_equal(D3,B - Csp)


    def check_mu(self):
        self.arith_init()
        
        #basic tests
        assert_array_equal(self.A*self.B.T,(self.Asp*self.Bsp.T).todense())
        assert_array_equal(self.A*self.C.T,(self.Asp*self.Csp.T).todense())

        for x in self.dtypes:
            for y in self.dtypes:
                A = self.A.astype(x)
                B = self.B.astype(y)
                C = self.C.astype(y)
                
                Asp = self.spmatrix(A)
                Bsp = self.spmatrix(B)
                Csp = self.spmatrix(C)

                D1 = A * B.T
                D2 = A * C.T
                D3 = B * C.T

                S1 = Asp * Bsp.T
                S2 = Asp * Csp.T
                S3 = Bsp * Csp.T

                assert_array_equal(D1,S1.todense())
                assert_array_equal(D2,S2.todense())
                assert_array_equal(D3,S3.todense())
                assert_equal(D1.dtype,S1.dtype)
                assert_equal(D2.dtype,S2.dtype)
                assert_equal(D3.dtype,S3.dtype)




class test_csr(_test_cs, _test_horiz_slicing, _test_vert_slicing,
               _test_slicing, _test_arith, NumpyTestCase):
    spmatrix = csr_matrix

    def check_constructor1(self):
        b = matrix([[0,4,0],
                   [3,0,0],
                   [0,2,0]],'d')
        bsp = csr_matrix(b)
        assert_array_almost_equal(bsp.data,[4,3,2])
        assert_array_equal(bsp.indices,[1,0,1])
        assert_array_equal(bsp.indptr,[0,1,2,3])
        assert_equal(bsp.getnnz(),3)
        assert_equal(bsp.getformat(),'csr')
        assert_array_equal(bsp.todense(),b)
    
    def check_constructor2(self):
        b = zeros((6,6),'d')
        b[3,4] = 5
        bsp = csr_matrix(b)
        assert_array_almost_equal(bsp.data,[5])
        assert_array_equal(bsp.indices,[4])
        assert_array_equal(bsp.indptr,[0,0,0,0,1,1,1])
        assert_array_almost_equal(bsp.todense(),b)
    
    def check_constructor3(self):
        b = matrix([[1,0],
                   [0,2],
                   [3,0]],'d')
        bsp = csr_matrix(b)
        assert_array_almost_equal(bsp.data,[1,2,3])
        assert_array_equal(bsp.indices,[0,1,0])
        assert_array_equal(bsp.indptr,[0,1,2,3])
        assert_array_almost_equal(bsp.todense(),b)
    
    def check_constructor4(self):
        """try using int64 indices"""
        data = arange( 6 ) + 1
        col = array( [1, 2, 1, 0, 0, 2], dtype='int64' )
        ptr = array( [0, 2, 4, 6], dtype='int64' )

        a = csr_matrix( (data, col, ptr), dims = (3,3) )
        
        b = matrix([[0,1,2],
                    [4,3,0],
                    [5,0,6]],'d')
        
        assert_equal(a.indptr.dtype,numpy.dtype('int64'))
        assert_equal(a.indices.dtype,numpy.dtype('int64'))
        assert_array_equal(a.todense(),b)

    def check_constructor5(self):
        """using (data, ij) format"""
        row  = numpy.array([2, 3, 1, 3, 0, 1, 3, 0, 2, 1, 2])
        col  = numpy.array([0, 1, 0, 0, 1, 1, 2, 2, 2, 2, 1])
        data = numpy.array([  6.,  10.,   3.,   9.,   1.,   4.,
                              11.,   2.,   8.,   5.,   7.])
        
        ij = vstack((row,col))
        csr = csr_matrix((data,ij),(4,3))                
        assert_array_equal(arange(12).reshape(4,3),csr.todense())


    def check_empty(self):
        """Test manipulating empty matrices. Fails in SciPy SVN <= r1768
        """
        # This test should be made global (in _test_cs), but first we
        # need a uniform argument order / syntax for constructing an
        # empty sparse matrix. (coo_matrix is currently different).
        shape = (5, 5)
        for mytype in ['int32', 'float32', 'float64', 'complex64', 'complex128']:
            a = self.spmatrix(shape, dtype=mytype)
            b = a + a
            c = 2 * a
            d = a + a.tocsc()
            e = a * a.tocoo()
            assert_equal(e.A, a.A*a.A)
            # These fail in all revisions <= r1768:
            assert_equal(e.dtype,mytype)
            assert_equal(e.A.dtype,mytype)

    def check_ensure_sorted_indices(self):
        #print 'sorting CSR indices'
        data = arange( 5 )
        col = array( [7, 2, 1, 5, 4] )
        ptr = [0, 3, 5]
        asp = csr_matrix( (data, col, ptr), dims = (2,10) )
        bsp = asp.copy()
        #print 'in\n', asp
        asp.ensure_sorted_indices( inplace = True )
        #print 'out\n', asp
        assert_array_equal(asp.indices,[1, 2, 7, 4, 5])
        for ir in range( asp.shape[0] ):
            for ic in range( asp.shape[1] ):
                assert_equal( asp[ir, ic], bsp[ir, ic] )

    def check_get_submatrix(self):
        a = csr_matrix( array([[1,2,3,4],[1,2,3,5],[0,2,0,1]]) )
        i0 = slice( 0, 2 )
        i1 = ( 1, 3 )
        b = a.get_submatrix( i0, i1 )

        aa = a.toarray()
        ab = b.toarray()

        assert b.dtype == a.dtype
        assert b.shape == (2,2)
        assert_equal( ab, aa[i0,i1[0]:i1[1]] )

class test_csc(_test_cs, _test_horiz_slicing, _test_vert_slicing,
               _test_slicing, _test_arith, NumpyTestCase):
    spmatrix = csc_matrix

    def check_constructor1(self):
        b = matrix([[1,0,0],[3,0,1],[0,2,0]],'d')
        bsp = csc_matrix(b)
        assert_array_almost_equal(bsp.data,[1,3,2,1])
        assert_array_equal(bsp.indices,[0,1,2,1])
        assert_array_equal(bsp.indptr,[0,2,3,4])
        assert_equal(bsp.getnnz(),4)
        assert_equal(bsp.getformat(),'csc')

    def check_constructor2(self):
        b = zeros((6,6),'d')
        b[2,4] = 5
        bsp = csc_matrix(b)
        assert_array_almost_equal(bsp.data,[5])
        assert_array_equal(bsp.indices,[2])
        assert_array_equal(bsp.indptr,[0,0,0,0,0,1,1])

    def check_constructor3(self):
        b = matrix([[1,0],[0,0],[0,2]],'d')
        bsp = csc_matrix(b)
        assert_array_almost_equal(bsp.data,[1,2])
        assert_array_equal(bsp.indices,[0,2])
        assert_array_equal(bsp.indptr,[0,1,2])

    def check_constructor5(self):
        """using (data, ij) format"""
        row  = numpy.array([2, 3, 1, 3, 0, 1, 3, 0, 2, 1, 2])
        col  = numpy.array([0, 1, 0, 0, 1, 1, 2, 2, 2, 2, 1])
        data = numpy.array([  6.,  10.,   3.,   9.,   1.,   4.,
                              11.,   2.,   8.,   5.,   7.])
       
        ij = vstack((row,col))
        csc = csc_matrix((data,ij),(4,3))                
        assert_array_equal(arange(12).reshape(4,3),csc.todense())

    def check_empty(self):
        """Test manipulating empty matrices. Fails in SciPy SVN <= r1768
        """
        # This test should be made global (in _test_cs), but first we
        # need a uniform argument order / syntax for constructing an
        # empty sparse matrix. (coo_matrix is currently different).
        shape = (5, 5)
        for mytype in ['int32', 'float32', 'float64', 'complex64', 'complex128']:
            a = self.spmatrix(shape, dtype=mytype)
            b = a + a
            c = 2 * a
            d = a + a.tocsc()
            e = a * a.tocoo()
            assert_equal(e.A, a.A*a.A)
            assert_equal(e.dtype, mytype)
            assert_equal(e.A.dtype, mytype)

    def check_ensure_sorted_indices(self):
        #print 'sorting CSC indices'
        data = arange( 5 )
        row = array( [7, 2, 1, 5, 4] )
        ptr = [0, 3, 5]
        asp = csc_matrix( (data, row, ptr), dims = (10,2) )
        bsp = asp.copy()
        #print 'in\n', asp
        asp.ensure_sorted_indices( inplace = True )
        #print 'out\n', asp
        assert_array_equal(asp.indices,[1, 2, 7, 4, 5])
        for ir in range( asp.shape[0] ):
            for ic in range( asp.shape[1] ):
                assert_equal( asp[ir, ic], bsp[ir, ic] )

    def check_get_submatrix(self):
        a = csc_matrix( array([[1,2,3,4],[1,2,3,5],[0,2,0,1]]) )
        i0 = slice( 0, 2 )
        i1 = ( 1, 3 )
        b = a.get_submatrix( i0, i1 )

        aa = a.toarray()
        ab = b.toarray()

        assert_equal(b.dtype, a.dtype)
        assert_equal(b.shape, (2,2))
        assert_equal( ab, aa[i0,i1[0]:i1[1]] )

class test_dok(_test_cs, NumpyTestCase):
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
        A = A + 10
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
        assert_equal(b.shape, (m, n))
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
        assert_equal(caught,5)


class test_lil(_test_cs, _test_horiz_slicing, NumpyTestCase,
               ParametricTestCase):
    spmatrix = lil_matrix

    B = lil_matrix((4,3))
    B[0,0] = 2
    B[1,2] = 7
    B[2,1] = 3
    B[3,0] = 10

    def check_dot(self):
        A = matrix(zeros((10,10)))
        A[0,3] = 10
        A[5,6] = 20

        B = lil_matrix((10,10))
        B[0,3] = 10
        B[5,6] = 20
        assert_array_equal(A * A.T, (B * B.T).todense())
        assert_array_equal(A * A.H, (B * B.H).todense())

    def check_scalar_mul(self):
        x = lil_matrix((3,3))
        x[0,0] = 2

        x = x*2
        assert_equal(x[0,0],4)

        x = x*0
        assert_equal(x[0,0],0)

    def check_reshape(self):
        x = lil_matrix((4,3))
        x[0,0] = 1
        x[2,1] = 3
        x[3,2] = 5
        x[0,2] = 7

        for s in [(12,1),(1,12)]:
            assert_array_equal(x.reshape(s).todense(),
                               x.todense().reshape(s))

    def check_lil_lil_assignment(self):
        """ Tests whether a row of one lil_matrix can be assigned to
        another.
        """
        B = self.B.copy()
        A = B / 10
        B[0,:] = A[0,:]
        assert_array_equal(A[0,:].A, B[0,:].A)

    def tst_inplace_op(self,op,arr,other,result):
        cpy = arr
        getattr(arr,"__i%s__" % op)(other)

        assert_array_equal(cpy.todense(),arr.todense())
        assert_array_equal(arr.todense(),result)

    def testip_inplace_ops(self):
        B = self.B[:3,:3].copy()
        B[:,:] = B-B
        C = B.todense()

        data = {'add':(B,C+C),
                'sub':(B,zeros(B.shape)),
                'mul':(3,C*3)}

        return [(self.tst_inplace_op,op,B,other,result)
                for op,(other,result) in data.iteritems()]

    def check_lil_slice_assignment(self):
        B = lil_matrix((4,3))
        B[0,0] = 5
        B[1,2] = 3
        B[2,1] = 7

        expected = array([[10,0,0],
                          [0,0,6],
                          [0,14,0],
                          [0,0,0]])

        B[:,:] = B+B
        assert_array_equal(B.todense(),expected)

        block = [[1,0],[0,4]]
        B[:2,:2] = csc_matrix(array(block))
        assert_array_equal(B.todense()[:2,:2],block)

    def check_lil_sequence_assignement(self):
        A = lil_matrix((4,3))
        B = lil_eye((3,4))

        i0 = [0,1,2]
        i1 = (0,1,2)
        i2 = array( i0 )

        A[0,i0] = B[i0,0]
        A[1,i1] = B[i1,1]
        A[2,i2] = B[i2,2]
        assert_array_equal(A.todense(),B.T.todense())

    def check_lil_iteration(self):
        row_data = [[1,2,3],[4,5,6]]
        B = lil_matrix(array(row_data))
        for r,row in enumerate(B):
            assert_array_equal(row.todense(),array(row_data[r],ndmin=2))

    def check_lil_from_csr(self):
        """ Tests whether a lil_matrix can be constructed from a
        csr_matrix.
        """
        B = lil_matrix((10,10))
        B[0,3] = 10
        B[5,6] = 20
        B[8,3] = 30
        B[3,8] = 40
        B[8,9] = 50
        C = B.tocsr()
        D = lil_matrix(C)
        assert_array_equal(C.A, D.A)

    def check_point_wise_multiply(self):
        l = lil_matrix((4,3))
        l[0,0] = 1
        l[1,1] = 2
        l[2,2] = 3
        l[3,1] = 4

        m = lil_matrix((4,3))
        m[0,0] = 1
        m[0,1] = 2
        m[2,2] = 3
        m[3,1] = 4
        m[3,2] = 4

        assert_array_equal(l.multiply(m).todense(),
                           m.multiply(l).todense())

        assert_array_equal(l.multiply(m).todense(),
                           [[1,0,0],
                            [0,0,0],
                            [0,0,9],
                            [0,16,0]])

    def check_lil_eye(self):
        for dim in [(3,5),(5,3)]:
            for k in range(-5,5):
                r,c = dim
                assert_array_equal(lil_eye(dim,k).todense(),
                                   speye(r,c,k).todense())

    def check_lil_diags(self):
        assert_array_equal(lil_diags([[1,2,3],[4,5],[6]],
                                     [0,1,2],(3,3)).todense(),
                           [[1,4,6],
                            [0,2,5],
                            [0,0,3]])

        assert_array_equal(lil_diags([[6],[4,5],[1,2,3]],
                                     [2,1,0],(3,3)).todense(),
                           [[1,4,6],
                            [0,2,5],
                            [0,0,3]])

        assert_array_equal(lil_diags([[6,7,8],[4,5],[1,2,3]],
                                     [2,1,0],(3,3)).todense(),
                           [[1,4,6],
                            [0,2,5],
                            [0,0,3]])

        assert_array_equal(lil_diags([[1,2,3],[4,5],[6]],
                                     [0,-1,-2],(3,3)).todense(),
                           [[1,0,0],
                            [4,2,0],
                            [6,5,3]])

        assert_array_equal(lil_diags([[6,7,8],[4,5]],
                                     [-2,-1],(3,3)).todense(),
                           [[0,0,0],
                            [4,0,0],
                            [6,5,0]])

class TestConstructUtils(NumpyTestCase):
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

    def check_spkron(self):
        from numpy import kron

        cases = []

        cases.append(array([[ 0]]))
        cases.append(array([[-1]]))
        cases.append(array([[ 4]]))
        cases.append(array([[10]]))
        cases.append(array([[0],[0]]))
        cases.append(array([[0,0]]))
        cases.append(array([[1,2],[3,4]]))
        cases.append(array([[0,2],[5,0]]))
        cases.append(array([[0,2,-6],[8,0,14]]))
        cases.append(array([[5,4],[0,0],[6,0]]))
        cases.append(array([[5,4,4],[1,0,0],[6,0,8]]))
        cases.append(array([[0,1,0,2,0,5,8]]))

        for a in cases:
            for b in cases:
                result = spkron(csr_matrix(a),csr_matrix(b)).todense()
                expected = kron(a,b)
        
                assert_array_equal(result,expected)

class TestCOO(NumpyTestCase):
    def check_constructor1(self):
        """unsorted triplet format"""
        row  = numpy.array([2, 3, 1, 3, 0, 1, 3, 0, 2, 1, 2])
        col  = numpy.array([0, 1, 0, 0, 1, 1, 2, 2, 2, 2, 1])
        data = numpy.array([  6.,  10.,   3.,   9.,   1.,   4.,
                              11.,   2.,   8.,   5.,   7.])
        
        coo = coo_matrix((data,(row,col)),(4,3))
        
        assert_array_equal(arange(12).reshape(4,3),coo.todense())

    def check_constructor2(self):
        """unsorted triplet format with duplicates (which are summed)"""
        row  = numpy.array([0,1,2,2,2,2,0,0,2,2])
        col  = numpy.array([0,2,0,2,1,1,1,0,0,2])
        data = numpy.array([2,9,-4,5,7,0,-1,2,1,-5])
        coo = coo_matrix((data,(row,col)),(3,3))

        mat = matrix([[4,-1,0],[0,0,9],[-3,7,0]])
        
        assert_array_equal(mat,coo.todense())

    def check_constructor3(self):
        """empty matrix"""
        coo = coo_matrix(None,dims=(4,3))
        
        assert_array_equal(zeros((4,3)),coo.todense())

    def check_constructor4(self):
        """from dense matrix"""
        mat = numpy.array([[0,1,0,0],
                           [7,0,3,0],
                           [0,4,0,0]])
        coo = coo_matrix(mat)
        assert_array_equal(mat,coo.todense())
        
##    def check_normalize( self ):
##        row  = numpy.array([2, 3, 1, 3, 0, 1, 3, 0, 2, 1, 2])
##        col  = numpy.array([0, 1, 0, 0, 1, 1, 2, 2, 2, 2, 1])
##        data = numpy.array([  6.,  10.,   3.,   9.,   1.,   4.,
##                              11.,   2.,   8.,   5.,   7.])
##
##        # coo.todense()
##        #    matrix([[  0.,   1.,   2.],
##        #            [  3.,   4.,   5.],
##        #            [  6.,   7.,   8.],
##        #            [  9.,  10.,  11.]])
##        coo = coo_matrix((data,(row,col)),(4,3))
##
##        ndata,nrow,ncol = coo._normalize(rowfirst=True)
##        sorted_rcd = zip(row, col, data)
##        sorted_rcd.sort()
##        assert(zip(nrow,ncol,ndata) == sorted_rcd) #should sort by rows, then cols
##        assert_array_equal(coo.data, data)                        #coo.data has not changed
##        assert_array_equal(coo.row, row)                          #coo.row has not changed
##        assert_array_equal(coo.col, col)                          #coo.col has not changed
##
##
##        ndata,nrow,ncol = coo._normalize(rowfirst=False)
##        assert(zip(ncol,nrow,ndata) == sorted(zip(col,row,data))) #should sort by cols, then rows
##        assert_array_equal(coo.data, ndata)                       #coo.data has changed
##        assert_array_equal(coo.row, nrow)                         #coo.row has changed
##        assert_array_equal(coo.col, ncol)                         #coo.col has changed
##
##        assert_array_equal(coo.tocsr().todense(), coo.todense())
##        assert_array_equal(coo.tocsc().todense(), coo.todense())


if __name__ == "__main__":
    NumpyTest().run()
