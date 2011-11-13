#
# Authors: Travis Oliphant, Ed Schofield, Robert Cimrman, Nathan Bell, and others

""" Test functions for sparse matrices

"""
__usage__ = """
Build sparse:
  python setup.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.sparse.test()'
Run tests if sparse is not installed:
  python tests/test_sparse.py
"""

import sys
import warnings

import numpy as np
from numpy import arange, zeros, array, dot, matrix, asmatrix, asarray, \
                  vstack, ndarray, transpose, diag, kron, inf, conjugate, \
                  int8, ComplexWarning

import random
from numpy.testing import assert_raises, assert_equal, assert_array_equal, \
        assert_array_almost_equal, assert_almost_equal, assert_, \
        dec, TestCase, run_module_suite

import scipy.sparse as sparse
from scipy.sparse import csc_matrix, csr_matrix, dok_matrix, \
        coo_matrix, lil_matrix, dia_matrix, bsr_matrix, \
        eye, isspmatrix, SparseEfficiencyWarning
from scipy.sparse.sputils import supported_dtypes
from scipy.sparse.linalg import splu


warnings.simplefilter('ignore', SparseEfficiencyWarning)
warnings.simplefilter('ignore', ComplexWarning)


#TODO check that spmatrix( ... , copy=X ) is respected
#TODO test prune
#TODO test has_sorted_indices
class _TestCommon:
    """test common functionality shared by all sparse formats"""

    def setUp(self):
        self.dat = matrix([[1,0,0,2],[3,0,1,0],[0,2,0,0]],'d')
        self.datsp = self.spmatrix(self.dat)

    def test_empty(self):
        """create empty matrices"""

        assert_equal(self.spmatrix((3,3)).todense(), np.zeros((3,3)))
        assert_equal(self.spmatrix((3,3)).nnz, 0)

    def test_invalid_shapes(self):
        assert_raises(ValueError, self.spmatrix, (-1,3) )
        assert_raises(ValueError, self.spmatrix, (3,-1) )
        assert_raises(ValueError, self.spmatrix, (-1,-1) )

    def test_repr(self):
        repr(self.datsp)

    def test_str(self):
        str(self.datsp)

    def test_empty_arithmetic(self):
        """Test manipulating empty matrices. Fails in SciPy SVN <= r1768
        """
        shape = (5, 5)
        for mytype in [np.dtype('int32'), np.dtype('float32'),
                np.dtype('float64'), np.dtype('complex64'),
                np.dtype('complex128')]:
            a = self.spmatrix(shape, dtype=mytype)
            b = a + a
            c = 2 * a
            d = a * a.tocsc()
            e = a * a.tocsr()
            f = a * a.tocoo()
            for m in [a,b,c,d,e,f]:
                assert_equal(m.A, a.A*a.A)
                # These fail in all revisions <= r1768:
                assert_equal(m.dtype,mytype)
                assert_equal(m.A.dtype,mytype)

    def test_abs(self):
        A = matrix([[-1, 0, 17],[0, -5, 0],[1, -4, 0],[0,0,0]],'d')
        assert_equal(abs(A),abs(self.spmatrix(A)).todense())

    def test_neg(self):
        A = matrix([[-1, 0, 17],[0, -5, 0],[1, -4, 0],[0,0,0]],'d')
        assert_equal(-A,(-self.spmatrix(A)).todense())

    def test_real(self):
        D = matrix([[1 + 3j, 2 - 4j]])
        A = self.spmatrix(D)
        assert_equal(A.real.todense(),D.real)

    def test_imag(self):
        D = matrix([[1 + 3j, 2 - 4j]])
        A = self.spmatrix(D)
        assert_equal(A.imag.todense(),D.imag)

    def test_diagonal(self):
        """Does the matrix's .diagonal() method work?
        """
        mats = []
        mats.append( [[1,0,2]] )
        mats.append( [[1],[0],[2]] )
        mats.append( [[0,1],[0,2],[0,3]] )
        mats.append( [[0,0,1],[0,0,2],[0,3,0]] )

        mats.append( kron(mats[0],[[1,2]]) )
        mats.append( kron(mats[0],[[1],[2]]) )
        mats.append( kron(mats[1],[[1,2],[3,4]]) )
        mats.append( kron(mats[2],[[1,2],[3,4]]) )
        mats.append( kron(mats[3],[[1,2],[3,4]]) )
        mats.append( kron(mats[3],[[1,2,3,4]]) )

        for m in mats:
            assert_equal(self.spmatrix(m).diagonal(),diag(m))


    def test_nonzero(self):
        A   = array([[1, 0, 1],[0, 1, 1],[ 0, 0, 1]])
        Asp = self.spmatrix(A)

        A_nz   = set( [tuple(ij) for ij in transpose(A.nonzero())] )
        Asp_nz = set( [tuple(ij) for ij in transpose(Asp.nonzero())] )

        assert_equal(A_nz, Asp_nz)


    def test_getrow(self):
        assert_array_equal(self.datsp.getrow(1).todense(), self.dat[1,:])
        assert_array_equal(self.datsp.getrow(-1).todense(), self.dat[-1,:])

    def test_getcol(self):
        assert_array_equal(self.datsp.getcol(1).todense(), self.dat[:,1])
        assert_array_equal(self.datsp.getcol(-1).todense(), self.dat[:,-1])

    def test_sum(self):
        """Does the matrix's .sum(axis=...) method work?
        """
        assert_array_equal(self.dat.sum(), self.datsp.sum())
        assert_array_equal(self.dat.sum(axis=None), self.datsp.sum(axis=None))
        assert_array_equal(self.dat.sum(axis=0), self.datsp.sum(axis=0))
        assert_array_equal(self.dat.sum(axis=1), self.datsp.sum(axis=1))

    def test_mean(self):
        """Does the matrix's .mean(axis=...) method work?
        """
        assert_array_equal(self.dat.mean(), self.datsp.mean())
        assert_array_equal(self.dat.mean(axis=None), self.datsp.mean(axis=None))
        assert_array_equal(self.dat.mean(axis=0), self.datsp.mean(axis=0))
        assert_array_equal(self.dat.mean(axis=1), self.datsp.mean(axis=1))

    def test_from_array(self):
        A = array([[1,0,0],[2,3,4],[0,5,0],[0,0,0]])
        assert_array_equal(self.spmatrix(A).toarray(), A)

        A = array([[1.0 + 3j,       0,      0],
                   [       0, 2.0 + 5,      0],
                   [       0,       0,      0]])
        assert_array_equal(self.spmatrix(A).toarray(), A)
        assert_array_equal(self.spmatrix(A, dtype='int16').toarray(), A.astype('int16'))

    def test_from_matrix(self):
        A = matrix([[1,0,0],[2,3,4],[0,5,0],[0,0,0]])
        assert_array_equal(self.spmatrix(A).todense(), A)

        A = matrix([[1.0 + 3j,       0,      0],
                    [       0, 2.0 + 5,      0],
                    [       0,       0,      0]])
        assert_array_equal(self.spmatrix(A).toarray(), A)
        assert_array_equal(self.spmatrix(A, dtype='int16').toarray(), A.astype('int16'))

    def test_from_list(self):
        A = [[1,0,0],[2,3,4],[0,5,0],[0,0,0]]
        assert_array_equal(self.spmatrix(A).todense(), A)

        A = [[1.0 + 3j,       0,      0],
             [       0, 2.0 + 5,      0],
             [       0,       0,      0]]
        assert_array_equal(self.spmatrix(A).toarray(), array(A))
        assert_array_equal(self.spmatrix(A, dtype='int16').todense(), array(A).astype('int16'))

    def test_from_sparse(self):
        D = array([[1,0,0],[2,3,4],[0,5,0],[0,0,0]])
        S = csr_matrix(D)
        assert_array_equal(self.spmatrix(S).toarray(), D)
        S = self.spmatrix(D)
        assert_array_equal(self.spmatrix(S).toarray(), D)


        D = array([[1.0 + 3j,       0,      0],
                   [       0, 2.0 + 5,      0],
                   [       0,       0,      0]])
        S = csr_matrix(D)
        assert_array_equal(self.spmatrix(S).toarray(), D)
        assert_array_equal(self.spmatrix(S, dtype='int16').toarray(), D.astype('int16'))
        S = self.spmatrix(D)
        assert_array_equal(self.spmatrix(S).toarray(), D)
        assert_array_equal(self.spmatrix(S, dtype='int16').toarray(), D.astype('int16'))

    #def test_array(self):
    #    """test array(A) where A is in sparse format"""
    #    assert_equal( array(self.datsp), self.dat )

    def test_todense(self):
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

    def test_toarray(self):
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

    def test_astype(self):
        D = array([[1.0 + 3j,       0,      0],
                   [       0, 2.0 + 5,      0],
                   [       0,       0,      0]])
        S = self.spmatrix(D)

        for x in supported_dtypes:
            assert_equal(S.astype(x).dtype,     D.astype(x).dtype)  # correct type
            assert_equal(S.astype(x).toarray(), D.astype(x))        # correct values
            assert_equal(S.astype(x).format,    S.format)           # format preserved

    def test_asfptype(self):
        A = self.spmatrix( arange(6,dtype='int32').reshape(2,3) )

        assert_equal( A.dtype , np.dtype('int32') )
        assert_equal( A.asfptype().dtype, np.dtype('float64') )
        assert_equal( A.asfptype().format, A.format )
        assert_equal( A.astype('int16').asfptype().dtype , np.dtype('float32') )
        assert_equal( A.astype('complex128').asfptype().dtype , np.dtype('complex128') )

        B = A.asfptype()
        C = B.asfptype()
        assert_( B is C )


    def test_mul_scalar(self):
        assert_array_equal(self.dat*2,(self.datsp*2).todense())
        assert_array_equal(self.dat*17.3,(self.datsp*17.3).todense())

    def test_rmul_scalar(self):
        assert_array_equal(2*self.dat,(2*self.datsp).todense())
        assert_array_equal(17.3*self.dat,(17.3*self.datsp).todense())

    def test_add(self):
        a = self.dat.copy()
        a[0,2] = 2.0
        b = self.datsp
        c = b + a
        assert_array_equal(c,[[2,0,2,4],[6,0,2,0],[0,4,0,0]])

    def test_radd(self):
        a = self.dat.copy()
        a[0,2] = 2.0
        b = self.datsp
        c = a + b
        assert_array_equal(c,[[2,0,2,4],[6,0,2,0],[0,4,0,0]])

    def test_sub(self):
        assert_array_equal((self.datsp - self.datsp).todense(),[[0,0,0,0],[0,0,0,0],[0,0,0,0]])

        A = self.spmatrix(matrix([[1,0,0,4],[-1,0,0,0],[0,8,0,-5]],'d'))
        assert_array_equal((self.datsp - A).todense(),self.dat - A.todense())
        assert_array_equal((A - self.datsp).todense(),A.todense() - self.dat)

    def test_rsub(self):
        assert_array_equal((self.dat - self.datsp),[[0,0,0,0],[0,0,0,0],[0,0,0,0]])
        assert_array_equal((self.datsp - self.dat),[[0,0,0,0],[0,0,0,0],[0,0,0,0]])

        A = self.spmatrix(matrix([[1,0,0,4],[-1,0,0,0],[0,8,0,-5]],'d'))
        assert_array_equal((self.dat - A),self.dat - A.todense())
        assert_array_equal((A - self.dat),A.todense() - self.dat)
        assert_array_equal(A.todense() - self.datsp,A.todense() - self.dat)
        assert_array_equal(self.datsp - A.todense(),self.dat - A.todense())

    def test_elementwise_multiply(self):
        # real/real
        A = array([[4,0,9],[2,-3,5]])
        B = array([[0,7,0],[0,-4,0]])
        Asp = self.spmatrix(A)
        Bsp = self.spmatrix(B)
        assert_almost_equal( Asp.multiply(Bsp).todense(), A*B) #sparse/sparse
        assert_almost_equal( Asp.multiply(B),             A*B) #sparse/dense

        # complex/complex
        C = array([[1-2j,0+5j,-1+0j],[4-3j,-3+6j,5]])
        D = array([[5+2j,7-3j,-2+1j],[0-1j,-4+2j,9]])
        Csp = self.spmatrix(C)
        Dsp = self.spmatrix(D)
        assert_almost_equal( Csp.multiply(Dsp).todense(), C*D) #sparse/sparse
        assert_almost_equal( Csp.multiply(D),             C*D) #sparse/dense

        # real/complex
        assert_almost_equal( Asp.multiply(Dsp).todense(), A*D) #sparse/sparse
        assert_almost_equal( Asp.multiply(D),             A*D) #sparse/dense


    def test_elementwise_divide(self):
        expected = [[1,0,0,1],[1,0,1,0],[0,1,0,0]]
        assert_array_equal((self.datsp / self.datsp).todense(),expected)

        denom = self.spmatrix(matrix([[1,0,0,4],[-1,0,0,0],[0,8,0,-5]],'d'))
        res = matrix([[1,0,0,0.5],[-3,0,inf,0],[0,0.25,0,0]],'d')
        assert_array_equal((self.datsp / denom).todense(),res)

        # complex
        A = array([[1-2j,0+5j,-1+0j],[4-3j,-3+6j,5]])
        B = array([[5+2j,7-3j,-2+1j],[0-1j,-4+2j,9]])
        Asp = self.spmatrix(A)
        Bsp = self.spmatrix(B)
        assert_almost_equal( (Asp / Bsp).todense(), A/B)

    def test_pow(self):
        A = matrix([[1,0,2,0],[0,3,4,0],[0,5,0,0],[0,6,7,8]])
        B = self.spmatrix( A )

        for exponent in [0,1,2,3]:
            assert_array_equal((B**exponent).todense(),A**exponent)

        #invalid exponents
        for exponent in [-1, 2.2, 1 + 3j]:
            self.assertRaises( Exception, B.__pow__, exponent )

        #nonsquare matrix
        B = self.spmatrix(A[:3,:])
        self.assertRaises( Exception, B.__pow__, 1 )


    def test_rmatvec(self):
        M = self.spmatrix(matrix([[3,0,0],[0,1,0],[2,0,3.0],[2,3,0]]))
        assert_array_almost_equal([1,2,3,4]*M, dot([1,2,3,4], M.toarray()))
        row = matrix([[1,2,3,4]])
        assert_array_almost_equal(row*M, row*M.todense())

    def test_small_multiplication(self):
        """test that A*x works for x with shape () (1,) and (1,1)
        """
        A = self.spmatrix([[1],[2],[3]])

        assert_(isspmatrix(A * array(1)))
        assert_equal((A * array(1)).todense(), [[1],[2],[3]])
        assert_equal(A * array([1]), array([1,2,3]))
        assert_equal(A * array([[1]]), array([[1],[2],[3]]))

    def test_matvec(self):
        M = self.spmatrix(matrix([[3,0,0],[0,1,0],[2,0,3.0],[2,3,0]]))
        col = matrix([1,2,3]).T
        assert_array_almost_equal(M * col, M.todense() * col)

        #check result dimensions (ticket #514)
        assert_equal((M * array([1,2,3])).shape,(4,))
        assert_equal((M * array([[1],[2],[3]])).shape,(4,1))
        assert_equal((M * matrix([[1],[2],[3]])).shape,(4,1))

        #check result type
        assert_(isinstance( M * array([1,2,3]), ndarray))
        assert_(isinstance( M * matrix([1,2,3]).T, matrix))

        #ensure exception is raised for improper dimensions
        bad_vecs = [array([1,2]), array([1,2,3,4]), array([[1],[2]]),
                    matrix([1,2,3]), matrix([[1],[2]])]
        for x in bad_vecs:
            assert_raises(ValueError, M.__mul__, x)

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

    def test_matmat_sparse(self):
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
        assert_array_almost_equal( asp*c, a*c)

        assert_array_almost_equal( a*csp, a*c)
        assert_array_almost_equal( a2*csp, a*c)
        csp = bsp.tocsr()
        assert_array_almost_equal((asp*csp).todense(), a*c)
        assert_array_almost_equal( asp*c, a*c)

        assert_array_almost_equal( a*csp, a*c)
        assert_array_almost_equal( a2*csp, a*c)
        csp = bsp.tocoo()
        assert_array_almost_equal((asp*csp).todense(), a*c)
        assert_array_almost_equal( asp*c, a*c)

        assert_array_almost_equal( a*csp, a*c)
        assert_array_almost_equal( a2*csp, a*c)

        # Test provided by Andy Fraser, 2006-03-26
        L = 30
        frac = .3
        random.seed(0) # make runs repeatable
        A = zeros((L,2))
        for i in xrange(L):
            for j in xrange(2):
                r = random.random()
                if r < frac:
                    A[i,j] = r/frac

        A = self.spmatrix(A)
        B = A*A.T
        assert_array_almost_equal(B.todense(), A.todense() * A.T.todense())
        assert_array_almost_equal(B.todense(), A.todense() * A.todense().T)


        # check dimension mismatch  2x2 times 3x2
        A = self.spmatrix( [[1,2],[3,4]] )
        B = self.spmatrix( [[1,2],[3,4],[5,6]] )
        assert_raises(ValueError, A.__mul__, B)

    def test_matmat_dense(self):
        a = matrix([[3,0,0],[0,1,0],[2,0,3.0],[2,3,0]])
        asp = self.spmatrix(a)

        # check both array and matrix types
        bs = [ array([[1,2],[3,4],[5,6]]), matrix([[1,2],[3,4],[5,6]]) ]

        for b in bs:
            result = asp*b
            assert_( isinstance(result, type(b)) )
            assert_equal( result.shape, (4,2) )
            assert_equal( result, dot(a,b) )

    def test_sparse_format_conversions(self):
        A = sparse.kron( [[1,0,2],[0,3,4],[5,0,0]], [[1,2],[0,3]] )
        D = A.todense()
        A = self.spmatrix(A)

        for format in ['bsr','coo','csc','csr','dia','dok','lil']:
            a = A.asformat(format)
            assert_equal(a.format,format)
            assert_array_equal(a.todense(), D)

            b = self.spmatrix(D+3j).asformat(format)
            assert_equal(b.format,format)
            assert_array_equal(b.todense(), D+3j)

            c = eval(format + '_matrix')(A)
            assert_equal(c.format,format)
            assert_array_equal(c.todense(), D)


    def test_tobsr(self):
        x = array([[1,0,2,0],[0,0,0,0],[0,0,4,5]])
        y = array([[0,1,2],[3,0,5]])
        A = kron(x,y)
        Asp = self.spmatrix(A)
        for format in ['bsr']:
            fn = getattr(Asp, 'to' + format )

            for X in [ 1, 2, 3, 6 ]:
                for Y in [ 1, 2, 3, 4, 6, 12]:
                    assert_equal( fn(blocksize=(X,Y)).todense(), A)


    def test_transpose(self):
        a = self.datsp.transpose()
        b = self.dat.transpose()
        assert_array_equal(a.todense(), b)
        assert_array_equal(a.transpose().todense(), self.dat)

        assert_array_equal( self.spmatrix((3,4)).T.todense(), zeros((4,3)) )


    def test_add_dense(self):
        """ adding a dense matrix to a sparse matrix
        """
        sum1 = self.dat + self.datsp
        assert_array_equal(sum1, 2*self.dat)
        sum2 = self.datsp + self.dat
        assert_array_equal(sum2, 2*self.dat)

    def test_sub_dense(self):
        """ subtracting a dense matrix to/from a sparse matrix
        """
        sum1 = 3*self.dat - self.datsp
        assert_array_equal(sum1, 2*self.dat)
        sum2 = 3*self.datsp - self.dat
        assert_array_equal(sum2, 2*self.dat)


    def test_copy(self):
        """ Check whether the copy=True and copy=False keywords work
        """
        A = self.datsp

        #check that copy preserves format
        assert_equal(A.copy().format, A.format)
        assert_equal(A.__class__(A,copy=True).format,  A.format)
        assert_equal(A.__class__(A,copy=False).format, A.format)

        assert_equal(A.copy().todense(), A.todense())
        assert_equal(A.__class__(A,copy=True).todense(),  A.todense())
        assert_equal(A.__class__(A,copy=False).todense(), A.todense())

        #check that XXX_matrix.toXXX() works
        toself = getattr(A,'to' + A.format)
        assert_equal(toself().format, A.format)
        assert_equal(toself(copy=True).format, A.format)
        assert_equal(toself(copy=False).format, A.format)

        assert_equal(toself().todense(), A.todense())
        assert_equal(toself(copy=True).todense(), A.todense())
        assert_equal(toself(copy=False).todense(), A.todense())


        # check whether the data is copied?
        # TODO: deal with non-indexable types somehow
        B = A.copy()
        try:
            B[0,0] += 1
            assert_(B[0,0] != A[0,0])
        except NotImplementedError:
            # not all sparse matrices can be indexed
            pass
        except TypeError:
            # not all sparse matrices can be indexed
            pass

    # Eventually we'd like to allow matrix products between dense
    # and sparse matrices using the normal dot() function:
    #def test_dense_dot_sparse(self):
    #    a = array([1.,2.,3.])
    #    dense_dot_dense = dot(a, self.dat)
    #    dense_dot_sparse = dot(a, self.datsp)
    #    assert_array_equal(dense_dot_dense, dense_dot_sparse)

    #def test_sparse_dot_dense(self):
    #    b = array([1.,2.,3.,4.])
    #    dense_dot_dense = dot(self.dat, b)
    #    dense_dot_sparse = dot(self.datsp, b)
    #    assert_array_equal(dense_dot_dense, dense_dot_sparse)



class _TestInplaceArithmetic:
    def test_imul_scalar(self):
        a = self.datsp.copy()
        a *= 2
        assert_array_equal(self.dat*2,a.todense())

        a = self.datsp.copy()
        a *= 17.3
        assert_array_equal(self.dat*17.3,a.todense())

    def test_idiv_scalar(self):
        a = self.datsp.copy()
        a /= 2
        assert_array_equal(self.dat/2,a.todense())

        a = self.datsp.copy()
        a /= 17.3
        assert_array_equal(self.dat/17.3,a.todense())


class _TestGetSet:
    def test_setelement(self):
        A = self.spmatrix((3,4))
        A[ 0, 0] = 0 # bug 870
        A[ 1, 2] = 4.0
        A[ 0, 1] = 3
        A[ 2, 0] = 2.0
        A[ 0,-1] = 8
        A[-1,-2] = 7
        A[ 0, 1] = 5
        assert_array_equal(A.todense(),[[0,5,0,8],[0,0,4,0],[2,0,7,0]])

        for ij in [(0,4),(-1,4),(3,0),(3,4),(3,-1)]:
            assert_raises(IndexError, A.__setitem__, ij, 123.0)

        for v in [[1,2,3], array([1,2,3])]:
            assert_raises(ValueError, A.__setitem__, (0,0), v)

        for v in [3j]:
            assert_raises(TypeError, A.__setitem__, (0,0), v)

    def test_getelement(self):
        D = array([[1,0,0],
                   [4,3,0],
                   [0,2,0],
                   [0,0,0]])
        A = self.spmatrix(D)

        M,N = D.shape

        for i in range(-M, M):
            for j in range(-N, N):
                assert_equal(A[i,j], D[i,j])

        for ij in [(0,3),(-1,3),(4,0),(4,3),(4,-1)]:
            assert_raises(IndexError, A.__getitem__, ij)

class _TestSolve:
    def test_solve(self):
        """ Test whether the lu_solve command segfaults, as reported by Nils
        Wagner for a 64-bit machine, 02 March 2005 (EJS)
        """
        n = 20
        np.random.seed(0) #make tests repeatable
        A = zeros((n,n), dtype=complex)
        x = np.random.rand(n)
        y = np.random.rand(n-1)+1j*np.random.rand(n-1)
        r = np.random.rand(n)
        for i in range(len(x)):
            A[i,i] = x[i]
        for i in range(len(y)):
            A[i,i+1] = y[i]
            A[i+1,i] = conjugate(y[i])
        A = self.spmatrix(A)
        x = splu(A).solve(r)
        assert_almost_equal(A*x,r)


class _TestHorizSlicing:
    """Tests horizontal slicing (e.g. [0, :]).  Tests for individual sparse
    matrix types that implement this should derive from this class.
    """
    def test_get_horiz_slice(self):
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
        assert_(caught == 2)


class _TestVertSlicing:
    """Tests vertical slicing (e.g. [:, 0]).  Tests for individual sparse
    matrix types that implement this should derive from this class.
    """
    def test_get_vert_slice(self):
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
        assert_(caught == 2)



class _TestBothSlicing:
    """Tests vertical and horizontal slicing (e.g. [:,0:2]). Tests for
    individual sparse matrix types that implement this should derive from this
    class.
    """
    def test_get_slices(self):
        B = asmatrix(arange(50.).reshape(5,10))
        A = self.spmatrix(B)
        assert_array_equal(A[2:5,0:3].todense(), B[2:5,0:3])
        assert_array_equal(A[1:,:-1].todense(),  B[1:,:-1])
        assert_array_equal(A[:-1,1:].todense(),  B[:-1,1:])

        # Now test slicing when a column contains only zeros
        E = matrix([[1, 0, 1], [4, 0, 0], [0, 0, 0], [0, 0, 1]])
        F = self.spmatrix(E)
        assert_array_equal(E[1:2, 1:2], F[1:2, 1:2].todense())
        assert_array_equal(E[:, 1:], F[:, 1:].todense())

class _TestFancyIndexing:
    """Tests fancy indexing features.  The tests for any matrix formats
    that implement these features should derive from this class.
    """
    def test_fancy_indexing_set(self):
        n, m = (5, 10)
        def _test_set(i, j, nitems):
            A = self.spmatrix((n, m))
            A[i, j] = 1
            assert_almost_equal(A.sum(), nitems)
            assert_almost_equal(A[i, j], 1)

        # [i,j]
        for i, j in [(2, 3), (-1, 8), (-1, -2), (array(-1), -2), (-1, array(-2)),
                     (array(-1), array(-2))]:
            _test_set(i, j, 1)

        # [i,1:2]
        for i, j in [(2, slice(m)), (2, slice(5, -2)), (array(2), slice(5, -2))]:
            _test_set(i, j, 3)

    def test_fancy_indexing(self):
        B = asmatrix(arange(50).reshape(5,10))
        A = self.spmatrix( B )

        # [i,j]
        assert_equal(A[2,3],  B[2,3])
        assert_equal(A[-1,8], B[-1,8])
        assert_equal(A[-1,-2],B[-1,-2])
        assert_equal(A[array(-1),-2],B[-1,-2])
        assert_equal(A[-1,array(-2)],B[-1,-2])
        assert_equal(A[array(-1),array(-2)],B[-1,-2])

        # [i,1:2]
        assert_equal(A[2,:].todense(),   B[2,:])
        assert_equal(A[2,5:-2].todense(),B[2,5:-2])
        assert_equal(A[array(2),5:-2].todense(),B[2,5:-2])

        # [i,[1,2]]
        assert_equal(A[3,[1,3]].todense(),  B[3,[1,3]])
        assert_equal(A[-1,[2,-5]].todense(),B[-1,[2,-5]])
        assert_equal(A[array(-1),[2,-5]].todense(),B[-1,[2,-5]])
        assert_equal(A[-1,array([2,-5])].todense(),B[-1,[2,-5]])
        assert_equal(A[array(-1),array([2,-5])].todense(),B[-1,[2,-5]])

        # [1:2,j]
        assert_equal(A[:,2].todense(),   B[:,2])
        assert_equal(A[3:4,9].todense(), B[3:4,9])
        assert_equal(A[1:4,-5].todense(),B[1:4,-5])
        assert_equal(A[2:-1,3].todense(),B[2:-1,3])
        assert_equal(A[2:-1,array(3)].todense(),B[2:-1,3])

        # [1:2,1:2]
        assert_equal(A[1:2,1:2].todense(),B[1:2,1:2])
        assert_equal(A[4:,3:].todense(),  B[4:,3:])
        assert_equal(A[:4,:5].todense(),  B[:4,:5])
        assert_equal(A[2:-1,:5].todense(),B[2:-1,:5])

        # [1:2,[1,2]]
        assert_equal(A[:,[2,8,3,-1]].todense(),B[:,[2,8,3,-1]])
        assert_equal(A[3:4,[9]].todense(),     B[3:4,[9]])
        assert_equal(A[1:4,[-1,-5]].todense(), B[1:4,[-1,-5]])
        assert_equal(A[1:4,array([-1,-5])].todense(), B[1:4,[-1,-5]])

        # [[1,2],j]
        assert_equal(A[[1,3],3].todense(),   B[[1,3],3])
        assert_equal(A[[2,-5],-4].todense(), B[[2,-5],-4])
        assert_equal(A[array([2,-5]),-4].todense(), B[[2,-5],-4])
        assert_equal(A[[2,-5],array(-4)].todense(), B[[2,-5],-4])
        assert_equal(A[array([2,-5]),array(-4)].todense(), B[[2,-5],-4])

        # [[1,2],1:2]
        assert_equal(A[[1,3],:].todense(),    B[[1,3],:])
        assert_equal(A[[2,-5],8:-1].todense(),B[[2,-5],8:-1])
        assert_equal(A[array([2,-5]),8:-1].todense(),B[[2,-5],8:-1])

        # [[1,2],[1,2]]
        assert_equal(A[[1,3],[2,4]],   B[[1,3],[2,4]])
        assert_equal(A[[-1,-3],[2,-4]],B[[-1,-3],[2,-4]])
        assert_equal(A[array([-1,-3]),[2,-4]],B[[-1,-3],[2,-4]])
        assert_equal(A[[-1,-3],array([2,-4])],B[[-1,-3],[2,-4]])
        assert_equal(A[array([-1,-3]),array([2,-4])],B[[-1,-3],[2,-4]])

        # [[[1],[2]],[1,2]]
        assert_equal(A[[[1],[3]],[2,4]].todense(),        B[[[1],[3]],[2,4]])
        assert_equal(A[[[-1],[-3],[-2]],[2,-4]].todense(),B[[[-1],[-3],[-2]],[2,-4]])
        assert_equal(A[array([[-1],[-3],[-2]]),[2,-4]].todense(),B[[[-1],[-3],[-2]],[2,-4]])
        assert_equal(A[[[-1],[-3],[-2]],array([2,-4])].todense(),B[[[-1],[-3],[-2]],[2,-4]])
        assert_equal(A[array([[-1],[-3],[-2]]),array([2,-4])].todense(),B[[[-1],[-3],[-2]],[2,-4]])

        # [i]
        assert_equal(A[1,:].todense(), B[1,:])
        assert_equal(A[-2,:].todense(),B[-2,:])
        assert_equal(A[array(-2),:].todense(),B[-2,:])

        # [1:2]
        assert_equal(A[1:4].todense(), B[1:4])
        assert_equal(A[1:-2].todense(),B[1:-2])

        # [[1,2]]
        assert_equal(A[[1,3]].todense(),  B[[1,3]])
        assert_equal(A[[-1,-3]].todense(),B[[-1,-3]])
        assert_equal(A[array([-1,-3])].todense(),B[[-1,-3]])

        # [[1,2],:][:,[1,2]]
        assert_equal(A[[1,3],:][:,[2,4]].todense(),    B[[1,3],:][:,[2,4]]    )
        assert_equal(A[[-1,-3],:][:,[2,-4]].todense(), B[[-1,-3],:][:,[2,-4]] )
        assert_equal(A[array([-1,-3]),:][:,array([2,-4])].todense(), B[[-1,-3],:][:,[2,-4]] )

        # [:,[1,2]][[1,2],:]
        assert_equal(A[:,[1,3]][[2,4],:].todense(),    B[:,[1,3]][[2,4],:]    )
        assert_equal(A[:,[-1,-3]][[2,-4],:].todense(), B[:,[-1,-3]][[2,-4],:] )
        assert_equal(A[:,array([-1,-3])][array([2,-4]),:].todense(), B[:,[-1,-3]][[2,-4],:] )


        # Check bug reported by Robert Cimrman:
        # http://thread.gmane.org/gmane.comp.python.scientific.devel/7986
        s = slice(int8(2),int8(4),None)
        assert_equal(A[s,:].todense(), B[2:4,:])
        assert_equal(A[:,s].todense(), B[:,2:4])

    def test_fancy_indexing_randomized(self):
        random.seed(0) # make runs repeatable

        NUM_SAMPLES = 50
        M = 6
        N = 4

        D = np.asmatrix(np.random.rand(M,N))
        D = np.multiply(D, D > 0.5)

        I = np.random.random_integers(-M + 1, M - 1, size=NUM_SAMPLES)
        J = np.random.random_integers(-N + 1, N - 1, size=NUM_SAMPLES)

        S = self.spmatrix(D)

        assert_equal(S[I,J], D[I,J])

        I_bad = I + M
        J_bad = J - N

        assert_raises(IndexError, S.__getitem__, (I_bad,J))
        assert_raises(IndexError, S.__getitem__, (I,J_bad))


class _TestArithmetic:
    """
    Test real/complex arithmetic
    """
    def arith_init(self):
        #these can be represented exactly in FP (so arithmetic should be exact)
        self.A = matrix([[   -1.5,    6.5,       0,    2.25,  0,  0],
                         [  3.125, -7.875,   0.625,       0,  0,  0],
                         [      0,      0,  -0.125,     1.0,  0,  0],
                         [      0,      0,   8.375,       0,  0,  0]],'float64')
        self.B = matrix([[  0.375,       0,    0,   0,      -5,     2.5],
                         [  14.25,   -3.75,    0,   0,  -0.125,       0],
                         [      0,    7.25,    0,   0,       0,       0],
                         [   18.5, -0.0625,    0,   0,       0,       0]],'complex128')
        self.B.imag = matrix([[    1.25,     0,   0,   0,  6, -3.875],
                              [    2.25, 4.125,   0,   0,  0,   2.75],
                              [       0, 4.125,   0,   0,  0,      0],
                              [ -0.0625,     0,   0,   0,  0,      0]],'float64')

        #fractions are all x/16ths
        assert_array_equal((self.A*16).astype('int32'),16*self.A)
        assert_array_equal((self.B.real*16).astype('int32'),16*self.B.real)
        assert_array_equal((self.B.imag*16).astype('int32'),16*self.B.imag)

        self.Asp = self.spmatrix(self.A)
        self.Bsp = self.spmatrix(self.B)

    def test_add_sub(self):
        self.arith_init()

        #basic tests
        assert_array_equal((self.Asp+self.Bsp).todense(),self.A+self.B)

        #check conversions
        for x in supported_dtypes:
            A   = self.A.astype(x)
            Asp = self.spmatrix(A)
            for y in supported_dtypes:
                B   = self.B.astype(y)
                Bsp = self.spmatrix(B)

                #addition
                D1 = A + B
                S1 = Asp + Bsp

                assert_equal(S1.dtype,D1.dtype)
                assert_array_equal(S1.todense(),D1)
                assert_array_equal(Asp + B,D1)          #check sparse + dense
                assert_array_equal(A + Bsp,D1)          #check dense + sparse

                #subtraction
                D1 = A - B
                S1 = Asp - Bsp

                assert_equal(S1.dtype,D1.dtype)
                assert_array_equal(S1.todense(),D1)
                assert_array_equal(Asp - B,D1)          #check sparse - dense
                assert_array_equal(A - Bsp,D1)          #check dense - sparse


    def test_mu(self):
        self.arith_init()

        #basic tests
        assert_array_equal((self.Asp*self.Bsp.T).todense(),self.A*self.B.T)

        for x in supported_dtypes:
            A   = self.A.astype(x)
            Asp = self.spmatrix(A)
            for y in supported_dtypes:
                B   = self.B.astype(y)
                Bsp = self.spmatrix(B)

                D1 = A * B.T
                S1 = Asp * Bsp.T

                assert_array_equal(S1.todense(),D1)
                assert_equal(S1.dtype,D1.dtype)



class TestCSR(_TestCommon, _TestGetSet, _TestSolve,
        _TestInplaceArithmetic, _TestArithmetic,
        _TestHorizSlicing, _TestVertSlicing, _TestBothSlicing,
        _TestFancyIndexing, TestCase):
    spmatrix = csr_matrix

    @dec.knownfailureif(True, "Fancy indexing is known to be broken for CSR" \
                              " matrices")
    def test_fancy_indexing_set(self):
        _TestFancyIndexing.test_fancy_indexing_set(self)

    def test_constructor1(self):
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

    def test_constructor2(self):
        b = zeros((6,6),'d')
        b[3,4] = 5
        bsp = csr_matrix(b)
        assert_array_almost_equal(bsp.data,[5])
        assert_array_equal(bsp.indices,[4])
        assert_array_equal(bsp.indptr,[0,0,0,0,1,1,1])
        assert_array_almost_equal(bsp.todense(),b)

    def test_constructor3(self):
        b = matrix([[1,0],
                   [0,2],
                   [3,0]],'d')
        bsp = csr_matrix(b)
        assert_array_almost_equal(bsp.data,[1,2,3])
        assert_array_equal(bsp.indices,[0,1,0])
        assert_array_equal(bsp.indptr,[0,1,2,3])
        assert_array_almost_equal(bsp.todense(),b)

### currently disabled
##    def test_constructor4(self):
##        """try using int64 indices"""
##        data = arange( 6 ) + 1
##        col = array( [1, 2, 1, 0, 0, 2], dtype='int64' )
##        ptr = array( [0, 2, 4, 6], dtype='int64' )
##
##        a = csr_matrix( (data, col, ptr), shape = (3,3) )
##
##        b = matrix([[0,1,2],
##                    [4,3,0],
##                    [5,0,6]],'d')
##
##        assert_equal(a.indptr.dtype,numpy.dtype('int64'))
##        assert_equal(a.indices.dtype,numpy.dtype('int64'))
##        assert_array_equal(a.todense(),b)

    def test_constructor4(self):
        """using (data, ij) format"""
        row  = array([2, 3, 1, 3, 0, 1, 3, 0, 2, 1, 2])
        col  = array([0, 1, 0, 0, 1, 1, 2, 2, 2, 2, 1])
        data = array([  6.,  10.,   3.,   9.,   1.,   4.,
                              11.,   2.,   8.,   5.,   7.])

        ij = vstack((row,col))
        csr = csr_matrix((data,ij),(4,3))
        assert_array_equal(arange(12).reshape(4,3),csr.todense())

    def test_constructor5(self):
        """infer dimensions from arrays"""
        indptr  = array([0,1,3,3])
        indices = array([0,5,1,2])
        data    = array([1,2,3,4])
        csr = csr_matrix((data, indices, indptr))
        assert_array_equal(csr.shape,(3,6))

    def test_sort_indices(self):
        data    = arange( 5 )
        indices = array( [7, 2, 1, 5, 4] )
        indptr  = array( [0, 3, 5] )
        asp = csr_matrix( (data, indices, indptr), shape=(2,10) )
        bsp = asp.copy()
        asp.sort_indices( )
        assert_array_equal(asp.indices,[1, 2, 7, 4, 5])
        assert_array_equal(asp.todense(),bsp.todense())

    def test_eliminate_zeros(self):
        data    = array( [1, 0, 0, 0, 2, 0, 3, 0] )
        indices = array( [1, 2, 3, 4, 5, 6, 7, 8] )
        indptr  = array( [0, 3, 8] )
        asp = csr_matrix( (data, indices, indptr), shape=(2,10) )
        bsp = asp.copy()
        asp.eliminate_zeros( )
        assert_array_equal(asp.nnz, 3)
        assert_array_equal(asp.data,[1, 2, 3])
        assert_array_equal(asp.todense(),bsp.todense())

    def test_unsorted_arithmetic(self):
        data    = arange( 5 )
        indices = array( [7, 2, 1, 5, 4] )
        indptr  = array( [0, 3, 5] )
        asp = csr_matrix( (data, indices, indptr), shape=(2,10) )
        data    = arange( 6 )
        indices = array( [8, 1, 5, 7, 2, 4] )
        indptr  = array( [0, 2, 6] )
        bsp = csr_matrix( (data, indices, indptr), shape=(2,10) )
        assert_equal((asp + bsp).todense(), asp.todense() + bsp.todense())




class TestCSC(_TestCommon, _TestGetSet, _TestSolve,
        _TestInplaceArithmetic, _TestArithmetic,
        _TestHorizSlicing, _TestVertSlicing, _TestBothSlicing,
        _TestFancyIndexing, TestCase):
    spmatrix = csc_matrix

    @dec.knownfailureif(True, "Fancy indexing is known to be broken for CSC" \
                              " matrices")
    def test_fancy_indexing_set(self):
        _TestFancyIndexing.test_fancy_indexing_set(self)

    def test_constructor1(self):
        b = matrix([[1,0,0,0],[0,0,1,0],[0,2,0,3]],'d')
        bsp = csc_matrix(b)
        assert_array_almost_equal(bsp.data,[1,2,1,3])
        assert_array_equal(bsp.indices,[0,2,1,2])
        assert_array_equal(bsp.indptr,[0,1,2,3,4])
        assert_equal(bsp.getnnz(),4)
        assert_equal(bsp.shape,b.shape)
        assert_equal(bsp.getformat(),'csc')

    def test_constructor2(self):
        b = zeros((6,6),'d')
        b[2,4] = 5
        bsp = csc_matrix(b)
        assert_array_almost_equal(bsp.data,[5])
        assert_array_equal(bsp.indices,[2])
        assert_array_equal(bsp.indptr,[0,0,0,0,0,1,1])

    def test_constructor3(self):
        b = matrix([[1,0],[0,0],[0,2]],'d')
        bsp = csc_matrix(b)
        assert_array_almost_equal(bsp.data,[1,2])
        assert_array_equal(bsp.indices,[0,2])
        assert_array_equal(bsp.indptr,[0,1,2])

    def test_constructor4(self):
        """using (data, ij) format"""
        row  = array([2, 3, 1, 3, 0, 1, 3, 0, 2, 1, 2])
        col  = array([0, 1, 0, 0, 1, 1, 2, 2, 2, 2, 1])
        data = array([  6.,  10.,   3.,   9.,   1.,   4.,
                              11.,   2.,   8.,   5.,   7.])

        ij = vstack((row,col))
        csc = csc_matrix((data,ij),(4,3))
        assert_array_equal(arange(12).reshape(4,3),csc.todense())

    def test_constructor5(self):
        """infer dimensions from arrays"""
        indptr  = array([0,1,3,3])
        indices = array([0,5,1,2])
        data    = array([1,2,3,4])
        csc = csc_matrix((data, indices, indptr))
        assert_array_equal(csc.shape,(6,3))

    def test_eliminate_zeros(self):
        data    = array( [1, 0, 0, 0, 2, 0, 3, 0] )
        indices = array( [1, 2, 3, 4, 5, 6, 7, 8] )
        indptr  = array( [0, 3, 8] )
        asp = csc_matrix( (data, indices, indptr), shape=(10,2) )
        bsp = asp.copy()
        asp.eliminate_zeros( )
        assert_array_equal(asp.nnz, 3)
        assert_array_equal(asp.data,[1, 2, 3])
        assert_array_equal(asp.todense(),bsp.todense())

    def test_sort_indices(self):
        data = arange( 5 )
        row = array( [7, 2, 1, 5, 4] )
        ptr = [0, 3, 5]
        asp = csc_matrix( (data, row, ptr), shape=(10,2) )
        bsp = asp.copy()
        asp.sort_indices()
        assert_array_equal(asp.indices,[1, 2, 7, 4, 5])
        assert_array_equal(asp.todense(),bsp.todense())

    def test_unsorted_arithmetic(self):
        data    = arange( 5 )
        indices = array( [7, 2, 1, 5, 4] )
        indptr  = array( [0, 3, 5] )
        asp = csc_matrix( (data, indices, indptr), shape=(10,2) )
        data    = arange( 6 )
        indices = array( [8, 1, 5, 7, 2, 4] )
        indptr  = array( [0, 2, 6] )
        bsp = csc_matrix( (data, indices, indptr), shape=(10,2) )
        assert_equal((asp + bsp).todense(), asp.todense() + bsp.todense())

class TestDOK(_TestCommon, _TestGetSet, _TestSolve, TestCase):
    spmatrix = dok_matrix

    def test_mult(self):
        A = dok_matrix((10,10))
        A[0,3] = 10
        A[5,6] = 20
        D = A*A.T
        E = A*A.H
        assert_array_equal(D.A, E.A)

    def test_add(self):
        A = dok_matrix((3,2))
        A[0,1] = -10
        A[2,0] = 20
        A = A + 10
        B = matrix([[10, 0], [10, 10], [30, 10]])
        assert_array_equal(A.todense(), B)

    def test_convert(self):
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

    def test_set_slice(self):
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
        except ValueError:
            caught += 1
        try:
            A[0,0] = arange(100)
        except ValueError:
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

    def test_ctor(self):
        caught = 0
        # Empty ctor
        try:
            A = dok_matrix()
        except TypeError, e:
            caught+=1
        assert_equal(caught, 1)

        # Dense ctor
        b = matrix([[1,0,0,0],[0,0,1,0],[0,2,0,3]],'d')
        A = dok_matrix(b)
        assert_equal(A.todense(), b)

        # Sparse ctor
        c = csr_matrix(b)
        assert_equal(A.todense(), c.todense())

    def test_resize(self):
        """A couple basic tests of the resize() method.

        resize(shape) resizes the array in-place.
        """
        a = dok_matrix((5,5))
        a[:,0] = 1
        a.resize((2,2))
        expected1 = array([[1,0],[1,0]])
        assert_array_equal(a.todense(), expected1)
        a.resize((3,2))
        expected2 = array([[1,0],[1,0],[0,0]])
        assert_array_equal(a.todense(), expected2)


    def test_ticket1160(self):
        """Regression test for ticket #1160."""
        a = dok_matrix((3,3))
        a[0,0] = 0
        # This assert would fail, because the above assignment would
        # incorrectly call __set_item__ even though the value was 0.
        assert_((0,0) not in a.keys(), "Unexpected entry (0,0) in keys")

        # Slice assignments were also affected.
        b = dok_matrix((3,3))
        b[:,0] = 0
        assert_(len(b.keys())==0, "Unexpected entries in keys")

    # The following five tests are duplicates from _TestCommon, so they can be
    # marked as knownfail for Python 2.4.  Once 2.4 is no longer supported,
    # these duplicates can be removed again.

    @dec.knownfailureif(sys.version[:3] == '2.4', "See ticket 1559")
    def test_add_dense(self):
        """ adding a dense matrix to a sparse matrix
        """
        sum1 = self.dat + self.datsp
        assert_array_equal(sum1, 2*self.dat)
        sum2 = self.datsp + self.dat
        assert_array_equal(sum2, 2*self.dat)

    @dec.knownfailureif(sys.version[:3] == '2.4', "See ticket 1559")
    def test_radd(self):
        a = self.dat.copy()
        a[0,2] = 2.0
        b = self.datsp
        c = a + b
        assert_array_equal(c,[[2,0,2,4],[6,0,2,0],[0,4,0,0]])

    @dec.knownfailureif(sys.version[:3] == '2.4', "See ticket 1559")
    def test_rsub(self):
        assert_array_equal((self.dat - self.datsp),[[0,0,0,0],[0,0,0,0],[0,0,0,0]])
        assert_array_equal((self.datsp - self.dat),[[0,0,0,0],[0,0,0,0],[0,0,0,0]])

        A = self.spmatrix(matrix([[1,0,0,4],[-1,0,0,0],[0,8,0,-5]],'d'))
        assert_array_equal((self.dat - A),self.dat - A.todense())
        assert_array_equal((A - self.dat),A.todense() - self.dat)
        assert_array_equal(A.todense() - self.datsp,A.todense() - self.dat)
        assert_array_equal(self.datsp - A.todense(),self.dat - A.todense())

    @dec.knownfailureif(sys.version[:3] == '2.4', "See ticket 1559")
    def test_matmat_sparse(self):
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
        assert_array_almost_equal( asp*c, a*c)

        assert_array_almost_equal( a*csp, a*c)
        assert_array_almost_equal( a2*csp, a*c)
        csp = bsp.tocsr()
        assert_array_almost_equal((asp*csp).todense(), a*c)
        assert_array_almost_equal( asp*c, a*c)

        assert_array_almost_equal( a*csp, a*c)
        assert_array_almost_equal( a2*csp, a*c)
        csp = bsp.tocoo()
        assert_array_almost_equal((asp*csp).todense(), a*c)
        assert_array_almost_equal( asp*c, a*c)

        assert_array_almost_equal( a*csp, a*c)
        assert_array_almost_equal( a2*csp, a*c)

        # Test provided by Andy Fraser, 2006-03-26
        L = 30
        frac = .3
        random.seed(0) # make runs repeatable
        A = zeros((L,2))
        for i in xrange(L):
            for j in xrange(2):
                r = random.random()
                if r < frac:
                    A[i,j] = r/frac

        A = self.spmatrix(A)
        B = A*A.T
        assert_array_almost_equal(B.todense(), A.todense() * A.T.todense())
        assert_array_almost_equal(B.todense(), A.todense() * A.todense().T)

        # check dimension mismatch  2x2 times 3x2
        A = self.spmatrix( [[1,2],[3,4]] )
        B = self.spmatrix( [[1,2],[3,4],[5,6]] )
        assert_raises(ValueError, A.__mul__, B)

    @dec.knownfailureif(sys.version[:3] == '2.4', "See ticket 1559")
    def test_sub_dense(self):
        """ subtracting a dense matrix to/from a sparse matrix
        """
        sum1 = 3*self.dat - self.datsp
        assert_array_equal(sum1, 2*self.dat)
        sum2 = 3*self.datsp - self.dat
        assert_array_equal(sum2, 2*self.dat)


class TestLIL( _TestCommon, _TestHorizSlicing, _TestVertSlicing,
        _TestBothSlicing, _TestGetSet, _TestSolve,
        _TestArithmetic, _TestInplaceArithmetic, _TestFancyIndexing,
        TestCase):
    spmatrix = lil_matrix

    B = lil_matrix((4,3))
    B[0,0] = 2
    B[1,2] = 7
    B[2,1] = 3
    B[3,0] = 10


    @dec.knownfailureif(True, "Fancy indexing is known to be broken for LIL" \
                              " matrices")
    def test_fancy_indexing_set(self):
        _TestFancyIndexing.test_fancy_indexing_set(self)

    @dec.knownfailureif(True, "Fancy indexing is known to be broken for LIL" \
                              " matrices")
    def test_fancy_indexing_randomized(self):
        _TestFancyIndexing.test_fancy_indexing_randomized(self)

    def test_dot(self):
        A = matrix(zeros((10,10)))
        A[0,3] = 10
        A[5,6] = 20

        B = lil_matrix((10,10))
        B[0,3] = 10
        B[5,6] = 20
        assert_array_equal(A * A.T, (B * B.T).todense())
        assert_array_equal(A * A.H, (B * B.H).todense())

    def test_scalar_mul(self):
        x = lil_matrix((3,3))
        x[0,0] = 2

        x = x*2
        assert_equal(x[0,0],4)

        x = x*0
        assert_equal(x[0,0],0)

    def test_reshape(self):
        x = lil_matrix((4,3))
        x[0,0] = 1
        x[2,1] = 3
        x[3,2] = 5
        x[0,2] = 7

        for s in [(12,1),(1,12)]:
            assert_array_equal(x.reshape(s).todense(),
                               x.todense().reshape(s))

    def test_lil_lil_assignment(self):
        """ Tests whether a row of one lil_matrix can be assigned to
        another.
        """
        B = self.B.copy()
        A = B / 10
        B[0,:] = A[0,:]
        assert_array_equal(A[0,:].A, B[0,:].A)


    def test_inplace_ops(self):
        A = lil_matrix([[0,2,3],[4,0,6]])
        B = lil_matrix([[0,1,0],[0,2,3]])

        data = {'add': (B,A + B),
                'sub': (B,A - B),
                'mul': (3,A * 3)}

        for op,(other,expected) in data.iteritems():
            result = A.copy()
            getattr(result, '__i%s__' % op)(other)

            assert_array_equal(result.todense(), expected.todense())


    def test_lil_slice_assignment(self):
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

    def test_lil_sequence_assignment(self):
        A = lil_matrix((4,3))
        B = eye(3,4,format='lil')

        i0 = [0,1,2]
        i1 = (0,1,2)
        i2 = array( i0 )

        A[0,i0] = B[i0,0]
        A[1,i1] = B[i1,1]
        A[2,i2] = B[i2,2]
        assert_array_equal(A.todense(),B.T.todense())

        # column slice
        A = lil_matrix((2,3))
        A[1,1:3] = [10,20]
        assert_array_equal(A.todense(), [[0,0,0],[0,10,20]])

        # column slice
        A = lil_matrix((3,2))
        A[1:3,1] = [[10],[20]]
        assert_array_equal(A.todense(), [[0,0],[0,10],[0,20]])

    def test_lil_iteration(self):
        row_data = [[1,2,3],[4,5,6]]
        B = lil_matrix(array(row_data))
        for r,row in enumerate(B):
            assert_array_equal(row.todense(),array(row_data[r],ndmin=2))

    def test_lil_from_csr(self):
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

    def test_fancy_indexing(self):
        M = arange(25).reshape(5,5)
        A = lil_matrix( M )

        assert_equal(A[array([1,2,3]),2:3].todense(), M[array([1,2,3]),2:3])

    def test_point_wise_multiply(self):
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



class TestCOO(_TestCommon, TestCase):
    spmatrix = coo_matrix
    def test_constructor1(self):
        """unsorted triplet format"""
        row  = array([2, 3, 1, 3, 0, 1, 3, 0, 2, 1, 2])
        col  = array([0, 1, 0, 0, 1, 1, 2, 2, 2, 2, 1])
        data = array([  6.,  10.,   3.,   9.,   1.,   4.,
                              11.,   2.,   8.,   5.,   7.])

        coo = coo_matrix((data,(row,col)),(4,3))

        assert_array_equal(arange(12).reshape(4,3),coo.todense())

    def test_constructor2(self):
        """unsorted triplet format with duplicates (which are summed)"""
        row  = array([0,1,2,2,2,2,0,0,2,2])
        col  = array([0,2,0,2,1,1,1,0,0,2])
        data = array([2,9,-4,5,7,0,-1,2,1,-5])
        coo = coo_matrix((data,(row,col)),(3,3))

        mat = matrix([[4,-1,0],[0,0,9],[-3,7,0]])

        assert_array_equal(mat,coo.todense())

    def test_constructor3(self):
        """empty matrix"""
        coo = coo_matrix( (4,3) )

        assert_array_equal(coo.shape,(4,3))
        assert_array_equal(coo.row,[])
        assert_array_equal(coo.col,[])
        assert_array_equal(coo.data,[])
        assert_array_equal(coo.todense(),zeros((4,3)))

    def test_constructor4(self):
        """from dense matrix"""
        mat = array([[0,1,0,0],
                     [7,0,3,0],
                     [0,4,0,0]])
        coo = coo_matrix(mat)
        assert_array_equal(coo.todense(),mat)

        #upgrade rank 1 arrays to row matrix
        mat = array([0,1,0,0])
        coo = coo_matrix(mat)
        assert_array_equal(coo.todense(),mat.reshape(1,-1))


class TestDIA(_TestCommon, _TestArithmetic, TestCase):
    spmatrix = dia_matrix

    def test_constructor1(self):
        D = matrix([[1, 0, 3, 0],
                    [1, 2, 0, 4],
                    [0, 2, 3, 0],
                    [0, 0, 3, 4]])
        data    = np.array([[1,2,3,4]]).repeat(3,axis=0)
        offsets = np.array([0,-1,2])
        assert_equal(dia_matrix( (data,offsets), shape=(4,4)).todense(), D)



class TestBSR(_TestCommon, _TestArithmetic, _TestInplaceArithmetic, TestCase):
    spmatrix = bsr_matrix

    def test_constructor1(self):
        """check native BSR format constructor"""
        indptr  = array([0,2,2,4])
        indices = array([0,2,2,3])
        data    = zeros((4,2,3))

        data[0] = array([[ 0,  1,  2],
                         [ 3,  0,  5]])
        data[1] = array([[ 0,  2,  4],
                         [ 6,  0, 10]])
        data[2] = array([[ 0,  4,  8],
                         [12,  0, 20]])
        data[3] = array([[ 0,  5, 10],
                         [15,  0, 25]])

        A = kron( [[1,0,2,0],[0,0,0,0],[0,0,4,5]], [[0,1,2],[3,0,5]] )
        Asp = bsr_matrix((data,indices,indptr),shape=(6,12))
        assert_equal(Asp.todense(),A)

        #infer shape from arrays
        Asp = bsr_matrix((data,indices,indptr))
        assert_equal(Asp.todense(),A)

    def test_constructor2(self):
        """construct from dense"""

        #test zero mats
        for shape in [ (1,1), (5,1), (1,10), (10,4), (3,7), (2,1)]:
            A = zeros(shape)
            assert_equal(bsr_matrix(A).todense(),A)
        A = zeros((4,6))
        assert_equal(bsr_matrix(A,blocksize=(2,2)).todense(),A)
        assert_equal(bsr_matrix(A,blocksize=(2,3)).todense(),A)

        A = kron( [[1,0,2,0],[0,0,0,0],[0,0,4,5]], [[0,1,2],[3,0,5]] )
        assert_equal(bsr_matrix(A).todense(),A)
        assert_equal(bsr_matrix(A,shape=(6,12)).todense(),A)
        assert_equal(bsr_matrix(A,blocksize=(1,1)).todense(),A)
        assert_equal(bsr_matrix(A,blocksize=(2,3)).todense(),A)
        assert_equal(bsr_matrix(A,blocksize=(2,6)).todense(),A)
        assert_equal(bsr_matrix(A,blocksize=(2,12)).todense(),A)
        assert_equal(bsr_matrix(A,blocksize=(3,12)).todense(),A)
        assert_equal(bsr_matrix(A,blocksize=(6,12)).todense(),A)

        A = kron( [[1,0,2,0],[0,1,0,0],[0,0,0,0]], [[0,1,2],[3,0,5]] )
        assert_equal(bsr_matrix(A,blocksize=(2,3)).todense(),A)

    def test_eliminate_zeros(self):
        data = kron([1, 0, 0, 0, 2, 0, 3, 0], [[1,1],[1,1]]).T
        data = data.reshape(-1,2,2)
        indices = array( [1, 2, 3, 4, 5, 6, 7, 8] )
        indptr  = array( [0, 3, 8] )
        asp = bsr_matrix( (data, indices, indptr), shape=(4,20) )
        bsp = asp.copy()
        asp.eliminate_zeros()
        assert_array_equal(asp.nnz, 3*4)
        assert_array_equal(asp.todense(),bsp.todense())

    def test_bsr_matvec(self):
        A = bsr_matrix( arange(2*3*4*5).reshape(2*4,3*5), blocksize=(4,5) )
        x = arange(A.shape[1]).reshape(-1,1)
        assert_equal(A*x, A.todense()*x)

    def test_bsr_matvecs(self):
        A = bsr_matrix( arange(2*3*4*5).reshape(2*4,3*5), blocksize=(4,5) )
        x = arange(A.shape[1]*6).reshape(-1,6)
        assert_equal(A*x, A.todense()*x)


if __name__ == "__main__":
    run_module_suite()
