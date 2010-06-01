#!/usr/bin/env python
#
# Created by: Pearu Peterson, March 2002
#
""" Test functions for linalg.basic module

"""
"""
Bugs:
1) solve.check_random_sym_complex fails if a is complex
   and transpose(a) = conjugate(a) (a is Hermitian).
"""
__usage__ = """
Build linalg:
  python setup_linalg.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.linalg.test()'
Run tests if linalg is not installed:
  python tests/test_basic.py
"""

import warnings

from numpy import arange, array, dot, zeros, identity, conjugate, transpose, \
        float32, zeros_like
import numpy.linalg as linalg

from numpy.testing import *

from scipy.linalg import solve, inv, det, lstsq, pinv, pinv2, norm,\
        solve_banded, solveh_banded, cholesky_banded


def random(size):
    return rand(*size)


class TestSolveBanded(TestCase):

    def test_real(self):
        a = array([[ 1.0, 20,  0,  0],
                   [ -30,  4,  6,  0],
                   [   2,  1, 20,  2],
                   [   0, -1,  7, 14]])
        ab = array([[ 0.0, 20,  6,  2],
                    [   1,  4, 20, 14],
                    [ -30,  1,  7,  0],
                    [   2, -1,  0,  0]])
        l,u = 2,1
        b4 = array([10.0, 0.0, 2.0, 14.0])
        b4by1 = b4.reshape(-1,1)
        b4by2 = array([[ 2, 1],
                       [-30, 4],
                       [  2, 3],
                       [  1, 3]])
        b4by4 = array([[1, 0, 0, 0],
                       [0, 0, 0, 1],
                       [0, 1, 0, 0],
                       [0, 1, 0, 0]])
        for b in [b4, b4by1, b4by2, b4by4]:
            x = solve_banded((l, u), ab, b)
            assert_array_almost_equal(dot(a, x), b)

    def test_complex(self):
        a = array([[ 1.0, 20,  0,   0],
                   [ -30,  4,  6,   0],
                   [  2j,  1, 20,  2j],
                   [   0, -1,  7,  14]])
        ab = array([[ 0.0, 20,  6,  2j],
                    [   1,  4, 20,  14],
                    [ -30,  1,  7,   0],
                    [  2j, -1,  0,   0]])
        l,u = 2,1
        b4 = array([10.0, 0.0, 2.0, 14.0j])
        b4by1 = b4.reshape(-1,1)
        b4by2 = array([[ 2, 1],
                       [-30, 4],
                       [  2, 3],
                       [  1, 3]])
        b4by4 = array([[1, 0, 0, 0],
                       [0, 0, 0,1j],
                       [0, 1, 0, 0],
                       [0, 1, 0, 0]])
        for b in [b4, b4by1, b4by2, b4by4]:
            x = solve_banded((l, u), ab, b)
            assert_array_almost_equal(dot(a, x), b)

    def test_bad_shape(self):
        ab = array([[ 0.0, 20,  6,  2],
                    [   1,  4, 20, 14],
                    [ -30,  1,  7,  0],
                    [   2, -1,  0,  0]])
        l,u = 2,1
        bad = array([1.0, 2.0, 3.0, 4.0]).reshape(-1,4)
        assert_raises(ValueError, solve_banded, (l, u), ab, bad)
        assert_raises(ValueError, solve_banded, (l, u), ab, [1.0, 2.0])

        # Values of (l,u) are not compatible with ab.
        assert_raises(ValueError, solve_banded, (1, 1), ab, [1.0, 2.0])


class TestSolveHBanded(TestCase):
    # solveh_banded currently has a DeprecationWarning.  When the warning
    # is removed in scipy 0.9, the 'ignore' filters and the test for the
    # warning can be removed.

    def test_01_upper(self):
        warnings.simplefilter('ignore', category=DeprecationWarning)
        # Solve
        # [ 4 1 0]     [1]
        # [ 1 4 1] X = [4]
        # [ 0 1 4]     [1]
        # with the RHS as a 1D array.
        ab = array([[-99, 1.0, 1.0], [4.0, 4.0, 4.0]])
        b = array([1.0, 4.0, 1.0])
        c, x = solveh_banded(ab, b)
        assert_array_almost_equal(x, [0.0, 1.0, 0.0])
        # Remove the following part of this test in scipy 0.9.
        a = array([[4.0, 1.0, 0.0], [1.0, 4.0, 1.0], [0.0, 1.0, 4.0]])
        fac = zeros_like(a)
        fac[range(3),range(3)] = c[-1]
        fac[(0,1),(1,2)] = c[0,1:]
        assert_array_almost_equal(a, dot(fac.T, fac))

    def test_02_upper(self):
        warnings.simplefilter('ignore', category=DeprecationWarning)
        # Solve
        # [ 4 1 0]     [1 4]
        # [ 1 4 1] X = [4 2]
        # [ 0 1 4]     [1 4]
        #
        ab = array([[-99, 1.0, 1.0],
                    [4.0, 4.0, 4.0]])
        b = array([[1.0, 4.0],
                   [4.0, 2.0],
                   [1.0, 4.0]])
        c, x = solveh_banded(ab, b)
        expected = array([[0.0, 1.0],
                          [1.0, 0.0],
                          [0.0, 1.0]])
        assert_array_almost_equal(x, expected)

    def test_03_upper(self):
        warnings.simplefilter('ignore', category=DeprecationWarning)
        # Solve
        # [ 4 1 0]     [1]
        # [ 1 4 1] X = [4]
        # [ 0 1 4]     [1]
        # with the RHS as a 2D array with shape (3,1).
        ab = array([[-99, 1.0, 1.0], [4.0, 4.0, 4.0]])
        b = array([1.0, 4.0, 1.0]).reshape(-1,1)
        c, x = solveh_banded(ab, b)
        assert_array_almost_equal(x, array([0.0, 1.0, 0.0]).reshape(-1,1))

    def test_01_lower(self):
        warnings.simplefilter('ignore', category=DeprecationWarning)
        # Solve
        # [ 4 1 0]     [1]
        # [ 1 4 1] X = [4]
        # [ 0 1 4]     [1]
        #
        ab = array([[4.0, 4.0, 4.0],
                    [1.0, 1.0, -99]])
        b = array([1.0, 4.0, 1.0])
        c, x = solveh_banded(ab, b, lower=True)
        assert_array_almost_equal(x, [0.0, 1.0, 0.0])

    def test_02_lower(self):
        warnings.simplefilter('ignore', category=DeprecationWarning)
        # Solve
        # [ 4 1 0]     [1 4]
        # [ 1 4 1] X = [4 2]
        # [ 0 1 4]     [1 4]
        #
        ab = array([[4.0, 4.0, 4.0],
                    [1.0, 1.0, -99]])
        b = array([[1.0, 4.0],
                   [4.0, 2.0],
                   [1.0, 4.0]])
        c, x = solveh_banded(ab, b, lower=True)
        expected = array([[0.0, 1.0],
                          [1.0, 0.0],
                          [0.0, 1.0]])
        assert_array_almost_equal(x, expected)

    def test_01_float32(self):
        warnings.simplefilter('ignore', category=DeprecationWarning)
        # Solve
        # [ 4 1 0]     [1]
        # [ 1 4 1] X = [4]
        # [ 0 1 4]     [1]
        #
        ab = array([[-99, 1.0, 1.0], [4.0, 4.0, 4.0]], dtype=float32)
        b = array([1.0, 4.0, 1.0], dtype=float32)
        c, x = solveh_banded(ab, b)
        assert_array_almost_equal(x, [0.0, 1.0, 0.0])

    def test_02_float32(self):
        warnings.simplefilter('ignore', category=DeprecationWarning)
        # Solve
        # [ 4 1 0]     [1 4]
        # [ 1 4 1] X = [4 2]
        # [ 0 1 4]     [1 4]
        #
        ab = array([[-99, 1.0, 1.0],
                    [4.0, 4.0, 4.0]], dtype=float32)
        b = array([[1.0, 4.0],
                   [4.0, 2.0],
                   [1.0, 4.0]], dtype=float32)
        c, x = solveh_banded(ab, b)
        expected = array([[0.0, 1.0],
                          [1.0, 0.0],
                          [0.0, 1.0]])
        assert_array_almost_equal(x, expected)

    def test_01_complex(self):
        warnings.simplefilter('ignore', category=DeprecationWarning)
        # Solve
        # [ 4 -j 0]     [ -j]
        # [ j 4 -j] X = [4-j]
        # [ 0 j  4]     [4+j]
        #
        ab = array([[-99, -1.0j, -1.0j], [4.0, 4.0, 4.0]])
        b = array([-1.0j, 4.0-1j, 4+1j])
        c, x = solveh_banded(ab, b)
        assert_array_almost_equal(x, [0.0, 1.0, 1.0])

    def test_02_complex(self):
        warnings.simplefilter('ignore', category=DeprecationWarning)
        # Solve
        # [ 4 -j 0]     [ -j    4j]
        # [ j 4 -j] X = [4-j  -1-j]
        # [ 0 j  4]     [4+j   4  ]
        #
        ab = array([[-99, -1.0j, -1.0j],
                    [4.0, 4.0, 4.0]])
        b = array([[   -1j,    4.0j],
                   [4.0-1j, -1.0-1j],
                   [4.0+1j,     4.0]])
        c, x = solveh_banded(ab, b)
        expected = array([[0.0, 1.0j],
                          [1.0,  0.0],
                          [1.0,  1.0]])
        assert_array_almost_equal(x, expected)

    def test_bad_shapes(self):
        warnings.simplefilter('ignore', category=DeprecationWarning)

        ab = array([[-99, 1.0, 1.0],
                    [4.0, 4.0, 4.0]])
        b = array([[1.0, 4.0],
                   [4.0, 2.0]])
        assert_raises(ValueError, solveh_banded, ab, b)
        assert_raises(ValueError, solveh_banded, ab, [1.0, 2.0])
        assert_raises(ValueError, solveh_banded, ab, [1.0])

    def test_00_deprecation_warning(self):
        warnings.simplefilter('error', category=DeprecationWarning)
        ab = array([[-99, 1.0, 1.0], [4.0, 4.0, 4.0]])
        b = array([1.0, 4.0, 1.0])
        assert_raises(DeprecationWarning, solveh_banded, ab, b)


class TestSolve(TestCase):

    def test_20Feb04_bug(self):
        a = [[1,1],[1.0,0]] # ok
        x0 = solve(a,[1,0j])
        assert_array_almost_equal(dot(a,x0),[1,0])

        a = [[1,1],[1.2,0]] # gives failure with clapack.zgesv(..,rowmajor=0)
        b = [1,0j]
        x0 = solve(a,b)
        assert_array_almost_equal(dot(a,x0),[1,0])

    def test_simple(self):
        a = [[1,20],[-30,4]]
        for b in ([[1,0],[0,1]],[1,0],
                  [[2,1],[-30,4]]):
            x = solve(a,b)
            assert_array_almost_equal(dot(a,x),b)

    def test_simple_sym(self):
        a = [[2,3],[3,5]]
        for lower in [0,1]:
            for b in ([[1,0],[0,1]],[1,0]):
                x = solve(a,b,sym_pos=1,lower=lower)
                assert_array_almost_equal(dot(a,x),b)

    def test_simple_sym_complex(self):
        a = [[5,2],[2,4]]
        for b in [[1j,0],
                  [[1j,1j],
                   [0,2]],
                  ]:
            x = solve(a,b,sym_pos=1)
            assert_array_almost_equal(dot(a,x),b)

    def test_simple_complex(self):
        a = array([[5,2],[2j,4]],'D')
        for b in [[1j,0],
                  [[1j,1j],
                   [0,2]],
                  [1,0j],
                  array([1,0],'D'),
                  ]:
            x = solve(a,b)
            assert_array_almost_equal(dot(a,x),b)

    def test_nils_20Feb04(self):
        n = 2
        A = random([n,n])+random([n,n])*1j
        X = zeros((n,n),'D')
        Ainv = inv(A)
        R = identity(n)+identity(n)*0j
        for i in arange(0,n):
            r = R[:,i]
            X[:,i] = solve(A,r)
        assert_array_almost_equal(X,Ainv)

    def test_random(self):

        n = 20
        a = random([n,n])
        for i in range(n): a[i,i] = 20*(.1+a[i,i])
        for i in range(4):
            b = random([n,3])
            x = solve(a,b)
            assert_array_almost_equal(dot(a,x),b)

    def test_random_complex(self):
        n = 20
        a = random([n,n]) + 1j * random([n,n])
        for i in range(n): a[i,i] = 20*(.1+a[i,i])
        for i in range(2):
            b = random([n,3])
            x = solve(a,b)
            assert_array_almost_equal(dot(a,x),b)

    def test_random_sym(self):
        n = 20
        a = random([n,n])
        for i in range(n):
            a[i,i] = abs(20*(.1+a[i,i]))
            for j in range(i):
                a[i,j] = a[j,i]
        for i in range(4):
            b = random([n])
            x = solve(a,b,sym_pos=1)
            assert_array_almost_equal(dot(a,x),b)

    def test_random_sym_complex(self):
        n = 20
        a = random([n,n])
        #a  = a + 1j*random([n,n]) # XXX: with this the accuracy will be very low
        for i in range(n):
            a[i,i] = abs(20*(.1+a[i,i]))
            for j in range(i):
                a[i,j] = conjugate(a[j,i])
        b = random([n])+2j*random([n])
        for i in range(2):
            x = solve(a,b,sym_pos=1)
            assert_array_almost_equal(dot(a,x),b)


class TestInv(TestCase):

    def test_simple(self):
        a = [[1,2],[3,4]]
        a_inv = inv(a)
        assert_array_almost_equal(dot(a,a_inv),
                                  [[1,0],[0,1]])
        a = [[1,2,3],[4,5,6],[7,8,10]]
        a_inv = inv(a)
        assert_array_almost_equal(dot(a,a_inv),
                                  [[1,0,0],[0,1,0],[0,0,1]])

    def test_random(self):
        n = 20
        for i in range(4):
            a = random([n,n])
            for i in range(n): a[i,i] = 20*(.1+a[i,i])
            a_inv = inv(a)
            assert_array_almost_equal(dot(a,a_inv),
                                      identity(n))
    def test_simple_complex(self):
        a = [[1,2],[3,4j]]
        a_inv = inv(a)
        assert_array_almost_equal(dot(a,a_inv),
                                  [[1,0],[0,1]])

    def test_random_complex(self):
        n = 20
        for i in range(4):
            a = random([n,n])+2j*random([n,n])
            for i in range(n): a[i,i] = 20*(.1+a[i,i])
            a_inv = inv(a)
            assert_array_almost_equal(dot(a,a_inv),
                                      identity(n))


class TestDet(TestCase):

    def test_simple(self):
        a = [[1,2],[3,4]]
        a_det = det(a)
        assert_almost_equal(a_det,-2.0)

    def test_simple_complex(self):
        a = [[1,2],[3,4j]]
        a_det = det(a)
        assert_almost_equal(a_det,-6+4j)

    def test_random(self):
        basic_det = linalg.det
        n = 20
        for i in range(4):
            a = random([n,n])
            d1 = det(a)
            d2 = basic_det(a)
            assert_almost_equal(d1,d2)

    def test_random_complex(self):
        basic_det = linalg.det
        n = 20
        for i in range(4):
            a = random([n,n]) + 2j*random([n,n])
            d1 = det(a)
            d2 = basic_det(a)
            assert_almost_equal(d1,d2)


def direct_lstsq(a,b,cmplx=0):
    at = transpose(a)
    if cmplx:
        at = conjugate(at)
    a1 = dot(at, a)
    b1 = dot(at, b)
    return solve(a1, b1)

class TestLstsq(TestCase):
    def test_random_overdet_large(self):
        #bug report: Nils Wagner
        n = 200
        a = random([n,2])
        for i in range(2): a[i,i] = 20*(.1+a[i,i])
        b = random([n,3])
        x = lstsq(a,b)[0]
        assert_array_almost_equal(x,direct_lstsq(a,b))

    def test_simple_exact(self):
        a = [[1,20],[-30,4]]
        for b in ([[1,0],[0,1]],[1,0],
                  [[2,1],[-30,4]]):
            x = lstsq(a,b)[0]
            assert_array_almost_equal(dot(a,x),b)

    def test_simple_overdet(self):
        a = [[1,2],[4,5],[3,4]]
        b = [1,2,3]
        x,res,r,s = lstsq(a,b)
        #XXX: check defintion of res
        assert_array_almost_equal(x,direct_lstsq(a,b))

    def test_simple_underdet(self):
        a = [[1,2,3],[4,5,6]]
        b = [1,2]
        x,res,r,s = lstsq(a,b)
        #XXX: need independent check
        assert_array_almost_equal(x,[[-0.05555556],
                                     [0.11111111],[0.27777778]])

    def test_random_exact(self):

        n = 20
        a = random([n,n])
        for i in range(n): a[i,i] = 20*(.1+a[i,i])
        for i in range(4):
            b = random([n,3])
            x = lstsq(a,b)[0]
            assert_array_almost_equal(dot(a,x),b)

    def test_random_complex_exact(self):
        n = 20
        a = random([n,n]) + 1j * random([n,n])
        for i in range(n): a[i,i] = 20*(.1+a[i,i])
        for i in range(2):
            b = random([n,3])
            x = lstsq(a,b)[0]
            assert_array_almost_equal(dot(a,x),b)

    def test_random_overdet(self):
        n = 20
        m = 15
        a = random([n,m])
        for i in range(m): a[i,i] = 20*(.1+a[i,i])
        for i in range(4):
            b = random([n,3])
            x,res,r,s = lstsq(a,b)
            assert r==m,'unexpected efficient rank'
            #XXX: check definition of res
            assert_array_almost_equal(x,direct_lstsq(a,b))

    def test_random_complex_overdet(self):
        n = 20
        m = 15
        a = random([n,m]) + 1j * random([n,m])
        for i in range(m):
            a[i,i] = 20*(.1+a[i,i])
        for i in range(2):
            b = random([n,3])
            x,res,r,s = lstsq(a,b)
            assert r==m,'unexpected efficient rank'
            #XXX: check definition of res
            assert_array_almost_equal(x,direct_lstsq(a,b,1))


class TestPinv(TestCase):

    def test_simple(self):
        a=array([[1,2,3],[4,5,6.],[7,8,10]])
        a_pinv = pinv(a)
        assert_array_almost_equal(dot(a,a_pinv),[[1,0,0],[0,1,0],[0,0,1]])
        a_pinv = pinv2(a)
        assert_array_almost_equal(dot(a,a_pinv),[[1,0,0],[0,1,0],[0,0,1]])

    def test_simple_0det(self):
        a=array([[1,2,3],[4,5,6.],[7,8,9]])
        a_pinv = pinv(a)
        a_pinv2 = pinv2(a)
        assert_array_almost_equal(a_pinv,a_pinv2)

    def test_simple_cols(self):
        a=array([[1,2,3],[4,5,6.]])
        a_pinv = pinv(a)
        a_pinv2 = pinv2(a)
        assert_array_almost_equal(a_pinv,a_pinv2)

    def test_simple_rows(self):
        a=array([[1,2],[3,4],[5,6]])
        a_pinv = pinv(a)
        a_pinv2 = pinv2(a)
        assert_array_almost_equal(a_pinv,a_pinv2)

class TestNorm(object):
    def test_zero_norm(self):
        assert_equal(norm([1,0,3], 0), 2)
        assert_equal(norm([1,2,3], 0), 3)

if __name__ == "__main__":
    run_module_suite()
