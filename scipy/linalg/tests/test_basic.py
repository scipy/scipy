#!/usr/bin/env python
#
# Created by: Pearu Peterson, March 2002
#
""" Test functions for linalg.basic module

"""
from __future__ import division, print_function, absolute_import

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

import numpy as np
from numpy import (arange, array, dot, zeros, identity, conjugate, transpose,
        float32)
import numpy.linalg as linalg

from numpy.testing import (TestCase, rand, run_module_suite, assert_raises,
    assert_equal, assert_almost_equal, assert_array_almost_equal, assert_,
    assert_allclose, assert_array_equal)

from scipy.linalg import (solve, inv, det, lstsq, pinv, pinv2, pinvh, norm,
        solve_banded, solveh_banded, solve_triangular, solve_circulant,
        circulant, LinAlgError)

from scipy.linalg._testutils import assert_no_overwrite

REAL_DTYPES = [np.float32, np.float64]
COMPLEX_DTYPES = [np.complex64, np.complex128]
DTYPES = REAL_DTYPES + COMPLEX_DTYPES

def random(size):
    return rand(*size)


class TestSolveBanded(TestCase):

    def test_real(self):
        a = array([[1.0, 20, 0, 0],
                   [-30, 4, 6, 0],
                   [2, 1, 20, 2],
                   [0, -1, 7, 14]])
        ab = array([[0.0, 20, 6, 2],
                    [1, 4, 20, 14],
                    [-30, 1, 7, 0],
                    [2, -1, 0, 0]])
        l,u = 2,1
        b4 = array([10.0, 0.0, 2.0, 14.0])
        b4by1 = b4.reshape(-1,1)
        b4by2 = array([[2, 1],
                       [-30, 4],
                       [2, 3],
                       [1, 3]])
        b4by4 = array([[1, 0, 0, 0],
                       [0, 0, 0, 1],
                       [0, 1, 0, 0],
                       [0, 1, 0, 0]])
        for b in [b4, b4by1, b4by2, b4by4]:
            x = solve_banded((l, u), ab, b)
            assert_array_almost_equal(dot(a, x), b)

    def test_complex(self):
        a = array([[1.0, 20, 0, 0],
                   [-30, 4, 6, 0],
                   [2j, 1, 20, 2j],
                   [0, -1, 7, 14]])
        ab = array([[0.0, 20, 6, 2j],
                    [1, 4, 20, 14],
                    [-30, 1, 7, 0],
                    [2j, -1, 0, 0]])
        l,u = 2,1
        b4 = array([10.0, 0.0, 2.0, 14.0j])
        b4by1 = b4.reshape(-1,1)
        b4by2 = array([[2, 1],
                       [-30, 4],
                       [2, 3],
                       [1, 3]])
        b4by4 = array([[1, 0, 0, 0],
                       [0, 0, 0,1j],
                       [0, 1, 0, 0],
                       [0, 1, 0, 0]])
        for b in [b4, b4by1, b4by2, b4by4]:
            x = solve_banded((l, u), ab, b)
            assert_array_almost_equal(dot(a, x), b)

    def test_tridiag_real(self):
        ab = array([[0.0, 20, 6, 2],
                   [1, 4, 20, 14],
                   [-30, 1, 7, 0]])
        a = np.diag(ab[0,1:], 1) + np.diag(ab[1,:], 0) + np.diag(ab[2,:-1], -1)
        b4 = array([10.0, 0.0, 2.0, 14.0])
        b4by1 = b4.reshape(-1,1)
        b4by2 = array([[2, 1],
                       [-30, 4],
                       [2, 3],
                       [1, 3]])
        b4by4 = array([[1, 0, 0, 0],
                       [0, 0, 0, 1],
                       [0, 1, 0, 0],
                       [0, 1, 0, 0]])
        for b in [b4, b4by1, b4by2, b4by4]:
            x = solve_banded((1, 1), ab, b)
            assert_array_almost_equal(dot(a, x), b)

    def test_tridiag_complex(self):
        ab = array([[0.0, 20, 6, 2j],
                   [1, 4, 20, 14],
                   [-30, 1, 7, 0]])
        a = np.diag(ab[0,1:], 1) + np.diag(ab[1,:], 0) + np.diag(ab[2,:-1], -1)
        b4 = array([10.0, 0.0, 2.0, 14.0j])
        b4by1 = b4.reshape(-1,1)
        b4by2 = array([[2, 1],
                       [-30, 4],
                       [2, 3],
                       [1, 3]])
        b4by4 = array([[1, 0, 0, 0],
                       [0, 0, 0, 1],
                       [0, 1, 0, 0],
                       [0, 1, 0, 0]])
        for b in [b4, b4by1, b4by2, b4by4]:
            x = solve_banded((1, 1), ab, b)
            assert_array_almost_equal(dot(a, x), b)

    def test_check_finite(self):
        a = array([[1.0, 20, 0, 0],
                   [-30, 4, 6, 0],
                   [2, 1, 20, 2],
                   [0, -1, 7, 14]])
        ab = array([[0.0, 20, 6, 2],
                    [1, 4, 20, 14],
                    [-30, 1, 7, 0],
                    [2, -1, 0, 0]])
        l,u = 2,1
        b4 = array([10.0, 0.0, 2.0, 14.0])
        x = solve_banded((l, u), ab, b4, check_finite=False)
        assert_array_almost_equal(dot(a, x), b4)

    def test_bad_shape(self):
        ab = array([[0.0, 20, 6, 2],
                    [1, 4, 20, 14],
                    [-30, 1, 7, 0],
                    [2, -1, 0, 0]])
        l,u = 2,1
        bad = array([1.0, 2.0, 3.0, 4.0]).reshape(-1,4)
        assert_raises(ValueError, solve_banded, (l, u), ab, bad)
        assert_raises(ValueError, solve_banded, (l, u), ab, [1.0, 2.0])

        # Values of (l,u) are not compatible with ab.
        assert_raises(ValueError, solve_banded, (1, 1), ab, [1.0, 2.0])

    def test_1x1(self):
        x = solve_banded((1, 1), [[0], [1], [0]], [[1, 2, 3]])
        assert_array_equal(x, [[1.0, 2.0, 3.0]])
        assert_equal(x.dtype, np.dtype('f8'))


class TestSolveHBanded(TestCase):

    def test_01_upper(self):
        # Solve
        # [ 4 1 2 0]     [1]
        # [ 1 4 1 2] X = [4]
        # [ 2 1 4 1]     [1]
        # [ 0 2 1 4]     [2]
        # with the RHS as a 1D array.
        ab = array([[0.0, 0.0, 2.0, 2.0],
                    [-99, 1.0, 1.0, 1.0],
                    [4.0, 4.0, 4.0, 4.0]])
        b = array([1.0, 4.0, 1.0, 2.0])
        x = solveh_banded(ab, b)
        assert_array_almost_equal(x, [0.0, 1.0, 0.0, 0.0])

    def test_02_upper(self):
        # Solve
        # [ 4 1 2 0]     [1 6]
        # [ 1 4 1 2] X = [4 2]
        # [ 2 1 4 1]     [1 6]
        # [ 0 2 1 4]     [2 1]
        #
        ab = array([[0.0, 0.0, 2.0, 2.0],
                    [-99, 1.0, 1.0, 1.0],
                    [4.0, 4.0, 4.0, 4.0]])
        b = array([[1.0, 6.0],
                   [4.0, 2.0],
                   [1.0, 6.0],
                   [2.0, 1.0]])
        x = solveh_banded(ab, b)
        expected = array([[0.0, 1.0],
                          [1.0, 0.0],
                          [0.0, 1.0],
                          [0.0, 0.0]])
        assert_array_almost_equal(x, expected)

    def test_03_upper(self):
        # Solve
        # [ 4 1 2 0]     [1]
        # [ 1 4 1 2] X = [4]
        # [ 2 1 4 1]     [1]
        # [ 0 2 1 4]     [2]
        # with the RHS as a 2D array with shape (3,1).
        ab = array([[0.0, 0.0, 2.0, 2.0],
                    [-99, 1.0, 1.0, 1.0],
                    [4.0, 4.0, 4.0, 4.0]])
        b = array([1.0, 4.0, 1.0, 2.0]).reshape(-1,1)
        x = solveh_banded(ab, b)
        assert_array_almost_equal(x, array([0.0, 1.0, 0.0, 0.0]).reshape(-1,1))

    def test_01_lower(self):
        # Solve
        # [ 4 1 2 0]     [1]
        # [ 1 4 1 2] X = [4]
        # [ 2 1 4 1]     [1]
        # [ 0 2 1 4]     [2]
        #
        ab = array([[4.0, 4.0, 4.0, 4.0],
                    [1.0, 1.0, 1.0, -99],
                    [2.0, 2.0, 0.0, 0.0]])
        b = array([1.0, 4.0, 1.0, 2.0])
        x = solveh_banded(ab, b, lower=True)
        assert_array_almost_equal(x, [0.0, 1.0, 0.0, 0.0])

    def test_02_lower(self):
        # Solve
        # [ 4 1 2 0]     [1 6]
        # [ 1 4 1 2] X = [4 2]
        # [ 2 1 4 1]     [1 6]
        # [ 0 2 1 4]     [2 1]
        #
        ab = array([[4.0, 4.0, 4.0, 4.0],
                    [1.0, 1.0, 1.0, -99],
                    [2.0, 2.0, 0.0, 0.0]])
        b = array([[1.0, 6.0],
                   [4.0, 2.0],
                   [1.0, 6.0],
                   [2.0, 1.0]])
        x = solveh_banded(ab, b, lower=True)
        expected = array([[0.0, 1.0],
                          [1.0, 0.0],
                          [0.0, 1.0],
                          [0.0, 0.0]])
        assert_array_almost_equal(x, expected)

    def test_01_float32(self):
        # Solve
        # [ 4 1 2 0]     [1]
        # [ 1 4 1 2] X = [4]
        # [ 2 1 4 1]     [1]
        # [ 0 2 1 4]     [2]
        #
        ab = array([[0.0, 0.0, 2.0, 2.0],
                    [-99, 1.0, 1.0, 1.0],
                    [4.0, 4.0, 4.0, 4.0]], dtype=float32)
        b = array([1.0, 4.0, 1.0, 2.0], dtype=float32)
        x = solveh_banded(ab, b)
        assert_array_almost_equal(x, [0.0, 1.0, 0.0, 0.0])

    def test_02_float32(self):
        # Solve
        # [ 4 1 2 0]     [1 6]
        # [ 1 4 1 2] X = [4 2]
        # [ 2 1 4 1]     [1 6]
        # [ 0 2 1 4]     [2 1]
        #
        ab = array([[0.0, 0.0, 2.0, 2.0],
                    [-99, 1.0, 1.0, 1.0],
                    [4.0, 4.0, 4.0, 4.0]], dtype=float32)
        b = array([[1.0, 6.0],
                   [4.0, 2.0],
                   [1.0, 6.0],
                   [2.0, 1.0]], dtype=float32)
        x = solveh_banded(ab, b)
        expected = array([[0.0, 1.0],
                          [1.0, 0.0],
                          [0.0, 1.0],
                          [0.0, 0.0]])
        assert_array_almost_equal(x, expected)

    def test_01_complex(self):
        # Solve
        # [ 4 -j  2  0]     [2-j]
        # [ j  4 -j  2] X = [4-j]
        # [ 2  j  4 -j]     [4+j]
        # [ 0  2  j  4]     [2+j]
        #
        ab = array([[0.0, 0.0, 2.0, 2.0],
                    [-99, -1.0j, -1.0j, -1.0j],
                    [4.0, 4.0, 4.0, 4.0]])
        b = array([2-1.0j, 4.0-1j, 4+1j, 2+1j])
        x = solveh_banded(ab, b)
        assert_array_almost_equal(x, [0.0, 1.0, 1.0, 0.0])

    def test_02_complex(self):
        # Solve
        # [ 4 -j  2  0]     [2-j 2+4j]
        # [ j  4 -j  2] X = [4-j -1-j]
        # [ 2  j  4 -j]     [4+j 4+2j]
        # [ 0  2  j  4]     [2+j j]
        #
        ab = array([[0.0, 0.0, 2.0, 2.0],
                    [-99, -1.0j, -1.0j, -1.0j],
                    [4.0, 4.0, 4.0, 4.0]])
        b = array([[2-1j, 2+4j],
                   [4.0-1j, -1-1j],
                   [4.0+1j, 4+2j],
                   [2+1j, 1j]])
        x = solveh_banded(ab, b)
        expected = array([[0.0, 1.0j],
                          [1.0, 0.0],
                          [1.0, 1.0],
                          [0.0, 0.0]])
        assert_array_almost_equal(x, expected)

    def test_tridiag_01_upper(self):
        # Solve
        # [ 4 1 0]     [1]
        # [ 1 4 1] X = [4]
        # [ 0 1 4]     [1]
        # with the RHS as a 1D array.
        ab = array([[-99, 1.0, 1.0], [4.0, 4.0, 4.0]])
        b = array([1.0, 4.0, 1.0])
        x = solveh_banded(ab, b)
        assert_array_almost_equal(x, [0.0, 1.0, 0.0])

    def test_tridiag_02_upper(self):
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
        x = solveh_banded(ab, b)
        expected = array([[0.0, 1.0],
                          [1.0, 0.0],
                          [0.0, 1.0]])
        assert_array_almost_equal(x, expected)

    def test_tridiag_03_upper(self):
        # Solve
        # [ 4 1 0]     [1]
        # [ 1 4 1] X = [4]
        # [ 0 1 4]     [1]
        # with the RHS as a 2D array with shape (3,1).
        ab = array([[-99, 1.0, 1.0], [4.0, 4.0, 4.0]])
        b = array([1.0, 4.0, 1.0]).reshape(-1,1)
        x = solveh_banded(ab, b)
        assert_array_almost_equal(x, array([0.0, 1.0, 0.0]).reshape(-1,1))

    def test_tridiag_01_lower(self):
        # Solve
        # [ 4 1 0]     [1]
        # [ 1 4 1] X = [4]
        # [ 0 1 4]     [1]
        #
        ab = array([[4.0, 4.0, 4.0],
                    [1.0, 1.0, -99]])
        b = array([1.0, 4.0, 1.0])
        x = solveh_banded(ab, b, lower=True)
        assert_array_almost_equal(x, [0.0, 1.0, 0.0])

    def test_tridiag_02_lower(self):
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
        x = solveh_banded(ab, b, lower=True)
        expected = array([[0.0, 1.0],
                          [1.0, 0.0],
                          [0.0, 1.0]])
        assert_array_almost_equal(x, expected)

    def test_tridiag_01_float32(self):
        # Solve
        # [ 4 1 0]     [1]
        # [ 1 4 1] X = [4]
        # [ 0 1 4]     [1]
        #
        ab = array([[-99, 1.0, 1.0], [4.0, 4.0, 4.0]], dtype=float32)
        b = array([1.0, 4.0, 1.0], dtype=float32)
        x = solveh_banded(ab, b)
        assert_array_almost_equal(x, [0.0, 1.0, 0.0])

    def test_tridiag_02_float32(self):
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
        x = solveh_banded(ab, b)
        expected = array([[0.0, 1.0],
                          [1.0, 0.0],
                          [0.0, 1.0]])
        assert_array_almost_equal(x, expected)

    def test_tridiag_01_complex(self):
        # Solve
        # [ 4 -j 0]     [ -j]
        # [ j 4 -j] X = [4-j]
        # [ 0 j  4]     [4+j]
        #
        ab = array([[-99, -1.0j, -1.0j], [4.0, 4.0, 4.0]])
        b = array([-1.0j, 4.0-1j, 4+1j])
        x = solveh_banded(ab, b)
        assert_array_almost_equal(x, [0.0, 1.0, 1.0])

    def test_tridiag_02_complex(self):
        # Solve
        # [ 4 -j 0]     [ -j    4j]
        # [ j 4 -j] X = [4-j  -1-j]
        # [ 0 j  4]     [4+j   4  ]
        #
        ab = array([[-99, -1.0j, -1.0j],
                    [4.0, 4.0, 4.0]])
        b = array([[-1j, 4.0j],
                   [4.0-1j, -1.0-1j],
                   [4.0+1j, 4.0]])
        x = solveh_banded(ab, b)
        expected = array([[0.0, 1.0j],
                          [1.0, 0.0],
                          [1.0, 1.0]])
        assert_array_almost_equal(x, expected)

    def test_check_finite(self):
        # Solve
        # [ 4 1 0]     [1]
        # [ 1 4 1] X = [4]
        # [ 0 1 4]     [1]
        # with the RHS as a 1D array.
        ab = array([[-99, 1.0, 1.0], [4.0, 4.0, 4.0]])
        b = array([1.0, 4.0, 1.0])
        x = solveh_banded(ab, b, check_finite=False)
        assert_array_almost_equal(x, [0.0, 1.0, 0.0])

    def test_bad_shapes(self):
        ab = array([[-99, 1.0, 1.0],
                    [4.0, 4.0, 4.0]])
        b = array([[1.0, 4.0],
                   [4.0, 2.0]])
        assert_raises(ValueError, solveh_banded, ab, b)
        assert_raises(ValueError, solveh_banded, ab, [1.0, 2.0])
        assert_raises(ValueError, solveh_banded, ab, [1.0])

    def test_1x1(self):
        x = solveh_banded([[1]], [[1, 2, 3]])
        assert_array_equal(x, [[1.0, 2.0, 3.0]])
        assert_equal(x.dtype, np.dtype('f8'))


class TestSolve(TestCase):
    def setUp(self):
        np.random.seed(1234)

    def test_20Feb04_bug(self):
        a = [[1,1],[1.0,0]]  # ok
        x0 = solve(a,[1,0j])
        assert_array_almost_equal(dot(a,x0),[1,0])

        a = [[1,1],[1.2,0]]  # gives failure with clapack.zgesv(..,rowmajor=0)
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
        for i in range(n):
            a[i,i] = 20*(.1+a[i,i])
        for i in range(4):
            b = random([n,3])
            x = solve(a,b)
            assert_array_almost_equal(dot(a,x),b)

    def test_random_complex(self):
        n = 20
        a = random([n,n]) + 1j * random([n,n])
        for i in range(n):
            a[i,i] = 20*(.1+a[i,i])
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
        # a  = a + 1j*random([n,n]) # XXX: with this the accuracy will be very low
        for i in range(n):
            a[i,i] = abs(20*(.1+a[i,i]))
            for j in range(i):
                a[i,j] = conjugate(a[j,i])
        b = random([n])+2j*random([n])
        for i in range(2):
            x = solve(a,b,sym_pos=1)
            assert_array_almost_equal(dot(a,x),b)

    def test_check_finite(self):
        a = [[1,20],[-30,4]]
        for b in ([[1,0],[0,1]],[1,0],
                  [[2,1],[-30,4]]):
            x = solve(a,b, check_finite=False)
            assert_array_almost_equal(dot(a,x),b)


class TestSolveTriangular(TestCase):

    def test_simple(self):
        """
        solve_triangular on a simple 2x2 matrix.
        """
        A = array([[1,0], [1,2]])
        b = [1, 1]
        sol = solve_triangular(A, b, lower=True)
        assert_array_almost_equal(sol, [1, 0])

        # check that it works also for non-contiguous matrices
        sol = solve_triangular(A.T, b, lower=False)
        assert_array_almost_equal(sol, [.5, .5])

        # and that it gives the same result as trans=1
        sol = solve_triangular(A, b, lower=True, trans=1)
        assert_array_almost_equal(sol, [.5, .5])

        b = identity(2)
        sol = solve_triangular(A, b, lower=True, trans=1)
        assert_array_almost_equal(sol, [[1., -.5], [0, 0.5]])

    def test_simple_complex(self):
        """
        solve_triangular on a simple 2x2 complex matrix
        """
        A = array([[1+1j, 0], [1j, 2]])
        b = identity(2)
        sol = solve_triangular(A, b, lower=True, trans=1)
        assert_array_almost_equal(sol, [[.5-.5j, -.25-.25j], [0, 0.5]])

    def test_check_finite(self):
        """
        solve_triangular on a simple 2x2 matrix.
        """
        A = array([[1,0], [1,2]])
        b = [1, 1]
        sol = solve_triangular(A, b, lower=True, check_finite=False)
        assert_array_almost_equal(sol, [1, 0])


class TestInv(TestCase):
    def setUp(self):
        np.random.seed(1234)

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
            for i in range(n):
                a[i,i] = 20*(.1+a[i,i])
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
            for i in range(n):
                a[i,i] = 20*(.1+a[i,i])
            a_inv = inv(a)
            assert_array_almost_equal(dot(a,a_inv),
                                      identity(n))

    def test_check_finite(self):
        a = [[1,2],[3,4]]
        a_inv = inv(a, check_finite=False)
        assert_array_almost_equal(dot(a,a_inv),
                                  [[1,0],[0,1]])


class TestDet(TestCase):
    def setUp(self):
        np.random.seed(1234)

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
            assert_allclose(d1, d2, rtol=1e-13)

    def test_check_finite(self):
        a = [[1,2],[3,4]]
        a_det = det(a, check_finite=False)
        assert_almost_equal(a_det,-2.0)


def direct_lstsq(a,b,cmplx=0):
    at = transpose(a)
    if cmplx:
        at = conjugate(at)
    a1 = dot(at, a)
    b1 = dot(at, b)
    return solve(a1, b1)


class TestLstsq(TestCase):

    lapack_drivers = ('gelsd', 'gelss', 'gelsd')

    def setUp(self):
        np.random.seed(1234)

    def test_simple_exact(self):
        for dtype in REAL_DTYPES:
            a = np.array(((1,20),(-30,4)), dtype=dtype)
            for lapack_driver in TestLstsq.lapack_drivers:
                    for overwrite in (True,False):
                        for bt in (((1,0),(0,1)), (1,0),((2,1),(-30,4))):
                            # Store values in case they are overwritten
                            # later
                            a1 = a.copy()
                            b = np.array(bt, dtype=dtype)
                            b1 = b.copy()
                            out = lstsq(a1,b1,lapack_driver=lapack_driver,
                                       overwrite_a=overwrite,
                                       overwrite_b=overwrite)
                            x = out[0]
                            r = out[2]
                            assert_(r == 2, 'unexpected efficient rank')
                            assert_allclose(dot(a,x),b,
                                            atol=25*np.finfo(a1.dtype).eps,
                                            rtol=25*np.finfo(a1.dtype).eps,
                                        err_msg="driver: %s" % lapack_driver)

    def test_simple_overdet(self):

        for dtype in REAL_DTYPES:
            a = np.array(((1,2),(4,5),(3,4)), dtype=dtype)
            b = np.array((1,2,3), dtype=dtype)
            for lapack_driver in TestLstsq.lapack_drivers:
                for overwrite in (True,False):
                    # Store values in case they are overwritten
                    # later
                    a1 = a.copy()
                    b1 = b.copy()
                    out = lstsq(a1,b1,lapack_driver=lapack_driver,
                                overwrite_a=overwrite, overwrite_b=overwrite)
                    x = out[0]
                    residuals = out[1]
                    r = out[2]
                    assert_(r == 2, 'unexpected efficient rank')
                    assert_allclose(abs((dot(a,x) - b)**2).sum(axis=0),
                                    residuals,
                                    rtol=25*np.finfo(a1.dtype).eps,
                                    atol=25*np.finfo(a1.dtype).eps,
                                    err_msg="driver: %s" % lapack_driver)
                    assert_allclose(x,(-0.428571428571429, 0.85714285714285),
                                    rtol=25*np.finfo(a1.dtype).eps,
                                    atol=25*np.finfo(a1.dtype).eps,
                                    err_msg="driver: %s" % lapack_driver)

    def test_simple_overdet_complex(self):

        for dtype in COMPLEX_DTYPES:
            a = np.array(((1+2j,2),(4,5),(3,4)), dtype=dtype)
            b = np.array((1,2+4j,3), dtype=dtype)
            for lapack_driver in TestLstsq.lapack_drivers:
                for overwrite in (True,False):
                    # Store values in case they are overwritten
                    # later
                    a1 = a.copy()
                    b1 = b.copy()
                    out = lstsq(a1,b1,lapack_driver=lapack_driver,
                                overwrite_a=overwrite, overwrite_b=overwrite)
                    x = out[0]
                    residuals = out[1]
                    r = out[2]
                    assert_(r == 2, 'unexpected efficient rank')
                    assert_allclose(abs((dot(a,x) - b)**2).sum(axis=0),
                                        residuals,
                                        rtol=25*np.finfo(a1.dtype).eps,
                                        atol=25*np.finfo(a1.dtype).eps,
                                        err_msg="driver: %s" % lapack_driver)
                    assert_allclose(x,(-0.4831460674157303+0.258426966292135j,
                                        0.921348314606741+0.292134831460674j),
                                        rtol=25*np.finfo(a1.dtype).eps,
                                        atol=25*np.finfo(a1.dtype).eps,
                                        err_msg="driver: %s" % lapack_driver)

    def test_simple_underdet(self):
        for dtype in REAL_DTYPES:
            a = np.array(((1,2,3),(4,5,6)), dtype=dtype)
            b = np.array((1,2), dtype=dtype)
            for lapack_driver in TestLstsq.lapack_drivers:
                for overwrite in (True,False):
                    # Store values in case they are overwritten
                    # later
                    a1 = a.copy()
                    b1 = b.copy()
                    out = lstsq(a1,b1,lapack_driver=lapack_driver,
                                overwrite_a=overwrite, overwrite_b=overwrite)
                    x = out[0]
                    r = out[2]
                    assert_(r == 2, 'unexpected efficient rank')
                    assert_allclose(x,(-0.055555555555555, 0.111111111111111,
                            0.277777777777777),
                            rtol=25*np.finfo(a1.dtype).eps,
                            atol=25*np.finfo(a1.dtype).eps,
                            err_msg="driver: %s" % lapack_driver)

    def test_random_exact(self):
        for dtype in REAL_DTYPES:
            for n in (20, 200):
                for lapack_driver in TestLstsq.lapack_drivers:
                    for overwrite in (True,False):
                        a = np.asarray(random([n,n]), dtype=dtype)
                        for i in range(n):
                            a[i,i] = 20*(.1+a[i,i])
                        for i in range(4):
                            b = np.asarray(random([n,3]), dtype=dtype)
                            # Store values in case they are overwritten
                            # later
                            a1 = a.copy()
                            b1 = b.copy()
                            out = lstsq(a1,b1,lapack_driver=lapack_driver,
                                       overwrite_a=overwrite,
                                       overwrite_b=overwrite)
                            x = out[0]
                            r = out[2]
                            assert_(r == n, 'unexpected efficient rank')
                            if dtype is np.float32:
                                assert_allclose(dot(a,x),b,
                                          rtol=400*np.finfo(a1.dtype).eps,
                                          atol=400*np.finfo(a1.dtype).eps,
                                          err_msg="driver: %s" % lapack_driver)
                            else:
                                assert_allclose(dot(a,x),b,
                                          rtol=1000*np.finfo(a1.dtype).eps,
                                          atol=1000*np.finfo(a1.dtype).eps,
                                          err_msg="driver: %s" % lapack_driver)

    def test_random_complex_exact(self):
        for dtype in COMPLEX_DTYPES:
            for n in (20, 200):
                for lapack_driver in TestLstsq.lapack_drivers:
                    for overwrite in (True,False):
                        a = np.asarray(random([n,n]) + 1j*random([n,n]),
                                       dtype=dtype)
                        for i in range(n):
                            a[i,i] = 20*(.1+a[i,i])
                        for i in range(2):
                            b = np.asarray(random([n,3]), dtype=dtype)
                            # Store values in case they are overwritten
                            # later
                            a1 = a.copy()
                            b1 = b.copy()
                            out = lstsq(a1,b1,lapack_driver=lapack_driver,
                                       overwrite_a=overwrite,
                                       overwrite_b=overwrite)
                            x = out[0]
                            r = out[2]
                            assert_(r == n, 'unexpected efficient rank')
                            if dtype is np.complex64:
                                assert_allclose(dot(a,x),b,
                                          rtol=400*np.finfo(a1.dtype).eps,
                                          atol=400*np.finfo(a1.dtype).eps,
                                          err_msg="driver: %s" % lapack_driver)
                            else:
                                assert_allclose(dot(a,x),b,
                                          rtol=1000*np.finfo(a1.dtype).eps,
                                          atol=1000*np.finfo(a1.dtype).eps,
                                          err_msg="driver: %s" % lapack_driver)

    def test_random_overdet(self):

        for dtype in REAL_DTYPES:
            for (n,m) in ((20,15), (200,2)):
                for lapack_driver in TestLstsq.lapack_drivers:
                    for overwrite in (True,False):
                        a = np.asarray(random([n,m]),dtype=dtype)
                        print(a.dtype)
                        for i in range(m):
                            a[i,i] = 20*(.1+a[i,i])
                        for i in range(4):
                            b = np.asarray(random([n,3]),dtype=dtype)
                            # Store values in case they are overwritten
                            # later
                            a1 = a.copy()
                            b1 = b.copy()
                            out = lstsq(a1,b1,lapack_driver=lapack_driver,
                                       overwrite_a=overwrite,
                                       overwrite_b=overwrite)
                            x = out[0]
                            r = out[2]
                            assert_(r == m, 'unexpected efficient rank')
                            assert_allclose(x, direct_lstsq(a,b,cmplx=0),
                                          rtol=25*np.finfo(a1.dtype).eps,
                                          atol=25*np.finfo(a1.dtype).eps,
                                          err_msg="driver: %s" % lapack_driver)

    def test_random_complex_overdet(self):

        for dtype in COMPLEX_DTYPES:
            for (n,m) in ((20,15), (200,2)):
                    for lapack_driver in TestLstsq.lapack_drivers:
                        for overwrite in (True,False):
                            a = np.asarray(random([n,m]) + 1j*random([n,m]),
                                           dtype=dtype)
                            for i in range(m):
                                a[i,i] = 20*(.1+a[i,i])
                            for i in range(2):
                                b = np.asarray(random([n,3]), dtype=dtype)
                                # Store values in case they are overwritten
                                # later
                                a1 = a.copy()
                                b1 = b.copy()
                                out = lstsq(a1,b1,lapack_driver=lapack_driver,
                                           overwrite_a=overwrite,
                                           overwrite_b=overwrite)
                                x = out[0]
                                r = out[2]
                                assert_(r == m, 'unexpected efficient rank')
                                assert_allclose(x,direct_lstsq(a,b,cmplx=1),
                                          rtol=25*np.finfo(a1.dtype).eps,
                                          atol=25*np.finfo(a1.dtype).eps,
                                          err_msg="driver: %s" % lapack_driver)

    def test_check_finite(self):
        for dtype in REAL_DTYPES:
            a = np.array(((1,20),(-30,4)), dtype=dtype)
            for bt in (((1,0),(0,1)),(1,0),
                      ((2,1),(-30,4))):
                for lapack_driver in TestLstsq.lapack_drivers:
                        for overwrite in (True,False):
                            for check_finite in (True,False):
                                b = np.array(bt,dtype=dtype)
                                # Store values in case they are overwritten
                                # later
                                a1 = a.copy()
                                b1 = b.copy()

                                out = lstsq(a1,b1,lapack_driver=lapack_driver,
                                            check_finite=check_finite,
                                            overwrite_a=overwrite,
                                            overwrite_b=overwrite)
                                x = out[0]
                                r = out[2]
                                assert_(r == 2, 'unexpected efficient rank')
                                assert_allclose(dot(a,x),b,
                                          rtol=25*np.finfo(a.dtype).eps,
                                          atol=25*np.finfo(a.dtype).eps,
                                          err_msg="driver: %s" % lapack_driver)


class TestPinv(TestCase):

    def test_simple_real(self):
        a = array([[1, 2, 3], [4, 5, 6], [7, 8, 10]], dtype=float)
        a_pinv = pinv(a)
        assert_array_almost_equal(dot(a,a_pinv), np.eye(3))
        a_pinv = pinv2(a)
        assert_array_almost_equal(dot(a,a_pinv), np.eye(3))

    def test_simple_complex(self):
        a = (array([[1, 2, 3], [4, 5, 6], [7, 8, 10]], dtype=float)
             + 1j * array([[10, 8, 7], [6, 5, 4], [3, 2, 1]], dtype=float))
        a_pinv = pinv(a)
        assert_array_almost_equal(dot(a, a_pinv), np.eye(3))
        a_pinv = pinv2(a)
        assert_array_almost_equal(dot(a, a_pinv), np.eye(3))

    def test_simple_singular(self):
        a = array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=float)
        a_pinv = pinv(a)
        a_pinv2 = pinv2(a)
        assert_array_almost_equal(a_pinv,a_pinv2)

    def test_simple_cols(self):
        a = array([[1, 2, 3], [4, 5, 6]], dtype=float)
        a_pinv = pinv(a)
        a_pinv2 = pinv2(a)
        assert_array_almost_equal(a_pinv,a_pinv2)

    def test_simple_rows(self):
        a = array([[1, 2], [3, 4], [5, 6]], dtype=float)
        a_pinv = pinv(a)
        a_pinv2 = pinv2(a)
        assert_array_almost_equal(a_pinv,a_pinv2)

    def test_check_finite(self):
        a = array([[1,2,3],[4,5,6.],[7,8,10]])
        a_pinv = pinv(a, check_finite=False)
        assert_array_almost_equal(dot(a,a_pinv),[[1,0,0],[0,1,0],[0,0,1]])
        a_pinv = pinv2(a, check_finite=False)
        assert_array_almost_equal(dot(a,a_pinv),[[1,0,0],[0,1,0],[0,0,1]])


class TestPinvSymmetric(TestCase):

    def test_simple_real(self):
        a = array([[1, 2, 3], [4, 5, 6], [7, 8, 10]], dtype=float)
        a = np.dot(a, a.T)
        a_pinv = pinvh(a)
        assert_array_almost_equal(np.dot(a, a_pinv), np.eye(3))

    def test_nonpositive(self):
        a = array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=float)
        a = np.dot(a, a.T)
        u, s, vt = np.linalg.svd(a)
        s[0] *= -1
        a = np.dot(u * s, vt)  # a is now symmetric non-positive and singular
        a_pinv = pinv2(a)
        a_pinvh = pinvh(a)
        assert_array_almost_equal(a_pinv, a_pinvh)

    def test_simple_complex(self):
        a = (array([[1, 2, 3], [4, 5, 6], [7, 8, 10]], dtype=float)
             + 1j * array([[10, 8, 7], [6, 5, 4], [3, 2, 1]], dtype=float))
        a = np.dot(a, a.conj().T)
        a_pinv = pinvh(a)
        assert_array_almost_equal(np.dot(a, a_pinv), np.eye(3))


class TestVectorNorms(object):

    def test_types(self):
        for dtype in np.typecodes['AllFloat']:
            x = np.array([1,2,3], dtype=dtype)
            tol = max(1e-15, np.finfo(dtype).eps.real * 20)
            assert_allclose(norm(x), np.sqrt(14), rtol=tol)
            assert_allclose(norm(x, 2), np.sqrt(14), rtol=tol)

        for dtype in np.typecodes['Complex']:
            x = np.array([1j,2j,3j], dtype=dtype)
            tol = max(1e-15, np.finfo(dtype).eps.real * 20)
            assert_allclose(norm(x), np.sqrt(14), rtol=tol)
            assert_allclose(norm(x, 2), np.sqrt(14), rtol=tol)

    def test_overflow(self):
        # unlike numpy's norm, this one is
        # safer on overflow
        a = array([1e20], dtype=float32)
        assert_almost_equal(norm(a), a)

    def test_stable(self):
        # more stable than numpy's norm
        a = array([1e4] + [1]*10000, dtype=float32)
        try:
            # snrm in double precision; we obtain the same as for float64
            # -- large atol needed due to varying blas implementations
            assert_allclose(norm(a) - 1e4, 0.5, atol=1e-2)
        except AssertionError:
            # snrm implemented in single precision, == np.linalg.norm result
            msg = ": Result should equal either 0.0 or 0.5 (depending on " \
                  "implementation of snrm2)."
            assert_almost_equal(norm(a) - 1e4, 0.0, err_msg=msg)

    def test_zero_norm(self):
        assert_equal(norm([1,0,3], 0), 2)
        assert_equal(norm([1,2,3], 0), 3)


class TestMatrixNorms(object):

    def test_matrix_norms(self):
        # Not all of these are matrix norms in the most technical sense.
        np.random.seed(1234)
        for n, m in (1, 1), (1, 3), (3, 1), (4, 4), (4, 5), (5, 4):
            for t in np.single, np.double, np.csingle, np.cdouble, np.int64:
                A = 10 * np.random.randn(n, m).astype(t)
                if np.issubdtype(A.dtype, np.complexfloating):
                    A = (A + 10j * np.random.randn(n, m)).astype(t)
                    t_high = np.cdouble
                else:
                    t_high = np.double
                for order in (None, 'fro', 1, -1, 2, -2, np.inf, -np.inf):
                    actual = norm(A, ord=order)
                    desired = np.linalg.norm(A, ord=order)
                    # SciPy may return higher precision matrix norms.
                    # This is a consequence of using LAPACK.
                    if not np.allclose(actual, desired):
                        desired = np.linalg.norm(A.astype(t_high), ord=order)
                        np.assert_allclose(actual, desired)


class TestOverwrite(object):
    def test_solve(self):
        assert_no_overwrite(solve, [(3,3), (3,)])

    def test_solve_triangular(self):
        assert_no_overwrite(solve_triangular, [(3,3), (3,)])

    def test_solve_banded(self):
        assert_no_overwrite(lambda ab, b: solve_banded((2,1), ab, b),
                            [(4,6), (6,)])

    def test_solveh_banded(self):
        assert_no_overwrite(solveh_banded, [(2,6), (6,)])

    def test_inv(self):
        assert_no_overwrite(inv, [(3,3)])

    def test_det(self):
        assert_no_overwrite(det, [(3,3)])

    def test_lstsq(self):
        assert_no_overwrite(lstsq, [(3,2), (3,)])

    def test_pinv(self):
        assert_no_overwrite(pinv, [(3,3)])

    def test_pinv2(self):
        assert_no_overwrite(pinv2, [(3,3)])

    def test_pinvh(self):
        assert_no_overwrite(pinvh, [(3,3)])


class TestSolveCirculant(TestCase):

    def test_basic1(self):
        c = np.array([1, 2, 3, 5])
        b = np.array([1, -1, 1, 0])
        x = solve_circulant(c, b)
        y = solve(circulant(c), b)
        assert_allclose(x, y)

    def test_basic2(self):
        # b is a 2-d matrix.
        c = np.array([1, 2, -3, -5])
        b = np.arange(12).reshape(4, 3)
        x = solve_circulant(c, b)
        y = solve(circulant(c), b)
        assert_allclose(x, y)

    def test_basic3(self):
        # b is a 3-d matrix.
        c = np.array([1, 2, -3, -5])
        b = np.arange(24).reshape(4, 3, 2)
        x = solve_circulant(c, b)
        y = solve(circulant(c), b)
        assert_allclose(x, y)

    def test_complex(self):
        # Complex b and c
        c = np.array([1+2j, -3, 4j, 5])
        b = np.arange(8).reshape(4, 2) + 0.5j
        x = solve_circulant(c, b)
        y = solve(circulant(c), b)
        assert_allclose(x, y)

    def test_random_b_and_c(self):
        # Random b and c
        np.random.seed(54321)
        c = np.random.randn(50)
        b = np.random.randn(50)
        x = solve_circulant(c, b)
        y = solve(circulant(c), b)
        assert_allclose(x, y)

    def test_singular(self):
        # c gives a singular circulant matrix.
        c = np.array([1, 1, 0, 0])
        b = np.array([1, 2, 3, 4])
        x = solve_circulant(c, b, singular='lstsq')
        y, res, rnk, s = lstsq(circulant(c), b)
        assert_allclose(x, y)
        assert_raises(LinAlgError, solve_circulant, x, y)

    def test_axis_args(self):
        # Test use of caxis, baxis and outaxis.

        # c has shape (2, 1, 4)
        c = np.array([[[-1, 2.5, 3, 3.5]], [[1, 6, 6, 6.5]]])

        # b has shape (3, 4)
        b = np.array([[0, 0, 1, 1], [1, 1, 0, 0], [1, -1, 0, 0]])

        x = solve_circulant(c, b, baxis=1)
        assert_equal(x.shape, (4, 2, 3))
        expected = np.empty_like(x)
        expected[:, 0, :] = solve(circulant(c[0]), b.T)
        expected[:, 1, :] = solve(circulant(c[1]), b.T)
        assert_allclose(x, expected)

        x = solve_circulant(c, b, baxis=1, outaxis=-1)
        assert_equal(x.shape, (2, 3, 4))
        assert_allclose(np.rollaxis(x, -1), expected)

        # np.swapaxes(c, 1, 2) has shape (2, 4, 1); b.T has shape (4, 3).
        x = solve_circulant(np.swapaxes(c, 1, 2), b.T, caxis=1)
        assert_equal(x.shape, (4, 2, 3))
        assert_allclose(x, expected)


if __name__ == "__main__":
    run_module_suite()
