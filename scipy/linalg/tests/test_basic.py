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

from numpy import arange, add, array, dot, zeros, identity, conjugate, transpose
import numpy.linalg as linalg

from numpy.testing import *

from scipy.linalg import solve,inv,det,lstsq, toeplitz, hankel, tri, triu, \
     tril, pinv, pinv2, solve_banded


def random(size):
    return rand(*size)

def get_mat(n):
    data = arange(n)
    data = add.outer(data,data)
    return data

class TestSolveBanded(TestCase):

    def test_simple(self):

        a = [[1,20,0,0],[-30,4,6,0],[2,1,20,2],[0,-1,7,14]]
        ab = [[0,20,6,2],
              [1,4,20,14],
              [-30,1,7,0],
              [2,-1,0,0]]
        l,u = 2,1
        for b in ([[1,0,0,0],[0,0,0,1],[0,1,0,0],[0,1,0,0]],
                  [[2,1],[-30,4],[2,3],[1,3]]):
            x = solve_banded((l,u),ab,b)
            assert_array_almost_equal(dot(a,x),b)

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

class TestTri(TestCase):
    def test_basic(self):
        assert_equal(tri(4),array([[1,0,0,0],
                                   [1,1,0,0],
                                   [1,1,1,0],
                                   [1,1,1,1]]))
        assert_equal(tri(4,dtype='f'),array([[1,0,0,0],
                                                [1,1,0,0],
                                                [1,1,1,0],
                                                [1,1,1,1]],'f'))
    def test_diag(self):
        assert_equal(tri(4,k=1),array([[1,1,0,0],
                                       [1,1,1,0],
                                       [1,1,1,1],
                                       [1,1,1,1]]))
        assert_equal(tri(4,k=-1),array([[0,0,0,0],
                                        [1,0,0,0],
                                        [1,1,0,0],
                                        [1,1,1,0]]))
    def test_2d(self):
        assert_equal(tri(4,3),array([[1,0,0],
                                     [1,1,0],
                                     [1,1,1],
                                     [1,1,1]]))
        assert_equal(tri(3,4),array([[1,0,0,0],
                                     [1,1,0,0],
                                     [1,1,1,0]]))
    def test_diag2d(self):
        assert_equal(tri(3,4,k=2),array([[1,1,1,0],
                                         [1,1,1,1],
                                         [1,1,1,1]]))
        assert_equal(tri(4,3,k=-2),array([[0,0,0],
                                          [0,0,0],
                                          [1,0,0],
                                          [1,1,0]]))

class TestTril(TestCase):
    def test_basic(self):
        a = (100*get_mat(5)).astype('l')
        b = a.copy()
        for k in range(5):
            for l in range(k+1,5):
                b[k,l] = 0
        assert_equal(tril(a),b)

    def test_diag(self):
        a = (100*get_mat(5)).astype('f')
        b = a.copy()
        for k in range(5):
            for l in range(k+3,5):
                b[k,l] = 0
        assert_equal(tril(a,k=2),b)
        b = a.copy()
        for k in range(5):
            for l in range(max((k-1,0)),5):
                b[k,l] = 0
        assert_equal(tril(a,k=-2),b)

class TestTriu(TestCase):
    def test_basic(self):
        a = (100*get_mat(5)).astype('l')
        b = a.copy()
        for k in range(5):
            for l in range(k+1,5):
                b[l,k] = 0
        assert_equal(triu(a),b)

    def test_diag(self):
        a = (100*get_mat(5)).astype('f')
        b = a.copy()
        for k in range(5):
            for l in range(max((k-1,0)),5):
                b[l,k] = 0
        assert_equal(triu(a,k=2),b)
        b = a.copy()
        for k in range(5):
            for l in range(k+3,5):
                b[l,k] = 0
        assert_equal(triu(a,k=-2),b)

class TestToeplitz(TestCase):
    def test_basic(self):
        y = toeplitz([1,2,3])
        assert_array_equal(y,[[1,2,3],[2,1,2],[3,2,1]])
        y = toeplitz([1,2,3],[1,4,5])
        assert_array_equal(y,[[1,4,5],[2,1,4],[3,2,1]])

class TestHankel(TestCase):
    def test_basic(self):
        y = hankel([1,2,3])
        assert_array_equal(y,[[1,2,3],[2,3,0],[3,0,0]])
        y = hankel([1,2,3],[3,4,5])
        assert_array_equal(y,[[1,2,3],[2,3,4],[3,4,5]])

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

if __name__ == "__main__":
    run_module_suite()
