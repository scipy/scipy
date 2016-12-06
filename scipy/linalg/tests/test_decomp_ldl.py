from __future__ import division, print_function, absolute_import

from numpy.testing import (TestCase, assert_array_almost_equal,
                           assert_raises, assert_)
from numpy import array, eye, zeros_like, any, tril_indices_from
from numpy.random import rand, randint, seed
from scipy.linalg import ldl


class TestLDL(TestCase):
    def setUp(self):
        seed(1234)

    def Test_args(self):
        A = eye(3)
        assert_raises(ValueError, ldl, A[:, :2])
        assert_raises(ValueError, ldl, A*1j) 

    def Test_simple(self):
        a = array([[-0.39-0.71j, 5.14-0.64j, -7.86-2.96j, 3.80+0.92j],
                   [5.14-0.64j, 8.86+1.81j, -3.52+0.58j, 5.32-1.59j],
                   [-7.86-2.96j, -3.52+0.58j, -2.83-0.03j, -1.54-2.86j],
                   [3.80+0.92j, 5.32-1.59j, -1.54-2.86j, -0.56+0.12j]])
        b = array([[5, 10, 1, 18],
                   [10, 2, 11, 1],
                   [1, 11, 19, 9],
                   [18, 1, 9, 0]])
        c = array([[52, 97, 112, 107, 50],
                   [97, 114, 89, 98, 13],
                   [112, 89, 64, 33, 6],
                   [107, 98, 33, 60, 73],
                   [50, 13, 6, 73, 77]])

        d = array([[2, 2, -4, 0, 4],
                   [2, -2, -2, 10, -8],
                   [-4, -2, 6, -8, -4],
                   [0, 10, -8, 6, -6],
                   [4, -8, -4, -6, 10]])
        for x in (b, c, d):
            l, d, p = ldl(x)
            assert_array_almost_equal(l.dot(d).dot(l.T)-x, zeros_like(x))
            l, d, p = ldl(x, lower=False)
            assert_array_almost_equal(l.dot(d).dot(l.T)-x, zeros_like(x))

        l, d, p = ldl(a, only_sym=True)
        assert_array_almost_equal(l.dot(d).dot(l.T)-a, zeros_like(a))

    def Test_rook_pivoting(self):
        a = array([[-0.39-0.71j, 5.14-0.64j, -7.86-2.96j, 3.80+0.92j],
                   [5.14-0.64j, 8.86+1.81j, -3.52+0.58j, 5.32-1.59j],
                   [-7.86-2.96j, -3.52+0.58j, -2.83-0.03j, -1.54-2.86j],
                   [3.80+0.92j, 5.32-1.59j, -1.54-2.86j, -0.56+0.12j]])
        b = array([[5, 10, 1, 18],
                   [10, 2, 11, 1],
                   [1, 11, 19, 9],
                   [18, 1, 9, 0]])
        c = array([[52, 97, 112, 107, 50],
                   [97, 114, 89, 98, 13],
                   [112, 89, 64, 33, 6],
                   [107, 98, 33, 60, 73],
                   [50, 13, 6, 73, 77]])
        d = array([[2, 2, -4, 0, 4],
                   [2, -2, -2, 10, -8],
                   [-4, -2, 6, -8, -4],
                   [0, 10, -8, 6, -6],
                   [4, -8, -4, -6, 10]])

        for x in (b, c, d):
            l, d, p = ldl(x, rook_pivoting=True)
            assert_array_almost_equal(l.dot(d).dot(l.T)-x, zeros_like(x))
            l, d, p = ldl(x, lower=False, rook_pivoting=True)
            assert_array_almost_equal(l.dot(d).dot(l.T)-x, zeros_like(x))
        l, d, p = ldl(a, only_sym=True, rook_pivoting=True)
        assert_array_almost_equal(l.dot(d).dot(l.T)-a, zeros_like(a))

    def Test_permutations(self):
        for _ in range(10):
            n = randint(1, 100)
            x = rand(n,n)
            x = x + x.conj().T
            x += eye(n)*randint(5, 1e6)
            lu_ind = tril_indices_from(x, k=-1)
            l, d, p = ldl(x, lower=0)
            assert_(not any(l[p, :][lu_ind]), 'Spin {} failed'.format(_))
