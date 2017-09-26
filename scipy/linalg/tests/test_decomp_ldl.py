from __future__ import division, absolute_import

from numpy.testing import assert_array_almost_equal, assert_
from numpy import (array, eye, zeros_like, zeros, tril_indices_from, tril,
                   empty_like, empty)
from numpy.random import rand, randint, seed
from scipy.linalg import ldl
from pytest import raises as assert_raises
from warnings import catch_warnings, simplefilter
from numpy import ComplexWarning


def test_args():
    A = eye(3)
    # Nonsquare array
    assert_raises(ValueError, ldl, A[:, :2])
    # Complex matrix with imaginary diagonal entries with "hermitian=True"
    with catch_warnings():
        simplefilter('error')
        assert_raises(ComplexWarning, ldl, A*1j)


def test_empty_array():
    a = empty((0, 0), dtype=complex)
    l, d, p = ldl(empty((0, 0)))
    assert_array_almost_equal(l, empty_like(a))
    assert_array_almost_equal(d, empty_like(a))
    assert_array_almost_equal(p, array([], dtype=int))


def test_simple():
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
    e = array([[-1.36+0.00j, 0+0j, 0+0j, 0+0j],
               [1.58-0.90j, -8.87+0j, 0+0j, 0+0j],
               [2.21+0.21j, -1.84+0.03j, -4.63+0j, 0+0j],
               [3.91-1.50j, -1.78-1.18j, 0.11-0.11j, -1.84+0.00j]])
    for x in (b, c, d):
        l, d, p = ldl(x)
        assert_array_almost_equal(l.dot(d).dot(l.T)-x, zeros_like(x))
        l, d, p = ldl(x, lower=False)
        assert_array_almost_equal(l.dot(d).dot(l.T)-x, zeros_like(x))

    l, d, p = ldl(a, hermitian=False)
    assert_array_almost_equal(l.dot(d).dot(l.T)-a, zeros_like(a))

    # Use upper part for the computation and use the lower part for comparison
    l, d, p = ldl(e.conj().T, lower=0)
    assert_array_almost_equal(tril(l.dot(d).dot(l.conj().T)-e), zeros((4, 4)))


def test_permutations():
    seed(1234)
    for _ in range(10):
        n = randint(1, 100)
        x = rand(n, n)
        x = x + x.conj().T
        x += eye(n)*randint(5, 1e6)
        lu_ind = tril_indices_from(x, k=-1)
        # Test whether permutations lead to a triangular array
        l, d, p = ldl(x, lower=0)
        assert_(not any(l[p, :][lu_ind]), 'Spin {} failed'.format(_))
