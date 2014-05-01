#!/usr/bin/env python
""" Test functions for the sparse.linalg.eigen.lobpcg module
"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (run_module_suite, assert_almost_equal, assert_equal,
        assert_allclose, assert_array_less)

from scipy import arange, ones, rand, r_, diag, linalg, eye
from scipy.linalg import eig, eigh, toeplitz
from scipy.sparse.linalg.eigen.lobpcg import lobpcg


def ElasticRod(n):
    # Fixed-free elastic rod
    L = 1.0
    le = L/n
    rho = 7.85e3
    S = 1.e-4
    E = 2.1e11
    mass = rho*S*le/6.
    k = E*S/le
    A = k*(diag(r_[2.*ones(n-1),1])-diag(ones(n-1),1)-diag(ones(n-1),-1))
    B = mass*(diag(r_[4.*ones(n-1),2])+diag(ones(n-1),1)+diag(ones(n-1),-1))
    return A,B


def MikotaPair(n):
    # Mikota pair acts as a nice test since the eigenvalues
    # are the squares of the integers n, n=1,2,...
    x = arange(1,n+1)
    B = diag(1./x)
    y = arange(n-1,0,-1)
    z = arange(2*n-1,0,-2)
    A = diag(z)-diag(y,-1)-diag(y,1)
    return A,B


def compare_solutions(A,B,m):
    n = A.shape[0]

    np.random.seed(0)

    V = rand(n,m)
    X = linalg.orth(V)

    eigs,vecs = lobpcg(A, X, B=B, tol=1e-5, maxiter=30)
    eigs.sort()

    w,v = eig(A,b=B)
    w.sort()

    assert_almost_equal(w[:int(m/2)],eigs[:int(m/2)],decimal=2)


def test_Small():
    A,B = ElasticRod(10)
    compare_solutions(A,B,10)
    A,B = MikotaPair(10)
    compare_solutions(A,B,10)


def test_ElasticRod():
    A,B = ElasticRod(100)
    compare_solutions(A,B,20)


def test_MikotaPair():
    A,B = MikotaPair(100)
    compare_solutions(A,B,20)


def test_trivial():
    n = 5
    X = ones((n, 1))
    A = eye(n)
    compare_solutions(A, None, n)


def _check_eigen(M, w, V, rtol=1e-8, atol=1e-14):
    assert_allclose(np.multiply(w, V), np.dot(M, V), rtol=rtol, atol=atol)


def _check_fiedler(n, p):
    # This is not necessarily the recommended way to find the Fiedler vector.
    np.random.seed(1234)
    col = np.zeros(n)
    col[1] = 1
    A = toeplitz(col)
    D = np.diag(A.sum(axis=1))
    L = D - A
    # Compute the full eigendecomposition using tricks, e.g.
    # http://www.cs.yale.edu/homes/spielman/561/2009/lect02-09.pdf
    tmp = np.pi * np.arange(n) / n
    analytic_w = 2 * (1 - np.cos(tmp))
    analytic_V = np.cos(np.outer(np.arange(n) + 1/2, tmp))
    _check_eigen(L, analytic_w, analytic_V)
    # Compute the full eigendecomposition using eigh.
    eigh_w, eigh_V = eigh(L)
    _check_eigen(L, eigh_w, eigh_V)
    # Check that the first eigenvalue is near zero and that the rest agree.
    assert_array_less(np.abs([eigh_w[0], analytic_w[0]]), 1e-14)
    assert_allclose(eigh_w[1:], analytic_w[1:])

    # Check small lobpcg eigenvalues.
    X = analytic_V[:, :p]
    lobpcg_w, lobpcg_V = lobpcg(L, X, largest=False)
    assert_equal(lobpcg_w.shape, (p,))
    assert_equal(lobpcg_V.shape, (n, p))
    _check_eigen(L, lobpcg_w, lobpcg_V)
    assert_array_less(np.abs(np.min(lobpcg_w)), 1e-14)
    assert_allclose(np.sort(lobpcg_w)[1:], analytic_w[1:p])

    # Check large lobpcg eigenvalues.
    X = analytic_V[:, -p:]
    lobpcg_w, lobpcg_V = lobpcg(L, X, largest=True)
    assert_equal(lobpcg_w.shape, (p,))
    assert_equal(lobpcg_V.shape, (n, p))
    _check_eigen(L, lobpcg_w, lobpcg_V)
    assert_allclose(np.sort(lobpcg_w), analytic_w[-p:])

    # Look for the Fiedler vector using some legitimate preconditioning.
    preconditioned_fiedler = np.concatenate((np.ones(n//2), -np.ones(n-n//2)))
    X = np.vstack((np.ones(n), preconditioned_fiedler)).T
    lobpcg_w, lobpcg_V = lobpcg(L, X, largest=False)
    # Mathematically, the smaller eigenvalue should be zero
    # and the larger should be the algebraic connectivity.
    lobpcg_w = np.sort(lobpcg_w)
    assert_allclose(lobpcg_w, analytic_w[:2], atol=1e-14)


def test_fiedler_small_8():
    # This triggers the dense path because 8 < 2*5.
    _check_fiedler(8, 2)


def test_fiedler_large_12():
    # This does not trigger the dense path, because 2*5 <= 12.
    _check_fiedler(12, 2)


if __name__ == "__main__":
    run_module_suite()
