"""
Copyright (C) 2010 David Fong and Michael Saunders
Distributed under the same license as Scipy

Testing Code for LSMR.

03 Jun 2010: First version release with lsmr.py

David Chin-lung Fong            clfong@stanford.edu
Institute for Computational and Mathematical Engineering
Stanford University

Michael Saunders                saunders@stanford.edu
Systems Optimization Laboratory
Dept of MS&E, Stanford University.

"""

from __future__ import division, print_function, absolute_import

from numpy import arange, concatenate, eye, zeros, ones, sqrt, \
                  transpose, hstack
from numpy.linalg import norm
from numpy.testing import run_module_suite, assert_almost_equal

from scipy.sparse import coo_matrix
from scipy.sparse.linalg.interface import aslinearoperator
from scipy.sparse.linalg import lsmr


class TestLSMR:
    def setUp(self):
        self.n = 10
        self.m = 10

    def assertCompatibleSystem(self, A, xtrue):
        Afun = aslinearoperator(A)
        b = Afun.matvec(xtrue)
        x = lsmr(A,b)[0]
        assert_almost_equal(norm(x - xtrue), 0, 6)

    def testIdentityACase1(self):
        A = eye(self.n)
        xtrue = zeros((self.n, 1))
        self.assertCompatibleSystem(A, xtrue)

    def testIdentityACase2(self):
        A = eye(self.n)
        xtrue = ones((self.n,1))
        self.assertCompatibleSystem(A, xtrue)

    def testIdentityACase3(self):
        A = eye(self.n)
        xtrue = transpose(arange(self.n,0,-1))
        self.assertCompatibleSystem(A, xtrue)

    def testBidiagonalA(self):
        A = lowerBidiagonalMatrix(20,self.n)
        xtrue = transpose(arange(self.n,0,-1))
        self.assertCompatibleSystem(A,xtrue)


class TestLSMRReturns:
    def setUp(self):
        self.n = 10
        self.A = lowerBidiagonalMatrix(20,self.n)
        self.xtrue = transpose(arange(self.n,0,-1))
        self.Afun = aslinearoperator(self.A)
        self.b = self.Afun.matvec(self.xtrue)
        self.returnValues = lsmr(self.A,self.b)

    def testNormr(self):
        x, istop, itn, normr, normar, normA, condA, normx = self.returnValues
        assert_almost_equal(normr, norm(self.b - self.Afun.matvec(x)))

    def testNormar(self):
        x, istop, itn, normr, normar, normA, condA, normx = self.returnValues
        assert_almost_equal(normar,
                norm(self.Afun.rmatvec(self.b - self.Afun.matvec(x))))

    def testNormx(self):
        x, istop, itn, normr, normar, normA, condA, normx = self.returnValues
        assert_almost_equal(normx, norm(x))


def lowerBidiagonalMatrix(m, n):
    # This is a simple example for testing LSMR.
    # It uses the leading m*n submatrix from
    # A = [ 1
    #       1 2
    #         2 3
    #           3 4
    #             ...
    #               n ]
    # suitably padded by zeros.
    #
    # 04 Jun 2010: First version for distribution with lsmr.py
    if m <= n:
        row = hstack((arange(m, dtype=int),
                      arange(1, m, dtype=int)))
        col = hstack((arange(m, dtype=int),
                      arange(m-1, dtype=int)))
        data = hstack((arange(1, m+1, dtype=float),
                       arange(1,m, dtype=float)))
        return coo_matrix((data, (row, col)), shape=(m,n))
    else:
        row = hstack((arange(n, dtype=int),
                      arange(1, n+1, dtype=int)))
        col = hstack((arange(n, dtype=int),
                      arange(n, dtype=int)))
        data = hstack((arange(1, n+1, dtype=float),
                       arange(1,n+1, dtype=float)))
        return coo_matrix((data,(row, col)), shape=(m,n))


def lsmrtest(m, n, damp):
    """Verbose testing of lsmr"""

    A = lowerBidiagonalMatrix(m,n)
    xtrue = arange(n,0,-1, dtype=float)
    Afun = aslinearoperator(A)

    b = Afun.matvec(xtrue)

    atol = 1.0e-7
    btol = 1.0e-7
    conlim = 1.0e+10
    itnlim = 10*n
    show = 1

    x, istop, itn, normr, normar, norma, conda, normx \
      = lsmr(A, b, damp, atol, btol, conlim, itnlim, show)

    j1 = min(n,5)
    j2 = max(n-4,1)
    print(' ')
    print('First elements of x:')
    str = ['%10.4f' % (xi) for xi in x[0:j1]]
    print(''.join(str))
    print(' ')
    print('Last  elements of x:')
    str = ['%10.4f' % (xi) for xi in x[j2-1:]]
    print(''.join(str))

    r = b - Afun.matvec(x)
    r2 = sqrt(norm(r)**2 + (damp*norm(x))**2)
    print(' ')
    str = 'normr (est.)  %17.10e' % (normr)
    str2 = 'normr (true)  %17.10e' % (r2)
    print(str)
    print(str2)
    print(' ')

if __name__ == "__main__":
    # Comment out the next line to run unit tests only
    lsmrtest(20,10,0)
    run_module_suite()
