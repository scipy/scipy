#!/usr/bin/env python
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
import unittest
from math import sqrt
from numpy import arange, concatenate, identity, zeros, ones
from numpy import transpose
from numpy.linalg import norm
from scipy.sparse import coo_matrix
from scipy.sparse.linalg.interface import aslinearoperator
from lsmr import lsmr

class lsmrTestCase(unittest.TestCase):
    def setUp(self):
        self.n = 10
        self.m = 10

    def assertCompatibleSystem(self, A, xtrue):
        Afun = aslinearoperator(A)
        b = Afun.matvec(xtrue)
        x = lsmr(A,b)[0]
        self.assertAlmostEqual(norm(x - xtrue), 0, 6)

    def testIdentityACase1(self):
        A = identity(self.n)
        xtrue = zeros((self.n, 1))
        self.assertCompatibleSystem(A, xtrue)

    def testIdentityACase2(self):
        A = identity(self.n)
        xtrue = ones((self.n,1))
        self.assertCompatibleSystem(A, xtrue)

    def testIdentityACase3(self):
        A = identity(self.n)
        xtrue = transpose(arange(self.n,0,-1))
        self.assertCompatibleSystem(A, xtrue)

    def testBidiagonalA(self):
        A = lowerBidiagonalMatrix(20,self.n)
        xtrue = transpose(arange(self.n,0,-1))
        self.assertCompatibleSystem(A,xtrue)

class lsmrReturnsTestCase(unittest.TestCase):
    def setUp(self):
        self.n = 10
        self.A = lowerBidiagonalMatrix(20,self.n)
        self.xtrue = transpose(arange(self.n,0,-1))
        self.Afun = aslinearoperator(self.A)
        self.b = self.Afun.matvec(self.xtrue)
        self.returnValues = lsmr(self.A,self.b)

    def testNormr(self):
        x, istop, itn, normr, normar, normA, condA, normx = self.returnValues;
        self.assertAlmostEqual(normr, norm(self.b - self.Afun.matvec(x)))

    def testNormar(self):
        x, istop, itn, normr, normar, normA, condA, normx = self.returnValues;
        self.assertAlmostEqual(normar, \
                norm(self.Afun.rmatvec(self.b - self.Afun.matvec(x))))

    def testNormx(self):
        x, istop, itn, normr, normar, normA, condA, normx = self.returnValues;
        self.assertAlmostEqual(normx, norm(x))

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
    row = concatenate((arange(m, dtype='intc'), \
        arange(1, m, dtype='intc')), axis = 1)
    col = concatenate((arange(m, dtype='intc'), \
        arange(m-1, dtype='intc')), axis = 1)
    data = concatenate((arange(1, m+1, dtype='double'), \
            arange(1,m, dtype='double')), axis = 1)
    return coo_matrix((data,(row, col)), shape=(m,n))
  else:
    row = concatenate((arange(n, dtype='intc'), \
        arange(1, n+1, dtype='intc')), axis = 1)
    col = concatenate((arange(n, dtype='intc'), \
        arange(n, dtype='intc')), axis = 1)
    data = concatenate((arange(1, n+1, dtype='double'), \
            arange(1,n+1, dtype='double')), axis = 1)
    return coo_matrix((data,(row, col)), shape=(m,n))

def lsmrtest(m, n, damp):
    """Verbose testing of lsmr"""

    A = lowerBidiagonalMatrix(m,n)
    xtrue = arange(n,0,-1, dtype='double')
    Afun = aslinearoperator(A)

    b = Afun.matvec(xtrue)

    atol      = 1.0e-7;
    btol      = 1.0e-7;
    conlim    = 1.0e+10;
    itnlim    = 10*n;
    show      = 1;

    x, istop, itn, normr, normar, norma, conda, normx \
      = lsmr(A, b, damp, atol, btol, conlim, itnlim, show )

    j1 = min(n,5);   j2 = max(n-4,1);
    print ' '
    print 'First elements of x:'
    str = [ '%10.4f' %(xi) for xi in x[0:j1] ]
    print ''.join(str)
    print ' '
    print 'Last  elements of x:'
    str = [ '%10.4f' %(xi) for xi in x[j2-1:] ]
    print ''.join(str)

    r    = b - Afun.matvec(x);
    r2   = sqrt(norm(r)**2 + (damp*norm(x))**2)
    print ' '
    str =  'normr (est.)  %17.10e' %(normr )
    str2 =  'normr (true)  %17.10e' %(r2 )
    print str
    print str2
    print ' '

if __name__ == "__main__":
    # Comment out the next line to run unit tests only
    lsmrtest(20,10,0)
    print("Continuing to run some unit tests...")
    unittest.main()
