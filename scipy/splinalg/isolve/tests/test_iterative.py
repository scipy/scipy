#!/usr/bin/env python
""" Test functions for the splinalg.isolve.iterative module
"""

import sys

from numpy import zeros, dot, diag, ones
from scipy.testing import *
from numpy.random import rand
#from numpy import arange, add, array, dot, zeros, identity, conjugate, transpose

from scipy.linalg import norm
from scipy.splinalg.isolve import cg, cgs, bicg, bicgstab, gmres, qmr


def callback(x):
    global A, b
    res = b-dot(A,x)
    #print "||A.x - b|| = " + str(norm(dot(A,x)-b))

class TestIterativeSolvers(TestCase):
    def __init__(self, *args, **kwds):
        TestCase.__init__(self, *args, **kwds)
        self.setUp()
    def setUp (self):
        global A, b
        n = 10
        self.tol = 1e-5
        self.x0 = zeros(n, float)
        A = rand(n, n)+diag(4*ones(n))
        self.A = 0.5 * (A+A.T)
        A = self.A
        self.b = rand(n)
        b = self.b

    def test_cg(self):
        bx0 = self.x0.copy()
        x, info = cg(self.A, self.b, self.x0, callback=callback)
        assert_array_equal(bx0, self.x0)
        assert norm(dot(self.A, x) - self.b) < 5*self.tol

    def test_bicg(self):
        bx0 = self.x0.copy()
        x, info = bicg(self.A, self.b, self.x0, callback=callback)
        assert_array_equal(bx0, self.x0)
        assert norm(dot(self.A, x) - self.b) < 5*self.tol

    def test_cgs(self):
        bx0 = self.x0.copy()
        x, info = cgs(self.A, self.b, self.x0, callback=callback)
        assert_array_equal(bx0, self.x0)
        assert norm(dot(self.A, x) - self.b) < 5*self.tol

    def test_bicgstab(self):
        bx0 = self.x0.copy()
        x, info = bicgstab(self.A, self.b, self.x0, callback=callback)
        assert_array_equal(bx0, self.x0)
        assert norm(dot(self.A, x) - self.b) < 5*self.tol

    def test_gmres(self):
        bx0 = self.x0.copy()
        x, info = gmres(self.A, self.b, self.x0, callback=callback)
        assert_array_equal(bx0, self.x0)
        assert norm(dot(self.A, x) - self.b) < 5*self.tol

    def test_qmr(self):
        bx0 = self.x0.copy()
        x, info = qmr(self.A, self.b, self.x0, callback=callback)
        assert_array_equal(bx0, self.x0)
        assert norm(dot(self.A, x) - self.b) < 5*self.tol

if __name__ == "__main__":
    nose.run(argv=['', __file__])
