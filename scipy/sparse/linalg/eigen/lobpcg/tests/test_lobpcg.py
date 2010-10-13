#!/usr/bin/env python
""" Test functions for the sparse.linalg.eigen.lobpcg module
"""

import numpy
from numpy.testing import assert_almost_equal, run_module_suite

from scipy import arange, ones, rand, set_printoptions, r_, diag, linalg, eye
from scipy.linalg import eig
from scipy.sparse.linalg.eigen.lobpcg import lobpcg


set_printoptions(precision=3,linewidth=90)



def ElasticRod(n):
    # Fixed-free elastic rod
    L = 1.0
    le=L/n
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

    numpy.random.seed(0)

    V = rand(n,m)
    X = linalg.orth(V)

    eigs,vecs = lobpcg(A, X, B=B, tol=1e-5, maxiter=30)
    eigs.sort()

    #w,v = symeig(A,B)
    w,v = eig(A,b=B)
    w.sort()

    assert_almost_equal(w[:m/2],eigs[:m/2],decimal=2)

    #from pylab import plot, show, legend, xlabel, ylabel
    #plot(arange(0,len(w[:m])),w[:m],'bx',label='Results by symeig')
    #plot(arange(0,len(eigs)),eigs,'r+',label='Results by lobpcg')
    #legend()
    #xlabel(r'Eigenvalue $i$')
    #ylabel(r'$\lambda_i$')
    #show()

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

if __name__ == "__main__":
    run_module_suite()
