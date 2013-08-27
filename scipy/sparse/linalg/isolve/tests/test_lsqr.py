from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_
from scipy.lib.six import xrange

import scipy.sparse
import scipy.sparse.linalg
from scipy.sparse.linalg import lsqr
from time import time

# Set up a test problem
n = 35
G = np.eye(n)
normal = np.random.normal
norm = np.linalg.norm

for jj in xrange(5):
    gg = normal(size=n)
    hh = gg * gg.T
    G += (hh + hh.T) * 0.5
    G += normal(size=n) * normal(size=n)

b = normal(size=n)

tol = 1e-10
show = False
maxit = None


def test_basic():
    svx = np.linalg.solve(G, b)
    X = lsqr(G, b, show=show, atol=tol, btol=tol, iter_lim=maxit)
    xo = X[0]
    assert_(norm(svx - xo) < 1e-5)

if __name__ == "__main__":
    svx = np.linalg.solve(G, b)

    tic = time()
    X = lsqr(G, b, show=show, atol=tol, btol=tol, iter_lim=maxit)
    xo = X[0]
    phio = X[3]
    psio = X[7]
    k = X[2]
    chio = X[8]
    mg = np.amax(G - G.T)
    if mg > 1e-14:
        sym = 'No'
    else:
        sym = 'Yes'

    print('LSQR')
    print("Is linear operator symmetric? " + sym)
    print("n: %3g  iterations:   %3g" % (n, k))
    print("Norms computed in %.2fs by LSQR" % (time() - tic))
    print(" ||x||  %9.4e  ||r|| %9.4e  ||Ar||  %9.4e " % (chio, phio, psio))
    print("Residual norms computed directly:")
    print(" ||x||  %9.4e  ||r|| %9.4e  ||Ar||  %9.4e" % (norm(xo),
                                                          norm(G*xo - b),
                                                          norm(G.T*(G*xo-b))))
    print("Direct solution norms:")
    print(" ||x||  %9.4e  ||r|| %9.4e " % (norm(svx), norm(G*svx - b)))
    print("")
    print(" || x_{direct} - x_{LSQR}|| %9.4e " % norm(svx-xo))
    print("")


def test_gh_2466():
    row = np.array([0, 0])
    col = np.array([0, 1])
    val = np.array([1, -1])
    A = scipy.sparse.coo_matrix((val, (row, col)), shape=(1, 2))
    b = np.asarray([4])
    lsqr(A, b)
