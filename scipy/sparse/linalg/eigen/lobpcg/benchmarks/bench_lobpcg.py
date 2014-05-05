from __future__ import division, print_function, absolute_import

from functools import partial
import time

import numpy as np
from numpy.testing import assert_allclose

from scipy import *
from scipy.linalg import eigh, orth, cho_factor, cho_solve
from scipy.sparse.linalg import lobpcg
from scipy.sparse.linalg.interface import LinearOperator


def _mikota_pair(n):
    # Mikota pair acts as a nice test since the eigenvalues
    # are the squares of the integers n, n=1,2,...
    x = arange(1,n+1)
    B = diag(1./x)
    y = arange(n-1,0,-1)
    z = arange(2*n-1,0,-2)
    A = diag(z)-diag(y,-1)-diag(y,1)
    return A,B


def _as2d(ar):
    if ar.ndim == 2:
        return ar
    else:  # Assume 1!
        aux = nm.array(ar, copy=False)
        aux.shape = (ar.shape[0], 1)
        return aux


def _precond(LorU, lower, x):
    y = cho_solve((LorU, lower), x)
    return _as2d(y)


def bench_lobpcg():
    print()
    print('                 lobpcg using mikota pairs')
    print('==============================================================')
    print('      shape      | blocksize |    operation   |   time   ')
    print('                                              | (seconds)')
    print('--------------------------------------------------------------')
    fmt = ' %15s |   %3d     |     %6s     | %6.2f '

    m = 10
    for n in 128, 256, 512, 1024, 2048:
        shape = (n, n)
        A, B = _mikota_pair(n)
        desired_eigenvalues = np.square(np.arange(1, m+1))

        tt = time.clock()
        X = rand(n, m)
        X = orth(X)
        LorU, lower = cho_factor(A, lower=0, overwrite_a=0)
        M = LinearOperator(shape,
                matvec=partial(_precond, LorU, lower),
                matmat=partial(_precond, LorU, lower))
        eigs, vecs = lobpcg(A, X, B, M, tol=1e-4, maxiter=40)
        eigs = sort(eigs)
        elapsed = time.clock() - tt
        assert_allclose(eigs, desired_eigenvalues)
        print(fmt % (shape, m, 'lobpcg', elapsed))

        tt = time.clock()
        w = eigh(A, B, eigvals_only=True, eigvals=(0, m-1))
        elapsed = time.clock() - tt
        assert_allclose(w, desired_eigenvalues)
        print(fmt % (shape, m, 'eigh', elapsed))

