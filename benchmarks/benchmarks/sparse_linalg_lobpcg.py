from functools import partial

import numpy as np
from .common import Benchmark, safe_import

with safe_import():
    from scipy.linalg import eigh, orth, cho_factor, cho_solve
    import scipy.sparse
    from scipy.sparse.linalg import lobpcg
    from scipy.sparse.linalg._interface import LinearOperator


def _sakurai(n):
    """ Example taken from
        T. Sakurai, H. Tadano, Y. Inadomi and U. Nagashima
        A moment-based method for large-scale generalized eigenvalue problems
        Appl. Num. Anal. Comp. Math. Vol. 1 No. 2 (2004) """

    A = scipy.sparse.eye(n, n)
    d0 = np.array(np.r_[5, 6 * np.ones(n-2), 5])
    d1 = -4 * np.ones(n)
    d2 = np.ones(n)
    B = scipy.sparse.spdiags([d2, d1, d0, d1, d2], [-2, -1, 0, 1, 2], n, n)

    k = np.arange(1, n+1)
    w_ex = np.sort(1. / (16.* np.power(np.cos(0.5*k*np.pi/(n+1)), 4)))  # exact eigenvalues

    return A, B, w_ex


def _mikota_pair(n):
    # Mikota pair acts as a nice test since the eigenvalues
    # are the squares of the integers n, n=1,2,...
    x = np.arange(1, n + 1)
    B = np.diag(1. / x)
    y = np.arange(n - 1, 0, -1)
    z = np.arange(2 * n - 1, 0, -2)
    A = np.diag(z) - np.diag(y, -1) - np.diag(y, 1)
    return A.astype(float), B.astype(float)


def _as2d(ar):
    if ar.ndim == 2:
        return ar
    else:  # Assume 1!
        aux = np.array(ar, copy=False)
        aux.shape = (ar.shape[0], 1)
        return aux


def _precond(LorU, lower, x):
    y = cho_solve((LorU, lower), x)
    return _as2d(y)


class Bench(Benchmark):
    params = [
        [],
        ['lobpcg', 'eigh']
    ]
    param_names = ['n', 'solver']

    def __init__(self):
        self.time_mikota.__func__.params = list(self.params)
        self.time_mikota.__func__.params[0] = [128, 256, 512, 1024, 2048]
        self.time_mikota.__func__.setup = self.setup_mikota

        self.time_sakurai.__func__.params = list(self.params)
        self.time_sakurai.__func__.params[0] = [50, 400]
        self.time_sakurai.__func__.setup = self.setup_sakurai

    def setup_mikota(self, n, solver):
        self.shape = (n, n)
        self.A, self.B = _mikota_pair(n)

        if solver == 'eigh' and n >= 512:
            # skip: slow, and not useful to benchmark
            raise NotImplementedError()

    def setup_sakurai(self, n, solver):
        self.shape = (n, n)
        self.A, self.B, all_eigenvalues = _sakurai(n)
        self.A_dense = self.A.A
        self.B_dense = self.B.A

    def time_mikota(self, n, solver):
        m = 10
        if solver == 'lobpcg':
            X = np.random.rand(n, m)
            X = orth(X)
            LorU, lower = cho_factor(self.A, lower=0, overwrite_a=0)
            M = LinearOperator(self.shape,
                               matvec=partial(_precond, LorU, lower),
                               matmat=partial(_precond, LorU, lower))
            eigs, vecs = lobpcg(self.A, X, self.B, M, tol=1e-4, maxiter=40)
        else:
            eigh(self.A, self.B, eigvals_only=True, eigvals=(0, m - 1))

    def time_sakurai(self, n, solver):
        m = 3
        if solver == 'lobpcg':
            X = np.random.rand(n, m)
            eigs, vecs, resnh = lobpcg(self.A, X, self.B, tol=1e-6, maxiter=500,
                                       retResidualNormsHistory=1)
        else:
            eigh(self.A_dense, self.B_dense, eigvals_only=True, eigvals=(0, m - 1))

    # Retain old benchmark results (remove this if changing the benchmark)
    time_mikota.version = "a1fb679758f7e5cf79d18cc4930afdff999fccc142fe7a4f63e73b39ab1f58bb"
    time_sakurai.version = "7c38d449924fb71f777bd408072ecc883b8b05e53a6544e97da3887fbc10b235"
