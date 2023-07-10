from functools import partial

import numpy as np
from .common import Benchmark, safe_import

with safe_import():
    from scipy.linalg import eigh, cho_factor, cho_solve
    from scipy.sparse import diags
    from scipy.sparse.linalg import lobpcg, eigsh, LinearOperator


def _sakurai_rev(n):
    """ Example turns a generalized eigenproblem for the matrix pair A and B
        T. Sakurai, H. Tadano, Y. Inadomi and U. Nagashima
        A moment-based method for large-scale generalized eigenvalue problems
        Appl. Num. Anal. Comp. Math. Vol. 1 No. 2 (2004)
        where A is the identity into an eigenpromem for the matrix B.
        The matrix B gets ill-conditioned with its size growing, leading to
        a lack of convergence especially for the original generalized
        eigenproblem used in earlier versions of this benchmark. The exact
        eigenvalues of B are given by
        k = np.arange(1, n+1)
        w_ex = np.sort(16.*np.power(np.cos(0.5*k*np.pi/(n+1)), 4))
        but unused in this benchmark. """

    d0 = np.r_[5, 6 * np.ones(n - 2), 5]
    d1 = -4 * np.ones(n - 1)
    d2 = np.ones(n - 2)
    B = diags([d2, d1, d0, d1, d2], [-2, -1, 0, 1, 2],
                           shape=(n, n))
    return B


def _mikota_pair(n):
    # Mikota pair acts as a nice test since the eigenvalues
    # are the squares of the integers n, n=1,2,...
    x = 1. / np.arange(1, n + 1)
    B = diags([x], [0], shape=(n, n))
    y = - np.arange(n - 1, 0, -1)
    z = np.arange(2 * n - 1, 0, -2)
    A = diags([y, z, y], [-1, 0, 1], shape=(n, n))
    return A.A.astype(float), B.A.astype(float)


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
        ['lobpcg', 'eigsh', 'eigh']
    ]
    param_names = ['n', 'solver']

    def __init__(self):
        self.time_mikota.__func__.params = list(self.params)
        self.time_mikota.__func__.params[0] = [128, 256, 512, 1024, 2048]
        self.time_mikota.__func__.setup = self.setup_mikota

        self.time_sakurai.__func__.params = list(self.params)
        self.time_sakurai.__func__.params[0] = [49, 99]
        self.time_sakurai.__func__.setup = self.setup_sakurai

    def setup_mikota(self, n, solver):
        self.shape = (n, n)
        self.A, self.B = _mikota_pair(n)

        if solver == 'eigh' and n >= 512:
            # skip: slow, and not useful to benchmark
            raise NotImplementedError()

    def setup_sakurai_rev(self, n, solver):
        self.shape = (n, n)
        self.A = _sakurai_rev(n)
        self.A_dense = self.A.A

    def time_mikota(self, n, solver):
        m = 10
        rng = np.random.default_rng(0)
        X =rng.normal(size=(n, m))
        if solver == 'lobpcg':
            LorU, lower = cho_factor(self.A, lower=0, overwrite_a=0)
            M = LinearOperator(self.shape,
                               matvec=partial(_precond, LorU, lower),
                               matmat=partial(_precond, LorU, lower))
            _, _ = lobpcg(self.A, X, self.B, M, tol=1e-4,
                                maxiter=40, largest=False)
        elif solver == 'eigsh':
            _, _ = eigsh(self.A, k=m, which='SA', tol=1e-9, maxiter=5000,
                                   v0=X[:, 0])
        else:
            _, _ = eigh(self.A, self.B, subset_by_index=(0, m - 1))

    def time_sakurai(self, n, solver):
        m = 3
        rng = np.random.default_rng(0)
        X =rng.normal(size=(n, m))
        if solver == 'lobpcg':
            _, _ = lobpcg(self.A, X, tol=1e-9, maxiter=5000, largest=False)
        elif solver == 'eigsh':
            _, _ = eigsh(self.A, k=m, which='SA', tol=1e-9, maxiter=5000,
                                   v0=X[:, 0])
        else:
            _, _ = eigh(self.A_dense, subset_by_index=(0, m - 1))
