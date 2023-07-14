from functools import partial

import numpy as np
from .common import Benchmark, safe_import

with safe_import():
    from scipy.linalg import eigh, cho_factor, cho_solve, sakurai
    from scipy.linalg import cholesky_banded, cho_solve_banded, eig_banded
    from scipy.sparse import diags
    from scipy.sparse.linalg import lobpcg, eigsh, LinearOperator


def _mikota_pair(n):
    """
    Mikota pair acts as a nice test since the eigenvalues
    are the squares of the integers n, n=1,2,... """
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
        self.time_sakurai.__func__.params[0] = [99, 999]
        self.time_sakurai.__func__.setup = self.setup_sakurai

    def setup_mikota(self, n, solver):
        self.shape = (n, n)
        self.A, self.B = _mikota_pair(n)

        if solver == 'eigh' and n >= 512:
            # skip: slow, and not useful to benchmark
            raise NotImplementedError()

    def setup_sakurai(self, n, solver):
        self.shape = (n, n)
        sakurai_obj = sakurai(n)
        self.A = sakurai_obj.banded

    
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
            c = cholesky_banded(self.A)
            a_f = lambda x: cho_solve_banded((c, False), x)
            _, _ = lobpcg(a_f, X, tol=1e-9, maxiter=10)
        elif solver == 'eigsh':
            c = cholesky_banded(self.A)
            a_f = lambda x: cho_solve_banded((c, False), x)
            a_l = LinearOperator((n, n), matvec=a_f, matmat=a_f, dtype='float64')
            _, _ = eigsh(a_l, k=m, which='LA', tol=1e-9, maxiter=10,
                                   v0=X[:, 0])
        else:
            _, _ = eig_banded(self.A, select='i', select_range=[0, m-1])
