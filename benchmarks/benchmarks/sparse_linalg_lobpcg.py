import numpy as np
from .common import Benchmark, safe_import

with safe_import():
    from scipy.linalg import (eigh, sakurai, mikota_pair,
                              cholesky_banded, cho_solve_banded, eig_banded)
    from scipy.sparse.linalg import lobpcg, eigsh, LinearOperator


class Bench(Benchmark):
    params = [
        [],
        ['lobpcg', 'eigsh', 'lapack']
    ]
    param_names = ['n', 'solver']

    def __init__(self):
        self.time_mikota.__func__.params = list(self.params)
        self.time_mikota.__func__.params[0] = [128, 256, 512, 1024, 2048]
        self.time_mikota.__func__.setup = self.setup_mikota

        self.time_sakurai.__func__.params = list(self.params)
        self.time_sakurai.__func__.params[0] = [50, 100]
        self.time_sakurai.__func__.setup = self.setup_sakurai

        self.time_sakuraii.__func__.params = list(self.params)
        self.time_sakuraii.__func__.params[0] = [500, 1000, 2000]
        self.time_sakuraii.__func__.setup = self.setup_sakuraii


    def setup_mikota(self, n, solver):
        self.shape = (n, n)
        mikota_pair_obj = mikota_pair(n)
        self.Ac = mikota_pair_obj.Kcallable
        self.Aa = mikota_pair_obj.Karray
        self.Bc = mikota_pair_obj.Mcallable
        self.Ba = mikota_pair_obj.Marray
        self.Ab = mikota_pair_obj.Kbanded
        self.eigenvalues = mikota_pair.eigenvalues

        # if solver == 'eigh' and n >= 512:
        #     # skip: slow, and not useful to benchmark
        #     raise NotImplementedError()


    def setup_sakurai(self, n, solver):
        self.shape = (n, n)
        sakurai_obj = sakurai(n)
        self.A = sakurai_obj.callable
        self.Aa = sakurai_obj.array
        self.eigenvalues = sakurai_obj.eigenvalues


    def setup_sakuraii(self, n, solver):
        self.shape = (n, n)
        sakurai_obj = sakurai(n)
        self.A = sakurai_obj.banded
        self.eigenvalues = sakurai_obj.eigenvalues

    
    def time_mikota(self, n, solver):
        def a(x):
            return cho_solve_banded((c, False), x)
        m = 10
        ee = self.eigenvalues[:m]
        tol = m * n * n * n* np.finfo(float).eps
        rng = np.random.default_rng(0)
        X =rng.normal(size=(n, m))
        if solver == 'lobpcg':
            c = cholesky_banded(self.Ab.astype(np.float32))
            el, _ = lobpcg(self.Ac, X, self.Bc, M=a, tol=1e-4,
                                maxiter=40, largest=False)
            accuracy = max(abs(ee - el) / ee)
            assert accuracy < tol
        elif solver == 'eigsh':
            # B = LinearOperator((n, n), matvec=self.Bc, matmat=self.Bc, dtype='float64')
            A = LinearOperator((n, n), matvec=self.Ac, matmat=self.Ac, dtype='float64')
            c = cholesky_banded(self.Ab)
            a_l = LinearOperator((n, n), matvec=a, matmat=a, dtype='float64')
            ea, _ = eigsh(self.Ba, k=m, M=A, Minv=a_l, which='LA', tol=1e-4, maxiter=50,
                                  v0=rng.normal(size=(n, 1)))
            accuracy = max(abs(ee - np.sort(1./ea)) / ee)
            assert accuracy < tol
        else:
            ed, _ = eigh(self.Aa, self.Ba, subset_by_index=(0, m - 1))
            accuracy = max(abs(ee - ed) / ee)
            assert accuracy < tol


    def time_sakurai(self, n, solver):
        m = 3
        ee = self.eigenvalues[:m]
        tol = 100 * n * n * n* np.finfo(float).eps
        rng = np.random.default_rng(0)
        X =rng.normal(size=(n, m))
        if solver == 'lobpcg':
            el, _ = lobpcg(self.A, X, tol=1e-9, maxiter=5000, largest=False)
            accuracy = max(abs(ee - el) / ee)
            assert accuracy < tol
        elif solver == 'eigsh':
            a_l = LinearOperator((n, n), matvec=self.A, matmat=self.A, dtype='float64')
            ea, _ = eigsh(a_l, k=m, which='SA', tol=1e-9, maxiter=15000,
                                   v0=rng.normal(size=(n, 1)))
            accuracy = max(abs(ee - ea) / ee)
            assert accuracy < tol
        else:
            ed, _ = eigh(self.Aa, subset_by_index=(0, m - 1))
            accuracy = max(abs(ee - ed) / ee)
            assert accuracy < tol


    def time_sakuraii(self, n, solver):
        def a(x):
            return cho_solve_banded((c, False), x)
        m = 3
        ee = self.eigenvalues[:m]
        tol = 10 * n * n * n* np.finfo(float).eps
        rng = np.random.default_rng(0)
        X =rng.normal(size=(n, m))
        if solver == 'lobpcg':
            c = cholesky_banded(self.A)
            el, _ = lobpcg(a, X, tol=1e-9, maxiter=8)
            accuracy = max(abs(ee - 1. / el) / ee)
            assert accuracy < tol
        elif solver == 'eigsh':
            c = cholesky_banded(self.A)
            a_l = LinearOperator((n, n), matvec=a, matmat=a, dtype='float64')
            ea, _ = eigsh(a_l, k=m, which='LA', tol=1e-9, maxiter=8,
                                   v0=rng.normal(size=(n, 1)))
            accuracy = max(abs(ee - np.sort(1./ea)) / ee)
            assert accuracy < tol
        else:
            ed, _ = eig_banded(self.A, select='i', select_range=[0, m-1])
            accuracy = max(abs(ee - ed) / ee)
            assert accuracy < tol
