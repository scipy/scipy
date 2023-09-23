import numpy as np
import warnings
from .common import Benchmark, safe_import

with safe_import():
    from scipy.linalg import (eigh,
                              cholesky_banded, cho_solve_banded, eig_banded)
    from scipy.sparse.linalg import lobpcg, eigsh, LinearOperator
    from scipy.sparse.linalg._special_sparse_arrays import (Sakurai,
                                                            MikotaPair)

# ensure that we are benchmarking a consistent outcome;
# (e.g. if the code wasn't able to find a solution at all,
# the timing of the benchmark would become useless).
msg = ("the benchmark code did not converge as expected, "
       "the timing is therefore useless")


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
        self.time_sakurai.__func__.params[0] = [50]
        self.time_sakurai.__func__.setup = self.setup_sakurai

        self.time_sakurai_inverse.__func__.params = list(self.params)
        self.time_sakurai_inverse.__func__.params[0] = [500, 1000]
        self.time_sakurai_inverse.__func__.setup = self.setup_sakurai_inverse


    def setup_mikota(self, n, solver):
        self.shape = (n, n)
        mik = MikotaPair(n)
        mik_k = mik.k
        mik_m = mik.m
        self.Ac = mik_k
        self.Aa = mik_k.toarray()
        self.Bc = mik_m
        self.Ba = mik_m.toarray()
        self.Ab = mik_k.tobanded()
        self.eigenvalues = mik.eigenvalues

        if solver == 'lapack' and n > 512:
            # skip: slow, and not useful to benchmark
            raise NotImplementedError()

    def setup_sakurai(self, n, solver):
        self.shape = (n, n)
        sakurai_obj = Sakurai(n, dtype='int')
        self.A = sakurai_obj
        self.Aa = sakurai_obj.toarray()
        self.eigenvalues = sakurai_obj.eigenvalues

    def setup_sakurai_inverse(self, n, solver):
        self.shape = (n, n)
        sakurai_obj = Sakurai(n)
        self.A = sakurai_obj.tobanded().astype(np.float64)
        self.eigenvalues = sakurai_obj.eigenvalues
    
    def time_mikota(self, n, solver):
        def a(x):
            return cho_solve_banded((c, False), x)

        m = 10
        ee = self.eigenvalues(m)
        tol = m * n * n * n* np.finfo(float).eps
        rng = np.random.default_rng(0)
        X = rng.normal(size=(n, m))
        if solver == 'lobpcg':
            # `lobpcg` allows callable parameters `Ac` and `Bc` directly
            # `lobpcg` solves ``Ax = lambda Bx`` and applies here a preconditioner
            # given by the matrix inverse in `np.float32` of 'Ab` that itself
            # is `np.float64`. 
            c = cholesky_banded(self.Ab.astype(np.float32))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                el, _ = lobpcg(self.Ac, X, self.Bc, M=a, tol=1e-4,
                               maxiter=40, largest=False)
            accuracy = max(abs(ee - el) / ee)
            assert accuracy < tol, msg
        elif solver == 'eigsh':
            # `eigsh` ARPACK here is called on ``Bx = 1/lambda Ax``
            # to get fast convergence speed similar to that of `lobpcg` above
            # requiring the inverse of the matrix ``A`` given by Cholesky on
            # the banded form `Ab` of ``A`` in full `np.float64` precision.
            # `eigsh` ARPACK does not allow the callable parameter `Bc` directly
            # requiring `LinearOperator` format for input in contrast to `lobpcg`
            B = LinearOperator((n, n), matvec=self.Bc, matmat=self.Bc, dtype='float64')
            A = LinearOperator((n, n), matvec=self.Ac, matmat=self.Ac, dtype='float64')
            c = cholesky_banded(self.Ab)
            a_l = LinearOperator((n, n), matvec=a, matmat=a, dtype='float64')
            ea, _ = eigsh(B, k=m, M=A, Minv=a_l, which='LA', tol=1e-4, maxiter=50,
                          v0 = rng.normal(size=(n, 1)))
            accuracy = max(abs(ee - np.sort(1./ea)) / ee)
            assert accuracy < tol, msg
        else:
            # `eigh` is the only dense eigensolver for generalized eigenproblems
            # ``Ax = lambda Bx`` and needs both matrices as dense arrays
            # making it very slow for large matrix sizes
            ed, _ = eigh(self.Aa, self.Ba, subset_by_index=(0, m - 1))
            accuracy = max(abs(ee - ed) / ee)
            assert accuracy < tol, msg

    def time_sakurai(self, n, solver):
        # the Sakurai matrix ``A`` is ill-conditioned so convergence of both
        # `lobpcg` and `eigsh` ARPACK on the matrix ``A`` itself is very slow
        # computing its smallest eigenvalues even from moderate sizes
        # requiring enormous numbers of iterations
        m = 3
        ee = self.eigenvalues(m)
        tol = 100 * n * n * n* np.finfo(float).eps
        rng = np.random.default_rng(0)
        X = rng.normal(size=(n, m))
        if solver == 'lobpcg':
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                el, _ = lobpcg(self.A, X, tol=1e-9, maxiter=5000, largest=False)
            accuracy = max(abs(ee - el) / ee)
            assert accuracy < tol, msg
        elif solver == 'eigsh':
            a_l = LinearOperator((n, n), matvec=self.A, matmat=self.A, dtype='float64')
            ea, _ = eigsh(a_l, k=m, which='SA', tol=1e-9, maxiter=15000,
                          v0 = rng.normal(size=(n, 1)))
            accuracy = max(abs(ee - ea) / ee)
            assert accuracy < tol, msg
        else:
            ed, _ = eigh(self.Aa, subset_by_index=(0, m - 1))
            accuracy = max(abs(ee - ed) / ee)
            assert accuracy < tol, msg

    def time_sakurai_inverse(self, n, solver):
        # apply inverse iterations in  `lobpcg` and `eigsh` ARPACK
        # using the Cholesky on the banded form in full `np.float64` precision
        # for fast convergence and compare to dense banded eigensolver `eig_banded`
        def a(x):
            return cho_solve_banded((c, False), x)
        m = 3
        ee = self.eigenvalues(m)
        tol = 10 * n * n * n* np.finfo(float).eps
        rng = np.random.default_rng(0)
        X = rng.normal(size=(n, m))
        if solver == 'lobpcg':
            c = cholesky_banded(self.A)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                el, _ = lobpcg(a, X, tol=1e-9, maxiter=8)
            accuracy = max(abs(ee - 1. / el) / ee)
            assert accuracy < tol, msg
        elif solver == 'eigsh':
            c = cholesky_banded(self.A)
            a_l = LinearOperator((n, n), matvec=a, matmat=a, dtype='float64')
            ea, _ = eigsh(a_l, k=m, which='LA', tol=1e-9, maxiter=8,
                          v0 = rng.normal(size=(n, 1)))
            accuracy = max(abs(ee - np.sort(1. / ea)) / ee)
            assert accuracy < tol, msg
        else:
            ed, _ = eig_banded(self.A, select='i', select_range=[0, m-1])
            accuracy = max(abs(ee - ed) / ee)
            assert accuracy < tol, msg
