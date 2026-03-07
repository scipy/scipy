"""
Check the speed of the conjugate gradient solver.
"""
import inspect

import numpy as np
from numpy.testing import assert_equal

from .common import Benchmark, SkipNotImplemented, XPBenchmark, safe_import

with safe_import():
    from scipy import linalg, sparse
    from scipy.sparse.linalg import cg, minres, gmres, tfqmr, spsolve, LinearOperator
    from scipy._lib._array_api import np_compat
with safe_import():
    from scipy.sparse.linalg import lgmres
with safe_import():
    from scipy.sparse.linalg import gcrotmk


def _create_sparse_poisson1d(n):
    # Make Gilbert Strang's favorite matrix
    # http://www-math.mit.edu/~gs/PIX/cupcakematrix.jpg
    P1d = sparse.diags(
        [[-1]*(n-1), [2]*n, [-1]*(n-1)],
        offsets=[-1, 0, 1],
        dtype=np.float64,
    )
    assert_equal(P1d.shape, (n, n))
    return P1d


def _create_sparse_poisson2d(n):
    P1d = _create_sparse_poisson1d(n)
    P2d = sparse.kronsum(P1d, P1d)
    assert_equal(P2d.shape, (n*n, n*n))
    return sparse.csr_array(P2d)


def _create_sparse_poisson2d_coo(n):
    P1d = _create_sparse_poisson1d(n)
    P2d = sparse.kronsum(P1d, P1d)
    assert_equal(P2d.shape, (n*n, n*n))
    return sparse.coo_array(P2d)


class Bench(Benchmark):
    params = [
        [4, 6, 10, 16, 25, 40, 64, 100],
        # ['dense', 'spsolve', 'cg', 'minres', 'gmres', 'lgmres', 'gcrotmk',
        #  'tfqmr']
        ['cg'],
        [False, True],
    ]
    mapping = {'spsolve': spsolve, 'cg': cg, 'minres': minres, 'gmres': gmres,
               'lgmres': lgmres, 'gcrotmk': gcrotmk, 'tfqmr': tfqmr}
    param_names = ['(n,n)', 'solver', 'pydata_sparse']

    def setup(self, n, solver, pydata_sparse):
        if solver == 'dense' and n >= 25:
            raise NotImplementedError()

        self.b = np.ones(n*n)
        P_sparse = _create_sparse_poisson2d(n)
        if pydata_sparse:
            with safe_import() as sparse_import:
                import sparse as pydata_sparse_module
            if sparse_import.error:
                raise SkipNotImplemented("pydata/sparse not available")
            P_sparse = pydata_sparse_module.GCXS.from_scipy_sparse(P_sparse)
        self.P_sparse = P_sparse

        if solver == 'dense':
            self.P_dense = self.P_sparse.toarray()

    def time_solve(self, n, solver, pydata_sparse):
        if solver == 'dense':
            linalg.solve(self.P_dense, self.b)
        else:
            self.mapping[solver](self.P_sparse, self.b)


class BatchedCGSparse(Benchmark):
    params = [
        # [2, 4, 6, 8, 16, 32, 64],
        # [1, 10, 100, 500, 1000, 5000, 10000],
        [2, 8],
        [1, 10],
        [False, True],
    ]
    param_names = ['(n,n)', 'batch_size', 'pydata_sparse']

    def setup(self, n, batch_size, pydata_sparse):
        if n >= 32 and batch_size >= 500:
            raise NotImplementedError()
        if n >= 16 and batch_size > 5000:
            raise NotImplementedError()
        rng = np.random.default_rng(42)

        self.batched = "xp" in inspect.signature(LinearOperator.__init__).parameters
        if self.batched:
            P_sparse = _create_sparse_poisson2d_coo(n)
            if batch_size > 1:
                P_sparse = sparse.vstack(
                    [P_sparse] * batch_size, format="coo"
                ).reshape(batch_size, n*n, n*n)
                self.b = rng.standard_normal((batch_size, n*n))
            else:
                P_sparse = P_sparse
                self.b = rng.standard_normal(n*n)
        else:
            P_sparse = _create_sparse_poisson2d(n)
            self.b = [rng.standard_normal(n*n) for _ in range(batch_size)]
        
        if pydata_sparse:
            with safe_import() as sparse_import:
                import sparse as pydata_sparse_module
            if sparse_import.error:
                raise SkipNotImplemented("pydata/sparse not available")
            P_sparse = pydata_sparse_module.GCXS.from_coo(
                pydata_sparse_module.COO.from_scipy_sparse(P_sparse)
            )
        
        self.P_sparse = P_sparse

    def time_solve(self, n, batch_size, pydata_sparse):
        if self.batched:
            cg(self.P_sparse, self.b)
        else:
            for i in range(batch_size):
                cg(self.P_sparse, self.b[i])


class BatchedCGDense(XPBenchmark):
    param_names = (*XPBenchmark.param_names, '(n,n)', 'batch_size')
    params = (
        *XPBenchmark.params,
        [2, 4, 8, 16, 24],
        [1, 10, 100, 500, 1000, 5000, 10000],
    )

    def setup(self, backend, n, batch_size):
        super().setup(backend, cg, static_argnames="method")
        self.xp = xp = np_compat if self.xp is np else self.xp

        if n >= 24 and batch_size > 100:
            raise NotImplementedError()
        if n >= 16 and batch_size > 500:
            raise NotImplementedError()
        rng = np.random.default_rng(42)

        self.batched = "xp" in inspect.signature(LinearOperator.__init__).parameters
        if self.batched:
            if batch_size > 1:
                self.A = _create_dense_random(n, batch_shape=(batch_size,), xp=self.xp)
                self.b = xp.asarray(rng.standard_normal((batch_size, n*n)))
            else:
                self.A = _create_dense_random(n, xp=self.xp)
                self.b = xp.asarray(rng.standard_normal(n*n))
        else:
            self.A = _create_dense_random(n)
            self.b = [xp.asarray(rng.standard_normal(n*n)) for _ in range(batch_size)]

        if self.warmup:
            if self.batched:
                cg(self.A, self.b)
            else:
                cg(self.A, self.b[0])

    def time_solve(self, backend, n, batch_size):
        if self.batched:
            cg(self.A, self.b)
        else:
            for i in range(batch_size):
                cg(self.A, self.b[i])


def _create_dense_random(n, batch_shape=None, xp=None):
    xp = np_compat if xp is None else xp
    rng = np.random.default_rng(42)
    M = rng.standard_normal((n*n, n*n))
    M = xp.asarray(M)
    reg = 1e-3
    if batch_shape:
        M = xp.broadcast_to(M[np.newaxis, ...], (*batch_shape, n*n, n*n))
    def matvec(x):
        return xp.squeeze(M.mT @ (M @ x[..., xp.newaxis]), axis=-1) + reg * x

    if "xp" in inspect.signature(LinearOperator.__init__).parameters:
        return LinearOperator(shape=M.shape, matvec=matvec, dtype=xp.float64, xp=xp)
    else:
        return LinearOperator(shape=M.shape, matvec=matvec, dtype=xp.float64)


class Lgmres(Benchmark):
    params = [
        [10, 50, 100, 1000, 10000],
        [10, 30, 60, 90, 180],
    ]
    param_names = ['n', 'm']

    def setup(self, n, m):
        rng = np.random.default_rng(1234)
        self.A = sparse.eye(n, n) + sparse.rand(n, n, density=0.01, random_state=rng)
        self.b = np.ones(n)

    def time_inner(self, n, m):
        lgmres(self.A, self.b, inner_m=m, maxiter=1)
