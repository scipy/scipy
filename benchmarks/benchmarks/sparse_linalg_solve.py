"""
Check the speed of the conjugate gradient solver.
"""
import numpy as np
from numpy.testing import assert_equal

from .common import Benchmark, safe_import

with safe_import():
    from scipy import linalg, sparse
    from scipy.sparse.linalg import cg, minres, gmres, tfqmr, spsolve
with safe_import():
    from scipy.sparse.linalg import lgmres
with safe_import():
    from scipy.sparse.linalg import gcrotmk


def _create_sparse_poisson1d(n):
    # Make Gilbert Strang's favorite matrix
    # http://www-math.mit.edu/~gs/PIX/cupcakematrix.jpg
    P1d = sparse.diags([[-1]*(n-1), [2]*n, [-1]*(n-1)], [-1, 0, 1])
    assert_equal(P1d.shape, (n, n))
    return P1d


def _create_sparse_poisson2d(n):
    P1d = _create_sparse_poisson1d(n)
    P2d = sparse.kronsum(P1d, P1d)
    assert_equal(P2d.shape, (n*n, n*n))
    return P2d.tocsr()


class Bench(Benchmark):
    params = [
        [4, 6, 10, 16, 25, 40, 64, 100],
        ['dense', 'spsolve', 'cg', 'minres', 'gmres', 'lgmres', 'gcrotmk',
         'tfqmr']
    ]
    mapping = {'spsolve': spsolve, 'cg': cg, 'minres': minres, 'gmres': gmres,
               'lgmres': lgmres, 'gcrotmk': gcrotmk, 'tfqmr': tfqmr}
    param_names = ['(n,n)', 'solver']

    def setup(self, n, solver):
        if solver == 'dense' and n >= 25:
            raise NotImplementedError()

        self.b = np.ones(n*n)
        self.P_sparse = _create_sparse_poisson2d(n)

        if solver == 'dense':
            self.P_dense = self.P_sparse.A

    def time_solve(self, n, solver):
        if solver == 'dense':
            linalg.solve(self.P_dense, self.b)
        else:
            self.mapping[solver](self.P_sparse, self.b)


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
