"""
Check the speed of the conjugate gradient solver.
"""
from __future__ import division, absolute_import, print_function

import time

import numpy as np
from numpy.testing import assert_allclose, assert_equal

try:
    from scipy import linalg, sparse
    import scipy.sparse.linalg
except ImportError:
    pass

from .common import Benchmark


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
    return P2d


class Bench(Benchmark):
    params = [
        [4, 6, 10, 16, 25, 40, 64, 100, 160, 250, 400, 640, 1000, 1600],
        ['sparse', 'dense']
    ]
    param_names = ['(n,n)', 'solver']

    def setup(self, n, solver):
        dense_is_active = (n**2 < 600)
        sparse_is_active = (n**2 < 20000)

        if solver == 'dense' and not dense_is_active:
            raise NotImplementedError()

        if solver == 'sparse' and not sparse_is_active:
            raise NotImplementedError()

        self.b = np.ones(n*n)
        self.P_sparse = _create_sparse_poisson2d(n)
        self.P_dense = self.P_sparse.A

    def time_cg(self, n, solver):
        if solver == 'dense':
            linalg.solve(self.P_dense, self.b)
        else:
            sparse.linalg.cg(self.P_sparse, self.b)

    def time_spsolve(self, n, solver):
        if solver == 'dense':
            linalg.solve(self.P_dense, self.b)
        else:
            sparse.linalg.cg(self.P_sparse, self.b)


if __name__ == '__main__':
    Tester().bench()
