"""
Check the speed of the conjugate gradient solver.
"""
from __future__ import division, absolute_import, print_function

import numpy as np
from numpy.testing import assert_equal

try:
    from scipy import linalg, sparse
    from scipy.sparse.linalg import cg, minres, spsolve
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
    return P2d.tocsr()


class Bench(Benchmark):
    params = [
        [4, 6, 10, 16, 25, 40, 64, 100],
        ['dense', 'spsolve', 'cg', 'minres']
    ]
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
        elif solver == 'cg':
            cg(self.P_sparse, self.b)
        elif solver == 'minres':
            minres(self.P_sparse, self.b)
        elif solver == 'spsolve':
            spsolve(self.P_sparse, self.b)
        else:
            raise ValueError('Unknown solver: %r' % solver)
