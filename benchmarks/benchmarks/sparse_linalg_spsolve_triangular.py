"""
Check the speed of the sparse triangular solve function.
"""
import numpy as np
from numpy.testing import assert_equal

from .common import Benchmark, safe_import

with safe_import():
    from scipy import sparse
    from scipy.sparse.linalg import spsolve, spsolve_triangular

def _create_sparse_poisson1d(n):
    # Make Gilbert Strang's favorite matrix
    # http://www-math.mit.edu/~gs/PIX/cupcakematrix.jpg
    # and take the lower triangular half
    P1d = sparse.diags([[-1]*(n-1), [2]*n, [-1]*(n-1)], [-1, 0, 1])
    assert_equal(P1d.shape, (n, n))
    return P1d


def _create_sparse_poisson2d_half(n):
    P1d = _create_sparse_poisson1d(n)
    P2d = sparse.kronsum(P1d, P1d)
    assert_equal(P2d.shape, (n*n, n*n))
    return sparse.tril(P2d).tocsr()


class Bench(Benchmark):
    params = [
        [100,1000],
        ["spsolve", "spsolve_triangular"],
    ]
    param_names = ['(n,n)',"method"]

    def setup(self, n, method):
        self.b = np.ones(n*n)
        self.P_sparse = _create_sparse_poisson2d_half(n)

    def time_solve(self, n, method):
        if method == "spsolve":
            spsolve(self.P_sparse, self.b)
        elif method == "spsolve_triangular":
            spsolve_triangular(self.P_sparse, self.b)
        else:
            raise NotImplementedError()

