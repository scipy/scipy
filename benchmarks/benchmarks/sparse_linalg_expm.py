"""benchmarks for the scipy.sparse.linalg._expm_multiply module"""
import math

import numpy as np
from .common import Benchmark, safe_import

with safe_import():
    import scipy.linalg
    from scipy.sparse.linalg import expm as sp_expm
    from scipy.sparse.linalg import expm_multiply


def random_sparse_csr(m, n, nnz_per_row):
    # Copied from the scipy.sparse benchmark.
    rows = np.arange(m).repeat(nnz_per_row)
    cols = np.random.randint(0, n, size=nnz_per_row*m)
    vals = np.random.random_sample(m*nnz_per_row)
    M = scipy.sparse.coo_matrix((vals, (rows, cols)), (m, n), dtype=float)
    return M.tocsr()


def random_sparse_csc(m, n, nnz_per_row, rng):
    # Copied from the scipy.sparse benchmark.
    rows = np.arange(m).repeat(nnz_per_row)
    cols = rng.integers(0, n, size=nnz_per_row*m)
    vals = rng.random(m*nnz_per_row)
    M = scipy.sparse.coo_matrix((vals, (rows, cols)), (m, n), dtype=float)
    # Use csc instead of csr, because sparse LU decomposition
    # raises a warning when I use csr.
    return M.tocsc()


class ExpmMultiply(Benchmark):
    def setup(self):
        self.n = 2000
        self.i = 100
        self.j = 200
        nnz_per_row = 25
        self.A = random_sparse_csr(self.n, self.n, nnz_per_row)

    def time_expm_multiply(self):
        # computing only column', j, 'of expm of the sparse matrix
        v = np.zeros(self.n, dtype=float)
        v[self.j] = 1
        A_expm_col_j = expm_multiply(self.A, v)
        A_expm_col_j[self.i]


class Expm(Benchmark):
    params = [
        [30, 100, 300],
        ['sparse', 'dense']
    ]
    param_names = ['n', 'format']

    def setup(self, n, format):
        rng = np.random.default_rng(1234)

        # Let the number of nonzero entries per row
        # scale like the log of the order of the matrix.
        nnz_per_row = int(math.ceil(math.log(n)))

        # time the sampling of a random sparse matrix
        self.A_sparse = random_sparse_csc(n, n, nnz_per_row, rng)

        # first format conversion
        self.A_dense = self.A_sparse.toarray()

    def time_expm(self, n, format):
        if format == 'sparse':
            sp_expm(self.A_sparse)
        elif format == 'dense':
            scipy.linalg.expm(self.A_dense)
