"""benchmarks for the scipy.sparse.linalg._expm_multiply module"""
from __future__ import division, print_function, absolute_import

import time
import math

import numpy as np
from numpy.testing import assert_allclose

try:
    import scipy.linalg
    from scipy.sparse.linalg import expm_multiply
except ImportError:
    pass

from .common import Benchmark


def random_sparse_csr(m, n, nnz_per_row):
    # Copied from the scipy.sparse benchmark.
    rows = np.arange(m).repeat(nnz_per_row)
    cols = np.random.random_integers(low=0, high=n-1, size=nnz_per_row*m)
    vals = np.random.random_sample(m*nnz_per_row)
    M = scipy.sparse.coo_matrix((vals,(rows,cols)), (m,n), dtype=float)
    return M.tocsr()


def random_sparse_csc(m, n, nnz_per_row):
    # Copied from the scipy.sparse benchmark.
    rows = np.arange(m).repeat(nnz_per_row)
    cols = np.random.random_integers(low=0, high=n-1, size=nnz_per_row*m)
    vals = np.random.random_sample(m*nnz_per_row)
    M = scipy.sparse.coo_matrix((vals,(rows,cols)), (m,n), dtype=float)
    # Use csc instead of csr, because sparse LU decomposition
    # raises a warning when I use csr.
    return M.tocsc()


class ExpmMultiply(Benchmark):
    params = [['sparse', 'full']]
    param_names = ['run format']

    def setup(self, *args):
        self.n = 2000
        self.i = 100
        self.j = 200
        shape = (self.n, self.n)
        nnz_per_row = 25
        self.A = random_sparse_csr(self.n, self.n, nnz_per_row)
        self.A_dense = self.A.toarray()

    def time_expm_multiply(self, format):
        if format == 'full':
            # computing full expm of the dense array...
            A_expm = scipy.linalg.expm(self.A_dense)
            full_expm_entry = A_expm[self.i, self.j]
        else:
            # computing only column', j, 'of expm of the sparse matrix...
            v = np.zeros(self.n, dtype=float)
            v[self.j] = 1
            A_expm_col_j = expm_multiply(self.A, v)
            expm_col_entry = A_expm_col_j[self.i]


class Expm(Benchmark):
    params = [
        [30, 100, 300],
        ['sparse', 'dense']
    ]
    param_names = ['n', 'format']
    goal_time = 0.5

    def setup(self, n, format):
        np.random.seed(1234)

        # Let the number of nonzero entries per row
        # scale like the log of the order of the matrix.
        nnz_per_row = int(math.ceil(math.log(n)))
        shape = (n, n)

        # time the sampling of a random sparse matrix
        self.A_sparse = random_sparse_csc(n, n, nnz_per_row)

        # first format conversion
        self.A_dense = self.A_sparse.toarray()

    def time_expm(self, n, format):
        if format == 'sparse':
            A_sparse_expm = scipy.linalg.expm(self.A_sparse)
        elif format == 'dense':
            A_dense_expm = scipy.linalg.expm(self.A_dense)
