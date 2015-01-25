"""benchmarks for the scipy.sparse.linalg._expm_multiply module"""
from __future__ import division, print_function, absolute_import

import time

import numpy as np
from numpy.testing import Tester, TestCase

import scipy.linalg
from scipy.sparse.linalg import expm_multiply


def random_sparse(m, n, nnz_per_row):
    # Copied from the scipy.sparse benchmark.
    rows = np.arange(m).repeat(nnz_per_row)
    cols = np.random.random_integers(low=0, high=n-1, size=nnz_per_row*m)
    vals = np.random.random_sample(m*nnz_per_row)
    M = scipy.sparse.coo_matrix((vals,(rows,cols)), (m,n), dtype=float)
    return M.tocsr()


class BenchmarkExpmMultiply(TestCase):

    def _help_bench_expm_multiply(self, A, i, j):
        n = A.shape[0]
        print('converting the sparse matrix to a dense array...')
        tm_start = time.clock()
        A_dense = A.toarray()
        tm_end = time.clock()
        print(tm_end - tm_start, ' seconds')
        print()
        print('computing full expm of the dense array...')
        tm_start = time.clock()
        A_expm = scipy.linalg.expm(A_dense)
        full_expm_entry = A_expm[i, j]
        tm_end = time.clock()
        print('expm(A)[%d, %d]:' % (i, j), full_expm_entry)
        print(tm_end - tm_start, ' seconds')
        print()
        print('computing only column', j, 'of expm of the sparse matrix...')
        tm_start = time.clock()
        v = np.zeros(n, dtype=float)
        v[j] = 1
        A_expm_col_j = expm_multiply(A, v)
        expm_col_entry = A_expm_col_j[i]
        tm_end = time.clock()
        print('expm(A)[%d, %d]:' % (i, j), expm_col_entry)
        print(tm_end - tm_start, ' seconds')
        print()
        if np.allclose(full_expm_entry, expm_col_entry):
            print('The two methods give the same answer.')
        else:
            print('!!! The two methods give different answers. !!!')
        print()

    def bench_expm_multiply(self):
        np.random.seed(1234)
        n = 2000
        i = 100
        j = 200
        shape = (n, n)
        nnz_per_row = 25
        print()
        print('expm multiply benchmarking')
        print('--------------------------')
        print()
        print('sampling a random sparse matrix...')
        print('shape:', shape)
        print('nnz per row:', nnz_per_row)
        tm_start = time.clock()
        A = random_sparse(n, n, nnz_per_row)
        tm_end = time.clock()
        print(tm_end - tm_start, ' seconds')
        print()
        self._help_bench_expm_multiply(A, i, j)
        print()


if __name__ == '__main__':
    Tester().bench()
