"""Compare the speed of expm for sparse vs. dense matrices.
"""
from __future__ import division, print_function, absolute_import

import time
import math

import numpy as np
from numpy.testing import (Tester, TestCase, assert_allclose)

import scipy.sparse
import scipy.linalg


def random_sparse(m, n, nnz_per_row):
    # Copied from the scipy.sparse benchmark.
    rows = np.arange(m).repeat(nnz_per_row)
    cols = np.random.random_integers(low=0, high=n-1, size=nnz_per_row*m)
    vals = np.random.random_sample(m*nnz_per_row)
    M = scipy.sparse.coo_matrix((vals,(rows,cols)), (m,n), dtype=float)
    # Use csc instead of csr, because sparse LU decomposition
    # raises a warning when I use csr.
    return M.tocsc()


class BenchmarkExpm(TestCase):

    def bench_expm_sparse_vs_dense(self):

        # print headers and define the column formats
        print()
        print('               Sparse and Dense Matrix Exponential')
        print('==============================================================')
        print('      shape      |   nnz   |      operation     |   time   ')
        print('                 | per row |                    | (seconds)')
        print('--------------------------------------------------------------')
        fmt = ' %15s |   %3d   | %18s | %6.2f '

        # check three roughly exponentially increasing matrix orders
        np.random.seed(1234)
        for n in (30, 100, 300):

            # Let the number of nonzero entries per row
            # scale like the log of the order of the matrix.
            nnz_per_row = int(math.ceil(math.log(n)))
            shape = (n, n)

            # time the sampling of a random sparse matrix
            tm_start = time.clock()
            A_sparse = random_sparse(n, n, nnz_per_row)
            tm_end = time.clock()
            tm_sampling = tm_end - tm_start

            # first format conversion
            tm_start = time.clock()
            A_dense = A_sparse.toarray()
            tm_end = time.clock()
            tm_first_fmt = tm_end - tm_start

            # sparse matrix exponential
            tm_start = time.clock()
            A_sparse_expm = scipy.linalg.expm(A_sparse)
            tm_end = time.clock()
            tm_sparse = tm_end - tm_start

            # dense matrix exponential
            tm_start = time.clock()
            A_dense_expm = scipy.linalg.expm(A_dense)
            tm_end = time.clock()
            tm_dense = tm_end - tm_start

            # second format conversion
            tm_start = time.clock()
            A_sparse_expm_as_dense = A_sparse_expm.toarray()
            tm_end = time.clock()
            tm_second_fmt = tm_end - tm_start

            # sum the format conversion times
            tm_fmt = tm_first_fmt + tm_second_fmt

            # check that the matrix exponentials are the same
            assert_allclose(A_sparse_expm_as_dense, A_dense_expm)

            # write the rows
            print(fmt % (shape, nnz_per_row, 'matrix sampling', tm_sampling))
            print(fmt % (shape, nnz_per_row, 'fmt conversions', tm_fmt))
            print(fmt % (shape, nnz_per_row, 'sparse expm', tm_sparse))
            print(fmt % (shape, nnz_per_row, 'dense expm', tm_dense))
        print()


if __name__ == '__main__':
    Tester().bench()
