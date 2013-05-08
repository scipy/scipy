"""benchmarks for the scipy.sparse.linalg._expm_multiply module"""
from __future__ import division, print_function, absolute_import

import time

import numpy as np
from numpy.testing import (assert_allclose, run_module_suite)

import scipy.linalg
from scipy.sparse.linalg import expm_multiply


class BenchmarkExpmMultiply(TestCase):

    def _help_bench_expm_multiply(self, A, i, j):
        print('computing full expm of A...')
        tm_start = time.clock()
        A_expm = scipy.linalg.expm(A)
        full_expm_entry = A_expm[i, j]
        tm_end = time.clock()
        print(tm_end - tm_start, ' seconds')
        print('expm(A)[%d, %d]:' % (i, j), full_expm_entry)
        print()
        print('computing only column', j, 'of expm of A...')
        tm_start = time.clock()
        v = np.zeros(n, dtype=float)
        v[j] = 1
        A_expm_col_j = expm_multiply(A, v)
        expm_col_entry = A_expm_col_j[i]
        tm_end = time.clock()
        print(tm_end - tm_start, ' seconds')
        print('expm(A)[%d, %d]:' % (i, j), expm_col_entry)
        print()
        if np.allclose(full_expm_entry, expm_col_entry):
            print('The two methods give the same answer.')
        else:
            print('!!! The two methods give different answers. !!!')
        print()

    def bench_expm_multiply(self):
        np.random.seed(1234)
        n = 1000
        i = 100
        j = 200
        print()
        print('expm multiply benchmarking')
        print('--------------------------')
        print()
        print('sampling a matrix')
        print('with %dx%d random standard normal entries...' % (n, n))
        tm_start = time.clock()
        A = np.random.randn(n, n)
        tm_end = time.clock()
        print(tm_end - tm_start, ' seconds')
        print('A.shape:', A.shape)
        print()
        self._help_bench_expm_multiply(A)
        print()
        print('sampling a matrix')
        print('with %dx%d random standard uniform entries...' % (n, n))
        tm_start = time.clock()
        A = np.random.uniform(size=(n, n))
        tm_end = time.clock()
        print(tm_end - tm_start, ' seconds')
        print('A.shape:', A.shape)
        print()
        self._help_bench_expm_multiply(A)
        print()


if __name__ == '__main__':
    run_module_suite()

