"""Compare the speed of exact one-norm calculation vs. its estimation.
"""
from __future__ import division, print_function, absolute_import

import time

import numpy as np
from numpy.testing import (Tester, TestCase, assert_allclose)

import scipy.sparse


class BenchmarkOneNormEst(TestCase):

    def bench_onenormest(self):

        # print headers and define the column formats
        print()
        print('   calculation and estimation of one-norm of matrix squaring')
        print('==============================================================')
        print('      shape      | repeats |      operation     |   time   ')
        print('                                                | (seconds)')
        print('--------------------------------------------------------------')
        fmt = ' %15s |   %3d   | %18s | %6.2f '

        np.random.seed(1234)
        nrepeats = 100
        for n in (2, 3, 5, 10, 30, 100, 300, 500, 1000):
            shape = (n, n)

            # Sample the matrices.
            tm_start = time.clock()
            matrices = []
            for i in range(nrepeats):
                M = np.random.randn(*shape)
                matrices.append(M)
            tm_end = time.clock()
            tm_sampling = tm_end - tm_start

            # Get the exact values of one-norms of squares.
            tm_start = time.clock()
            for M in matrices:
                M2 = M.dot(M)
                scipy.sparse.linalg.matfuncs._onenorm(M)
            tm_end = time.clock()
            tm_exact = tm_end - tm_start

            # Get the estimates of one-norms of squares.
            tm_start = time.clock()
            for M in matrices:
                scipy.sparse.linalg.matfuncs._onenormest_matrix_power(M, 2)
            tm_end = time.clock()
            tm_estimate = tm_end - tm_start

            # write the rows
            print(fmt % (shape, nrepeats, 'matrix sampling', tm_sampling))
            print(fmt % (shape, nrepeats, 'one-norm exact', tm_exact))
            print(fmt % (shape, nrepeats, 'one-norm estimate', tm_estimate))
        print()


if __name__ == '__main__':
    Tester().bench()
