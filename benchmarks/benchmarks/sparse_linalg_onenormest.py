"""Compare the speed of exact one-norm calculation vs. its estimation.
"""
from __future__ import division, print_function, absolute_import

import time

import numpy as np

import scipy.sparse.linalg


class BenchmarkOneNormEst(object):
    params = [
        [2, 3, 5, 10, 30, 100, 300, 500, 1000],
        ['exact', 'onenormest']
    ]
    param_names = ['n', 'solver']
    goal_time = 0.5

    def setup(self, n, solver):
        np.random.seed(1234)
        nrepeats = 100
        shape = (n, n)

        # Sample the matrices.
        self.matrices = []
        for i in range(nrepeats):
            M = np.random.randn(*shape)
            self.matrices.append(M)

    def time_onenormest(self, n, solver):
        if solver == 'exact':
            # Get the exact values of one-norms of squares.
            for M in self.matrices:
                M2 = M.dot(M)
                scipy.sparse.linalg.matfuncs._onenorm(M)
        elif solver == 'onenormest':
            # Get the estimates of one-norms of squares.
            for M in self.matrices:
                scipy.sparse.linalg.matfuncs._onenormest_matrix_power(M, 2)
