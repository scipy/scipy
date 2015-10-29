"""Benchmark the solve_toeplitz solver (Levinson recursion)
"""
from __future__ import division, absolute_import, print_function

import numpy as np

try:
    import scipy.linalg
except ImportError:
    pass

from .common import Benchmark


class SolveToeplitz(Benchmark):
    params = (
        ('float64', 'complex128'),
        (100, 300, 1000),
        ('toeplitz', 'generic')
    )
    param_names = ('dtype', 'n', 'solver')

    def __init__(self, T, y, r, c):
        self.T = scipy.linalg.toeplitz(c, r=r)
        self.y = y
        self.r = r
        self.c = c

    def setup(self, dtype, n, soltype):
        random = np.random.RandomState(1234)

        dtype = np.dtype(dtype)

        # Sample a random Toeplitz matrix representation and rhs.
        c = random.randn(n)
        r = random.randn(n)
        y = random.randn(n)
        if dtype == np.complex128:
            c = c + 1j*random.rand(n)
            r = r + 1j*random.rand(n)
            y = y + 1j*random.rand(n)

    def time_solve_toeplitz(self, dtype, n, soltype):
        if soltype == 'toeplitz':
            scipy.linalg.solve_toeplitz((self.c, self.r), self.y)
        else:
            scipy.linalg.solve(self.T, self.y)
