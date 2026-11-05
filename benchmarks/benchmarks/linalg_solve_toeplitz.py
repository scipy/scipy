"""Benchmark the solve_toeplitz solver (Levinson recursion)
"""
import numpy as np
from .common import Benchmark, safe_import

with safe_import():
    import scipy.linalg


class SolveToeplitz(Benchmark):
    params = (
        ('float64', 'complex128'),
        (100, 300, 1000),
        ('toeplitz', 'generic')
    )
    param_names = ('dtype', 'n', 'solver')

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

        self.c = c
        self.r = r
        self.y = y
        self.T = scipy.linalg.toeplitz(c, r=r)

    def time_solve_toeplitz(self, dtype, n, soltype):
        if soltype == 'toeplitz':
            scipy.linalg.solve_toeplitz((self.c, self.r), self.y)
        else:
            scipy.linalg.solve(self.T, self.y)


class SolveToeplitzMultipleRHS(Benchmark):
    params = (
        ('float64', 'complex128'),
        (300, 1000, 4000),
        (10, 100, 1000),
    )
    param_names = ('dtype', 'n', 'nrhs')

    def setup(self, dtype, n, nrhs):
        random = np.random.RandomState(1234)

        # Positive definite covariance of a Gaussian process on a regular grid
        c = np.exp(-(np.arange(n) / 10.)**2)
        c[0] += 1e-6
        y = random.randn(n, nrhs)
        if dtype == 'complex128':
            y = y + 1j*random.randn(n, nrhs)

        self.c = c
        self.y = y

    def time_solve_toeplitz(self, dtype, n, nrhs):
        scipy.linalg.solve_toeplitz(self.c, self.y)
