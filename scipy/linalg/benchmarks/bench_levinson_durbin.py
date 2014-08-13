"""
Benchmark the Levinson-Durbin implementation.

This algorithm solves Toeplitz systems.

"""
from __future__ import division, print_function, absolute_import

import time

import numpy as np
from numpy.testing import assert_allclose
import scipy.linalg


def bench_levinson_durbin():
    np.random.seed(1234)
    print()
    print('                    Levinson-Durbin vs. generic solver')
    print('==============================================================')
    print('      shape      |  solver   |        dtype       |   time   ')
    print('                 |           |                    | (seconds)')
    print('--------------------------------------------------------------')
    fmt = ' %15s |   %6s  | %18s | %6.2f '
    for n in (100, 300, 1000):
        for dtype in (np.float64, np.complex128):

            # Sample a random Toeplitz matrix representation and rhs.
            c = np.random.randn(n)
            r = np.random.randn(n)
            y = np.random.randn(n)
            if dtype == np.complex128:
                c = c + 1j*np.random.rand(n)
                r = r + 1j*np.random.rand(n)
                y = y + 1j*np.random.rand(n)

            # generic solver
            tm = time.clock()
            T = scipy.linalg.toeplitz(c, r=r)
            x_generic = scipy.linalg.solve(T, y)
            nseconds = time.clock() - tm
            print(fmt % (T.shape, 'generic', T.dtype, nseconds))

            # toeplitz-specific solver
            tm = time.clock()
            x_toeplitz = scipy.linalg.solve_toeplitz(c, r=r, y=y)
            nseconds = time.clock() - tm
            print(fmt % (T.shape, 'toeplitz', T.dtype, nseconds))

            # Check that the solutions are the same.
            assert_allclose(x_generic, x_toeplitz)
