"""Benchmark the solve_toeplitz solver (Levinson recursion)
"""
from __future__ import division, print_function, absolute_import

import time

import numpy as np
from numpy.testing import assert_array_almost_equal, Tester
import scipy.linalg


def bench_solve_toeplitz():
    random = np.random.RandomState(1234)
    print()
    print('                    solve_toeplitz vs. generic solver')
    print('==============================================================')
    print('      shape      |  solver     |        dtype       |   time   ')
    print('                 |             |                    | (seconds)')
    print('--------------------------------------------------------------')
    fmt = ' %15s |   %s  | %18s | %6.2f '
    for n in (100, 300, 1000):
        for dtype in (np.float64, np.complex128):

            # Sample a random Toeplitz matrix representation and rhs.
            c = random.randn(n)
            r = random.randn(n)
            y = random.randn(n)
            if dtype == np.complex128:
                c = c + 1j*random.rand(n)
                r = r + 1j*random.rand(n)
                y = y + 1j*random.rand(n)

            # generic solver
            tm = time.clock()
            T = scipy.linalg.toeplitz(c, r=r)
            x_generic = scipy.linalg.solve(T, y)
            nseconds = time.clock() - tm
            print(fmt % (T.shape, 'generic ', T.dtype, nseconds))

            # toeplitz-specific solver
            tm = time.clock()
            x_toeplitz = scipy.linalg.solve_toeplitz((c, r), y)
            nseconds = time.clock() - tm
            print(fmt % (T.shape, 'toeplitz', T.dtype, nseconds))

            # Check that the solutions are the sameself.
            assert_array_almost_equal(x_generic, x_toeplitz)


if __name__ == "__main__":
    Tester().bench()
