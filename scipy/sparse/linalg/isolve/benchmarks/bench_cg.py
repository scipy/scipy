"""Check the speed of the conjugate gradient solver.
"""
from __future__ import division, print_function, absolute_import

import time

import numpy as np
from numpy.testing import Tester, TestCase, assert_allclose, assert_equal

import scipy.linalg
import scipy.sparse

def _create_sparse_poisson1d(n):
    # Make Gilbert Strang's favorite matrix
    # http://www-math.mit.edu/~gs/PIX/cupcakematrix.jpg
    P1d = scipy.sparse.diags([[-1]*(n-1), [2]*n, [-1]*(n-1)], [-1, 0, 1])
    assert_equal(P1d.shape, (n, n))
    return P1d

def _create_sparse_poisson2d(n):
    P1d = _create_sparse_poisson1d(n)
    P2d = scipy.sparse.kronsum(P1d, P1d)
    assert_equal(P2d.shape, (n*n, n*n))
    return P2d

class BenchmarkConjuateGradientSolver(TestCase):

    def bench_cg(self):

        # print headers and define the column formats
        print()
        print('         generic solve vs. conjugate gradient solve')
        print('==============================================================')
        print('       shape       | repeats |      operation     |   time   ')
        print('                                                  | (seconds)')
        print('--------------------------------------------------------------')
        fmt = ' %17s |   %3d   | %18s | %6.2f '

        dense_is_active = True
        sparse_is_active = True
        repeats = 100
        for n in 4, 6, 10, 16, 25, 40, 64, 100, 160, 250, 400, 640, 1000, 1600:
            if not dense_is_active and not sparse_is_active:
                break
            b = np.ones(n*n)
            P_sparse = _create_sparse_poisson2d(n)

            # Optionally use the generic dense solver.
            if dense_is_active:
                P_dense = P_sparse.A
                tm_start = time.clock()
                for i in range(repeats):
                    x_dense = scipy.linalg.solve(P_dense, b)
                tm_end = time.clock()
                tm_dense = tm_end - tm_start

            # Optionally use the sparse conjugate gradient solver.
            if sparse_is_active:
                tm_start = time.clock()
                for i in range(repeats):
                    x_sparse, info = scipy.sparse.linalg.cg(P_sparse, b)
                tm_end = time.clock()
                tm_sparse = tm_end - tm_start

            # Check that the solutions are close to each other.
            if dense_is_active and sparse_is_active:
                assert_allclose(x_dense, x_sparse, rtol=1e-4)

            # Write the rows.
            shape = (n*n, n*n)
            if dense_is_active:
                print(fmt % (shape, repeats, 'dense solve', tm_dense))
            if sparse_is_active:
                print(fmt % (shape, repeats, 'sparse cg', tm_sparse))

            dense_is_active = (tm_dense < 5)
            sparse_is_active = (tm_sparse < 5)

        print()


if __name__ == '__main__':
    Tester().bench()
