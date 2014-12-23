"""Check the speed of the conjugate gradient solver.
"""
from __future__ import division, print_function, absolute_import

import time

import numpy as np
from numpy.testing import Tester, TestCase

import scipy.linalg
import scipy.sparse

def _create_sparse_poisson1d(n):
    # Make Gilbert Strang's favorite matrix
    # http://www-math.mit.edu/~gs/PIX/cupcakematrix.jpg
    P1d = scipy.sparse.diags([[-1]*(n-1), [2]*n, [-1]*(n-1)], [-1, 0, 1])
    assert_equal(P1d.shape, (n, n))
    return P1d

def _create_sparse_poisson2d(n):
    P1d = create_sparse_poisson1d(n)
    P2d = scipy.sparse.kronsum(P, P)
    assert_equal(P2d.shape, (n*n, n*n))
    return P2d

class BenchmarkConjuateGradientSolver(TestCase):

    def bench_cg(self):

        # plain solve vs sparse cg vs dense cg
        # solver (dense generic vs cg)
        # shape
        # time

        # print headers and define the column formats
        print()
        print('         generic solve vs. conjugate gradient solve')
        print('==============================================================')
        print('      shape      | repeats |      operation     |   time   ')
        print('                                                | (seconds)')
        print('--------------------------------------------------------------')
        fmt = ' %15s |   %3d   | %18s | %6.2f '

        np.random.seed(1234)
        for n in 2, 10, 30:

            P_sparse = _create_sparse_poisson2d(n)
            P_dense = P2d.A
            b = np.ones(n*n)

            # Use the generic dense solver.
            tm_start = time.clock()
            x_dense = scipy.linalg.solve(P_dense, b)
            tm_end = time.clock()
            tm_dense = tm_end - tm_start

            # Use the conjugate gradient solver with the sparse matrix.
            tm_start = time.clock()
            x_sparse = scipy.sparse.linalg.cg(P_sparse, b)
            tm_end = time.clock()
            tm_sparse = tm_end - tm_start

            # Check that the solutions are close to each other.
            assert_allclose(x_dense, x_sparse)

            # Write the rows.
            nrepeats = 1
            shape = (n*n, n*n)
            print(fmt % (shape, nrepeats, 'dense solve', tm_dense))
            print(fmt % (shape, nrepeats, 'sparse cg', tm_sparse))
        print()


if __name__ == '__main__':
    Tester().bench()
