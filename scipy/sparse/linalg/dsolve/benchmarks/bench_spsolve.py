"""
Check the speed of the sparse solver.

"""
from __future__ import division, print_function, absolute_import

import time

import numpy as np
from numpy.testing import Tester, TestCase, assert_allclose, assert_equal

from scipy import linalg, sparse
from scipy.sparse.linalg import spsolve


def _create_sparse_poisson1d(n):
    # Make Gilbert Strang's favorite matrix
    # http://www-math.mit.edu/~gs/PIX/cupcakematrix.jpg
    P1d = sparse.diags([[-1]*(n-1), [2]*n, [-1]*(n-1)], [-1, 0, 1])
    assert_equal(P1d.shape, (n, n))
    return P1d


def _create_sparse_poisson2d(n):
    P1d = _create_sparse_poisson1d(n)
    P2d = sparse.kronsum(P1d, P1d)
    assert_equal(P2d.shape, (n*n, n*n))
    return P2d


class BenchmarkConjuateGradientSolver(TestCase):

    def bench_spsolve(self):
        np.random.seed(1234)
        n = 40
        P_sparse = _create_sparse_poisson2d(n).tocsc()
        P_dense = P_sparse.A
        rows = P_dense.shape[0]

        # print headers and define the column formats
        print()
        print('    dense vs. sparse solve of Ax=b with %d rows' % rows)
        print('==============================================================')
        print('   columns | repeats |      operation     |   time   ')
        print('           |         |                    | (seconds)')
        print('--------------------------------------------------------------')
        fmt = '   %5s   |   %3d   | %18s | %6.2f '

        # Define the sizes of A and b in the equation Ax=b.
        # The shape of A will be (n*n, n*n), and the shape of b will
        # vary between (n*n, 1) and (n*n, 1000).
        # The time to construct these matrices will not count
        # towards the time reported in the benchmarks.
        b_column_counts = (10, 30, 100, 300, 1000)
        b_full = np.random.randn(rows, max(b_column_counts))
        b_dense_list = [b_full[:, :c] for c in b_column_counts]
        b_sparse_list = [sparse.csr_matrix(b) for b in b_dense_list]
        repeats = 10
        for i, b_column_count in enumerate(b_column_counts):
            b_dense = b_dense_list[i]
            b_sparse = b_sparse_list[i]

            # Use the generic dense solver.
            tm_start = time.clock()
            for i in range(repeats):
                x_dense = linalg.solve(P_dense, b_dense)
            tm_end = time.clock()
            tm_dense = tm_end - tm_start
            print(fmt % (b_column_count, repeats, 'dense solve', tm_dense))

            # Use the sparse solver.
            tm_start = time.clock()
            for i in range(repeats):
                x_sparse = spsolve(P_sparse, b_sparse)
            tm_end = time.clock()
            tm_sparse = tm_end - tm_start
            print(fmt % (b_column_count, repeats, 'sparse solve', tm_sparse))

            # Check that the solutions are close to each other.
            assert_allclose(x_dense, x_sparse.A, rtol=1e-4)

        print()


if __name__ == '__main__':
    Tester().bench()

