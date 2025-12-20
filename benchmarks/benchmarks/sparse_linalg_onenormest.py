"""Compare the speed of exact one-norm calculation vs. its estimation.
"""
import numpy as np

from .common import Benchmark, safe_import

with safe_import():
    import scipy.sparse
    import scipy.sparse.linalg
    from scipy.sparse.linalg import aslinearoperator, onenormest


class BenchmarkOneNormEst(Benchmark):
    params = [
        [2, 3, 5, 10, 30, 100, 300, 500, 1000, 1e4, 1e5, 1e6],
        ['exact', 'onenormest']
    ]
    param_names = ['n', 'solver']

    def setup(self, n, solver):
        self.rng = np.random.default_rng(1234)
        nrepeats = 100
        shape = (int(n), int(n))

        if solver == 'exact' and n >= 300:
            # skip: slow, and not useful to benchmark
            raise NotImplementedError()

        if n <= 1000:
            # Sample the matrices.
            self.matrices = []
            for i in range(nrepeats):
                M = self.rng.standard_normal(shape)
                self.matrices.append(M)
        else:
            max_nnz = 100000
            nrepeats = 1

            self.matrices = []
            for i in range(nrepeats):
                M = scipy.sparse.rand(
                    shape[0],
                    shape[1],
                    min(max_nnz/(shape[0]*shape[1]), 1e-5),
                    random_state=self.rng,
                )
                self.matrices.append(M)

    def time_onenormest(self, n, solver):
        if solver == 'exact':
            # Get the exact values of one-norms of squares.
            for M in self.matrices:
                M.dot(M)
                scipy.sparse.linalg._matfuncs._onenorm(M)
        elif solver == 'onenormest':
            # Get the estimates of one-norms of squares.
            for M in self.matrices:
                onenormest(aslinearoperator(M)**2, rng=self.rng)
