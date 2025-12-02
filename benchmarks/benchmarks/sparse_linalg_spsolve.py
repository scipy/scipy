"""Check the speed of the sparse solve function with many RHS."""

import numpy as np

from .common import Benchmark, safe_import

with safe_import():
    from scipy import sparse
    from scipy.sparse.linalg import LaplacianNd, spsolve


class SolveSparseRHS(Benchmark):
    param_names = ["Nsq", "K", "density", "batch_size", "use_umfpack"]
    params = [
        [500],  # (Nsq, Nsq) matrix for (N, N) grid
        [5_003],  # (Nsq, K) RHS matrix, arbitrary K
        [0.01, 0.1, 0.5],  # density of the sparse RHS
        [1, 10, 100, 6000],  # batch size
        [True, False],  # use umfpack or not
    ]

    def setup(self, Nsq, K, density, batch_size, use_umfpack):
        Ng = np.sqrt(Nsq).astype(int)
        self.A = -LaplacianNd((Ng, Ng), dtype=float).tosparse().tocsc()
        N = self.A.shape[0]
        self.A[-1, -1] += 1  # make A non-singular
        self.b = sparse.random_array((N, K), density=density, format="csc", dtype=float)
        self.batch_size = batch_size
        self.use_umfpack = use_umfpack

    def _run_solver(self, Nsq, K, density, batch_size, use_umfpack):
        try:
            spsolve(
                self.A, self.b, batch_size=self.batch_size, use_umfpack=self.use_umfpack
            )
        except TypeError:
            # older versions do not have batch_size
            if self.batch_size == 1:
                spsolve(self.A, self.b, use_umfpack=self.use_umfpack)
            else:
                pass  # skip benchmark

    def time_solve(self, Nsq, K, density, batch_size, use_umfpack):
        self._run_solver(Nsq, K, density, batch_size, use_umfpack)

    def peakmem_solve(self, Nsq, K, density, batch_size, use_umfpack):
        self._run_solver(Nsq, K, density, batch_size, use_umfpack)
