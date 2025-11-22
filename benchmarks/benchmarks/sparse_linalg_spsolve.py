"""Check the speed of the sparse solve function with many RHS."""

import numpy as np

from .common import Benchmark, safe_import

with safe_import():
    from scipy import sparse
    from scipy.sparse.linalg import LaplacianNd, spsolve


class SolveSparseRHS(Benchmark):
    param_names = ["Nsq", "K", "density", "use_umfpack"]
    params = [
        [5000],  # (Nsq, Nsq) matrix for (N, N) grid
        [2000],  # (Nsq, K) RHS
        [0.01, 0.1, 0.5],  # density of the sparse RHS
        # [1, 100, 1000],  # batch size
        # [True, False],  # use umfpack or not
        [False],  # use umfpack or not
    ]

    def setup(self, Nsq, K, density, use_umfpack):
        Ng = np.sqrt(Nsq).astype(int)
        self.A = -LaplacianNd((Ng, Ng), dtype=float).tosparse().tocsc()
        N = self.A.shape[0]
        self.A[-1, -1] += 1  # make A non-singular
        self.b = sparse.random_array((N, K), density=density, format="csc", dtype=float)
        # self.batch_size = batch_size
        self.use_umfpack = use_umfpack

    def time_solve(self, Nsq, K, density, use_umfpack):
        spsolve(
            self.A, self.b, use_umfpack=self.use_umfpack
        )

    def peakmem_solve(self, Nsq, K, density, use_umfpack):
        spsolve(
            self.A, self.b, use_umfpack=self.use_umfpack
        )
