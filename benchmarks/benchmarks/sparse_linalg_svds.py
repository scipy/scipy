import os
import numpy as np
from .common import Benchmark, safe_import

with safe_import():
    from scipy.linalg import svd
    from scipy.sparse.linalg import svds


class BenchSVDS(Benchmark):
    # Benchmark SVD using the MatrixMarket test matrices recommended by the
    # author of PROPACK at http://sun.stanford.edu/~rmunk/PROPACK/
    # tests silently faling by `lobpcg` removed: west0479
    params = [
        [25],
        ["abb313", "illc1033", "illc1850", "qh1484", "rbs480a", "tols4000",
         "well1033", "well1850", "west0479", "west2021"],
        # TODO: re-include propack
        ['arpack', 'lobpcg', 'svd']  # 'propack' failing (Aug. 2023)
    ]
    param_names = ['k', 'problem', 'solver']

    def setup(self, k, problem, solver):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        datafile = os.path.join(dir_path, "svds_benchmark_files",
                                "svds_benchmark_files.npz")
        matrices = np.load(datafile, allow_pickle=True)
        self.A = matrices[problem][()]

    def time_svds(self, k, problem, solver):
        rng = np.random.default_rng(0)
        if solver = 'svd':
            _, s, _ = svd(self.A.toarray(), full_matrices=False)
            else
            _, s, _ = svds(self.A, k=k, solver=solver, random_state=rng)
