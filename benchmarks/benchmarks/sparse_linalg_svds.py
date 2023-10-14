import os
import warnings
import numpy as np
from .common import Benchmark, safe_import

with safe_import():
    from scipy.linalg import svd
    from scipy.sparse.linalg import svds

# ensure that we are benchmarking a consistent outcome;
# (e.g. if the code wasn't able to find a solution accurately
# enough the timing of the benchmark would become useless).
msg = ("the benchmark code did not converge as expected, "
       "the timing is therefore useless")


class BenchSVDS(Benchmark):
    # Benchmark SVD using the MatrixMarket test matrices recommended by the
    # author of PROPACK at http://sun.stanford.edu/~rmunk/PROPACK/
    params = [
        [20],
        ["abb313", "illc1033", "illc1850", "qh1484", "rbs480a", "tols4000",
         "well1033", "well1850", "west0479", "west2021"],
        # TODO: re-include propack
        ['arpack', 'lobpcg', 'svd']  # 'propack' failing (Aug. 2023)
    ]
    param_names = ['k', 'problem', 'solver']

    def setup(self, k, problem, solver):
        if solver == 'svd' and problem == "tols4000":
            # skip: fails breaking the benchmark
            raise NotImplementedError()
        dir_path = os.path.dirname(os.path.realpath(__file__))
        datafile = os.path.join(dir_path, "svds_benchmark_files",
                                "svds_benchmark_files.npz")
        matrices = np.load(datafile, allow_pickle=True)
        self.A = matrices[problem][()]
        _, s, _ = svd(self.A.toarray(), full_matrices=False)
        self.top_singular_values = np.flip(s[:int(k/2)])
        self.tol = k * np.prod(self.A.shape) * np.finfo(float).eps
        self.rng = np.random.default_rng(0)

    def time_svds(self, k, problem, solver):
        if solver == 'svd':
            _ = svd(self.A.toarray(), full_matrices=False)
        else:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                _, s, _ = svds(self.A, k=k, solver=solver, random_state=self.rng,
                               maxiter = 50, tol=1e-6)
            accuracy = np.max(np.abs(1 - s[int(k/2):] / self.top_singular_values))
            assert accuracy < self.tol, msg
