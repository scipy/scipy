import os
import warnings
import numpy as np
from .common import Benchmark, safe_import
from asv_runner.benchmarks.mark import SkipNotImplemented

with safe_import():
    from scipy.linalg import svd
    from scipy.sparse.linalg import svds


class BenchSVDS(Benchmark):
    # Benchmark SVDS using the MatrixMarket test matrices recommended by the
    # author of PROPACK at http://sun.stanford.edu/~rmunk/PROPACK/
    # The dense SVD solve is benchmarked as the baseline.
    # `arpack` convergence uniformly for all singular values, while
    # `lobpcg` extreme singular values may converge much faster, so
    # the accuracy is only checked on the top k/2 singular values.
    # The check uses the relative error since for some tested matrices
    # the maximal singular values are very large due to poor matrix scaling.
    # On 7/19/24, all accuracy checks pass with ``maxiter = 200, tol=1e-6``.
    # If changes are made to a solver, the accuracy check may fail
    # resulting in the failure of the benchmark for the solver.

    # problems ``tols4000`` and ``west2021`` take too long and may fail
    # with some solvers so excluded from the benchmark.
    params = [
        [20],
        [
            "abb313",
            "illc1033",
            "illc1850",
            "qh1484",
            "rbs480a",
            "well1033",
            "well1850",
            "west0479",
        ],
        ["propack", "arpack", "lobpcg", "svd"],
    ]
    param_names = ['k', 'problem', 'solver']

    def __init__(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        datafile = os.path.join(dir_path, "svds_benchmark_files",
                                "svds_benchmark_files.npz")
        self.matrices = np.load(datafile, allow_pickle=True)

    def setup(self, k, problem, solver):
        self.A = self.matrices[problem][()]
        _, s, _ = svd(self.A.toarray(), full_matrices=False)
        self.top_singular_values = np.flip(s[:int(k/2)])
        self.tol = k * np.max(self.A.shape) * np.finfo(float).eps
        self.rng = np.random.default_rng(98360967947894649386)

    def time_svds(self, k, problem, solver):
        # The 'svd' solver find all ``m = np.min(self.A.shape) >> k``
        # singular pairs but may still be expected to outperform
        # the sparse solvers benchmarked here if m is small enough.
        # It is commonly thus used as a baseline for comparisons.
        if solver == 'svd':
            svd(self.A.toarray(), full_matrices=False)
        else:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # parameters `maxiter` and `tol` are tuned for fair comparison
                _, s, _ = svds(self.A, k=k, solver=solver, random_state=self.rng,
                               maxiter = 200, tol=1e-6)
            accuracy = np.sum(np.abs(1 - s[int(k/2):] / self.top_singular_values))
            # ensure that we are benchmarking a consistent outcome;
            # (e.g. if the code wasn't able to find a solution accurately
            # enough the timing of the benchmark would become useless).
            if accuracy > self.tol:
                # not convergence sufficiently for timing to be relevant
                raise SkipNotImplemented("Insufficient accuracy achieved")
