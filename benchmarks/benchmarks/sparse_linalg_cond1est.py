#!/usr/bin/env python3
# =============================================================================
#     File: sparse_linalg_cond1est.py
#  Created: 2025-09-02 12:41
#   Author: Bernie Roesler
#
"""Benchmark testing for scipy.sparse.linalg.cond1est."""
# =============================================================================

import numpy as np

from .common import Benchmark, safe_import

with safe_import():
    from scipy.sparse import random_array
    from scipy.sparse.linalg import splu


class Norm1EstBench(Benchmark):
    params = [
        [5000],
        [1e-5, 1e-3, 0.1],
    ]
    param_names = ["N", "density"]

    def setup(self, N, density):
        rng = np.random.default_rng(56)
        self.A = random_array(
            (N, N), density=density, format="csc", dtype=float, rng=rng
        )
        # Make the matrix non-singular
        self.A.setdiag(N * N)

    def time_splu(self, N, density):
        splu(self.A)

    def peakmem_splu(self, N, density):
        splu(self.A)


# =============================================================================
# =============================================================================
