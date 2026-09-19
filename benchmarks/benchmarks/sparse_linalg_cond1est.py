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
    from scipy.sparse.linalg import LaplacianNd, cond1est


class Cond1EstBench(Benchmark):
    params = [
        list(np.sqrt(np.r_[500, 1000, 2000, 5000, 10_000]).astype(int)),
    ]
    param_names = ["sqrtN"]

    def setup(self, sqrtN):
        self.A = LaplacianNd((sqrtN, sqrtN), dtype=float).tosparse().tocsc()

    def time_cond1est(self, sqrtN):
        cond1est(self.A)

    def peakmem_cond1est(self, sqrtN):
        cond1est(self.A)


# =============================================================================
# =============================================================================
