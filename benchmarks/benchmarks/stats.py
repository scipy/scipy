from __future__ import division, absolute_import, print_function

import warnings

import numpy as np

try:
    from scipy.stats import anderson_ksamp
except ImportError:
    pass

from .common import Benchmark


class Anderson_KSamp(Benchmark):
    def setup(self, *args):
        self.rand = [np.random.normal(loc=i, size=1000) for i in range(3)]

    def time_anderson_ksamp(self):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            anderson_ksamp(self.rand)
