from __future__ import division, absolute_import, print_function

import numpy as np

try:
    from scipy.special import ai_zeros, bi_zeros, erf
except ImportError:
    pass

from .common import Benchmark


class Airy(Benchmark):
    def time_ai_zeros(self):
        ai_zeros(100000)

    def time_bi_zeros(self):
        bi_zeros(100000)


class Erf(Benchmark):
    def setup(self, *args):
        self.rand = np.random.rand(1e5)

    def time_real(self, offset):
        erf(self.rand + offset)

    time_real.params = [0.0, 2.0]
    time_real.param_names = ['offset']
