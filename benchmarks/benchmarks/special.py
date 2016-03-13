from __future__ import division, absolute_import, print_function

import numpy as np

try:
    from scipy.special import ai_zeros, bi_zeros, erf
except ImportError:
    pass

try:
    # wasn't always in scipy.special, so import separately
    from scipy.special import comb
except ImportError:
    pass

from .common import Benchmark, with_attributes


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


class Comb(Benchmark):

    def setup(self, *args):
        self.N = np.arange(1, 1000, 50)
        self.k = np.arange(1, 1000, 50)

    @with_attributes(params=[(10, 100, 1000, 10000), (1, 10, 100)],
                     param_names=['N', 'k'])
    def time_comb_exact(self, N, k):
        comb(N, k, exact=True)

    def time_comb_float(self):
        comb(self.N[:,None], self.k[None,:])
