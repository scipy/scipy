from __future__ import division, absolute_import, print_function

import numpy as np

try:
    from scipy.special import ai_zeros, bi_zeros
except ImportError:
    pass

from .common import Benchmark


class Airy(Benchmark):
    def time_ai_zeros(self):
        ai_zeros(100000)

    def time_bi_zeros(self):
        bi_zeros(100000)

