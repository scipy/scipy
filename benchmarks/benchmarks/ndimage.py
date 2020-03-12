from __future__ import division, absolute_import, print_function

import numpy as np
import timeit

try:
    from scipy.ndimage import gabor_filter
except ImportError:
    pass

from .common import Benchmark


class GaborFilter(Benchmark):
    param_names = ['truncate']
    params = [1, 2, 3]

    def setup(self, truncate):
        np.random.seed(123456)
        self.image = np.random.randint(0, 256, size=[25, 25])
        self.volume = np.random.randint(0, 256, size=[25, 25, 25])

    def time_gabor_2d(self, truncate):
        gabor_filter(self.image, 2.5, 0, 0.1, 0, truncate=truncate)

    def time_gabor_3d(self, truncate):
        gabor_filter(self.volume, 5, 0, 0.1, 0, truncate=truncate)

