from __future__ import division, absolute_import, print_function

from itertools import product

import numpy as np

try:
    from scipy.signal import convolve2d, correlate2d, lti, lsim, lsim2
except ImportError:
    pass

from .common import Benchmark


class Convolve2D(Benchmark):
    def setup(self):
        np.random.seed(1234)

        # sample a bunch of pairs of 2d arrays
        pairs = []
        for ma, na, mb, nb in product((1, 2, 8, 13, 30), repeat=4):
            a = np.random.randn(ma, na)
            b = np.random.randn(mb, nb)
            pairs.append((a, b))
        self.pairs = pairs

    def time_convolutions(self):
        fns = (convolve2d, correlate2d)
        modes = ('full', 'valid', 'same')
        boundaries = ('fill', 'wrap', 'symm')

        # compute 2d convolutions and correlations for each 2d array pair
        for a, b in self.pairs:
            for fn, mode, boundary in product(fns, modes, boundaries):
                if mode == 'valid':
                    if b.shape[0] > a.shape[0] or b.shape[1] > a.shape[1]:
                        continue
                fn(a, b, mode=mode, boundary=boundary)


class LTI(Benchmark):
    def setup(self):
        self.system = lti(1.0, [1, 0, 1])
        self.t = np.arange(0, 1000, 0.5)
        self.u = np.sin(2 * self.t)

    def time_lsim(self):
        lsim(self.system, self.u, self.t)

    def time_lsim2(self):
        lsim2(self.system, self.u, self.t)
