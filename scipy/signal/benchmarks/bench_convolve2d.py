"""
Check the speed of 2d convolution.

"""
from __future__ import division, print_function, absolute_import

from itertools import product
import time

import numpy as np
from numpy.testing import Tester

from scipy.signal import convolve2d, correlate2d


def bench_convolve2d():
    np.random.seed(1234)

    # sample a bunch of pairs of 2d arrays
    tm_start = time.clock()
    pairs = []
    for ma, na, mb, nb in product((1, 2, 8, 13, 30), repeat=4):
        a = np.random.randn(ma, na)
        b = np.random.randn(mb, nb)
        pairs.append((a, b))
    tm_end = time.clock()
    tm_total = tm_end - tm_start
    print('time to sample random pairs of 2d arrays:')
    print(tm_total)
    print()

    fns = (convolve2d, correlate2d)
    modes = ('full', 'valid', 'same')
    boundaries = ('fill', 'wrap', 'symm')

    # compute 2d convolutions and correlations for each 2d array pair
    tm_start = time.clock()
    for a, b in pairs:
        for fn, mode, boundary in product(fns, modes, boundaries):
            if mode == 'valid':
                if b.shape[0] > a.shape[0] or b.shape[1] > a.shape[1]:
                    continue
            fn(a, b, mode=mode, boundary=boundary)
    tm_end = time.clock()
    tm_total = tm_end - tm_start
    print('time to compute 2d convolutions:')
    print(tm_total)
    print()


if __name__ == '__main__':
    Tester().bench()
