from __future__ import division, absolute_import, print_function

import sys
import numpy.linalg as nl

from numpy.testing import assert_
from numpy.random import rand

try:
    import scipy.linalg as sl
except ImportError:
    pass


def random(size):
    return rand(*size)


class Bench(object):
    params = [
        [20, 100, 500, 1000],
        ['contig', 'nocont'],
        ['numpy', 'scipy']
    ]
    param_names = ['size', 'contiguous', 'module']
    goal_time = 0.5

    def setup(self, size, contig, module):
        a = random([size,size])
        # larger diagonal ensures non-singularity:
        for i in range(size):
            a[i,i] = 10*(.1+a[i,i])
        b = random([size])

        if contig != 'contig':
            a = a[-1::-1,-1::-1]  # turn into a non-contiguous array
            assert_(not a.flags['CONTIGUOUS'])

        self.a = a
        self.b = b

    def time_solve(self, size, contig, module):
        if module == 'numpy':
            nl.solve(self.a, self.b)
        else:
            sl.solve(self.a, self.b)

    def time_inv(self, size, contig, module):
        if module == 'numpy':
            nl.inv(self.a)
        else:
            sl.inv(self.a)

    def time_det(self, size, contig, module):
        if module == 'numpy':
            nl.det(self.a)
        else:
            sl.det(self.a)

    def time_eigvals(self, size, contig, module):
        if module == 'numpy':
            nl.eigvals(self.a)
        else:
            sl.eigvals(self.a)

    def time_svd(self, size, contig, module):
        if module == 'numpy':
            nl.svd(self.a)
        else:
            sl.svd(self.a)
