from __future__ import division, print_function, absolute_import

from math import sqrt

# Import testing parameters
try:
    from scipy.optimize._tstutils import (methods, mstrings, functions,
                                          fstrings, description)
except ImportError:
    pass


class Zeros(object):
    params = [
        fstrings,
        mstrings
    ]
    param_names = ['test function', 'solver']
    goal_time = 0.5

    def setup(self, func, meth):
        self.a = .5
        self.b = sqrt(3)
        repeat = 2000

        self.func = functions[fstrings.index(func)]
        self.meth = methods[mstrings.index(meth)]

    def time_zeros(self, func, meth):
        self.meth(self.func, self.a, self.b)
