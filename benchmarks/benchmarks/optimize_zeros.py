from __future__ import division, print_function, absolute_import

from math import sqrt, exp, cos, sin

# Import testing parameters
try:
    from scipy.optimize._tstutils import methods, mstrings, functions, fstrings
except ImportError:
    pass
from scipy.optimize import newton  # newton predates benchmarks

from .common import Benchmark


class Zeros(Benchmark):
    params = [
        fstrings,
        mstrings
    ]
    param_names = ['test function', 'solver']

    def setup(self, func, meth):
        self.a = .5
        self.b = sqrt(3)

        self.func = functions[fstrings.index(func)]
        self.meth = methods[mstrings.index(meth)]

    def time_zeros(self, func, meth):
        self.meth(self.func, self.a, self.b)


class Newton(Benchmark):
    params = [
        ['f1', 'f2'],
        ['newton', 'secant', 'halley']
    ]
    param_names = ['test function', 'solver']

    def setup(self, func, meth):
        self.x0 = 3
        self.f_1 = None
        self.f_2 = None
        if func == 'f1':
            self.f = lambda x: x ** 2 - 2 * x - 1
            if meth in ('newton', 'halley'):
                self.f_1 = lambda x: 2 * x - 2
            if meth == 'halley':
                self.f_2 = lambda x: 2.0 + 0 * x
        else:
            self.f = lambda x: exp(x) - cos(x)
            if meth in ('newton', 'halley'):
                self.f_1 = lambda x: exp(x) + sin(x)
            if meth == 'halley':
                self.f_2 = lambda x: exp(x) + cos(x)

    def time_newton(self, func, meth):
        newton(self.f, self.x0, args=(), fprime=self.f_1, fprime2=self.f_2)
