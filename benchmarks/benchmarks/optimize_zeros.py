from math import sqrt, exp, cos, sin
import numpy as np

from .common import Benchmark, safe_import

# Import testing parameters
with safe_import():
    from scipy.optimize._tstutils import methods, mstrings, functions, fstrings
from scipy.optimize import newton  # newton predates benchmarks


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


class NewtonArray(Benchmark):
    params = [['loop', 'array'], ['newton', 'secant', 'halley']]
    param_names = ['vectorization', 'solver']

    def setup(self, vec, meth):
        if vec == 'loop':
            if meth == 'newton':
                self.fvec = lambda f, x0, args, fprime, fprime2: [
                    newton(f, x, args=(a0, a1) + args[2:], fprime=fprime)
                    for (x, a0, a1) in zip(x0, args[0], args[1])
                ]
            elif meth == 'halley':
                self.fvec = lambda f, x0, args, fprime, fprime2: [
                    newton(
                        f, x, args=(a0, a1) + args[2:], fprime=fprime,
                        fprime2=fprime2
                    ) for (x, a0, a1) in zip(x0, args[0], args[1])
                ]
            else:
                self.fvec = lambda f, x0, args, fprime, fprime2: [
                    newton(f, x, args=(a0, a1) + args[2:]) for (x, a0, a1)
                    in zip(x0, args[0], args[1])
                ]
        else:
            if meth == 'newton':
                self.fvec = lambda f, x0, args, fprime, fprime2: newton(
                    f, x0, args=args, fprime=fprime
                )
            elif meth == 'halley':
                self.fvec = newton
            else:
                self.fvec = lambda f, x0, args, fprime, fprime2: newton(
                    f, x0, args=args
                )

    def time_array_newton(self, vec, meth):

        def f(x, *a):
            b = a[0] + x * a[3]
            return a[1] - a[2] * (np.exp(b / a[5]) - 1.0) - b / a[4] - x

        def f_1(x, *a):
            b = a[3] / a[5]
            return -a[2] * np.exp(a[0] / a[5] + x * b) * b - a[3] / a[4] - 1

        def f_2(x, *a):
            b = a[3] / a[5]
            return -a[2] * np.exp(a[0] / a[5] + x * b) * b ** 2

        a0 = np.array([
            5.32725221, 5.48673747, 5.49539973,
            5.36387202, 4.80237316, 1.43764452,
            5.23063958, 5.46094772, 5.50512718,
            5.42046290
        ])
        a1 = (np.sin(range(10)) + 1.0) * 7.0
        args = (a0, a1, 1e-09, 0.004, 10, 0.27456)
        x0 = [7.0] * 10
        self.fvec(f, x0, args=args, fprime=f_1, fprime2=f_2)
