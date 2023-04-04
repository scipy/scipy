import numpy as np

from .common import Benchmark, with_attributes, safe_import

with safe_import():
    from scipy.special import (
        ai_zeros,
        bi_zeros,
        erf,
        expn,
        factorial,
        factorial2,
        factorialk,
    )
with safe_import():
    # wasn't always in scipy.special, so import separately
    from scipy.special import comb
with safe_import():
    from scipy.special import loggamma


class Airy(Benchmark):
    def time_ai_zeros(self):
        ai_zeros(100000)

    def time_bi_zeros(self):
        bi_zeros(100000)


class Erf(Benchmark):
    def setup(self, *args):
        self.rand = np.random.rand(100000)

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


class Loggamma(Benchmark):

    def setup(self):
        x, y = np.logspace(3, 5, 10), np.logspace(3, 5, 10)
        x, y = np.meshgrid(x, y)
        self.large_z = x + 1j*y

    def time_loggamma_asymptotic(self):
        loggamma(self.large_z)


class Expn(Benchmark):

    def setup(self):
        n, x = np.arange(50, 500), np.logspace(0, 20, 100)
        n, x = np.meshgrid(n, x)
        self.n, self.x = n, x

    def time_expn_large_n(self):
        expn(self.n, self.x)


class Factorial(Benchmark):
    def setup(self, *args):
        self.positive_ints = np.arange(10, 111, step=20, dtype=int)
        self.negative_ints = -1 * self.positive_ints
        self.positive_floats = np.linspace(100.2, 1000.8, num=10)
        self.negative_floats = -1 * self.positive_floats

    @with_attributes(params=[(100, 1000, 10000)],
                     param_names=['n'])
    def time_factorial_exact_false_scalar_positive_int(self, n):
        factorial(n, exact=False)

    def time_factorial_exact_false_scalar_negative_int(self):
        factorial(-10000, exact=False)

    @with_attributes(params=[(100.8, 1000.3, 10000.5)],
                     param_names=['n'])
    def time_factorial_exact_false_scalar_positive_float(self, n):
        factorial(n, exact=False)

    def time_factorial_exact_false_scalar_negative_float(self):
        factorial(-10000.8, exact=False)

    def time_factorial_exact_false_array_positive_int(self):
        factorial(self.positive_ints, exact=False)

    def time_factorial_exact_false_array_negative_int(self):
        factorial(self.negative_ints, exact=False)

    def time_factorial_exact_false_array_positive_float(self):
        factorial(self.positive_floats, exact=False)

    def time_factorial_exact_false_array_negative_float(self):
        factorial(self.negative_floats, exact=False)

    @with_attributes(params=[(100, 200, 400)],
                     param_names=['n'])
    def time_factorial_exact_true_scalar_positive_int(self, n):
        factorial(n, exact=True)

    def time_factorial_exact_true_scalar_negative_int(self):
        factorial(-10000, exact=True)

    def time_factorial_exact_true_scalar_negative_float(self):
        factorial(-10000.8, exact=True)

    def time_factorial_exact_true_array_positive_int(self):
        factorial(self.positive_ints, exact=True)

    def time_factorial_exact_true_array_negative_int(self):
        factorial(self.negative_ints, exact=True)

    def time_factorial_exact_true_array_negative_float(self):
        factorial(self.negative_floats, exact=True)


class Factorial2(Benchmark):
    def setup(self, *args):
        self.positive_ints = np.arange(100, 201, step=20, dtype=int)
        self.negative_ints = -1 * self.positive_ints

    @with_attributes(params=[(100, 200, 400)],
                     param_names=['n'])
    def time_factorial2_exact_false_scalar_positive_int(self, n):
        factorial2(n, exact=False)

    def time_factorial2_exact_false_scalar_negative_int(self):
        factorial2(-10000, exact=False)

    def time_factorial2_exact_false_array_positive_int(self):
        factorial2(self.positive_ints, exact=False)

    def time_factorial2_exact_false_array_negative_int(self):
        factorial2(self.negative_ints, exact=False)

    @with_attributes(params=[(100, 200, 400)],
                     param_names=['n'])
    def time_factorial2_exact_true_scalar_positive_int(self, n):
        factorial2(n, exact=True)

    def time_factorial2_exact_true_scalar_negative_int(self):
        factorial2(-10000, exact=True)


class FactorialK(Benchmark):
    @with_attributes(params=[(100, 500), range(1, 10)],
                     param_names=['n', 'k'])
    def time_factorialk_exact_true_scalar_positive_int(self, n, k):
        factorialk(n, k, exact=True)

    def time_factorialk_exact_false_scalar_negative_int(self):
        factorialk(-10000, 3, exact=True)
