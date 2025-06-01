import numpy as np

from .common import Benchmark, with_attributes, safe_import

with safe_import():
    from scipy.special import ai_zeros, bi_zeros, erf, expn
with safe_import():
    # wasn't always in scipy.special, so import separately
    from scipy.special import comb
with safe_import():
    from scipy.special import loggamma
with safe_import():
    from scipy.special import abs_sq


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


class AbsSq(Benchmark):

    def setup(self):
        rng = np.random.default_rng()
        self.d64 = rng.random((500, 500))
        self.d32 = self.d64.astype(np.float32)
        self.d128 = self.d64.astype(np.float128)
        self.z128 = self.d64 + 1j * rng.random((500, 500))
        self.z64 = self.z128.astype(np.complex64)
        self.z256 = self.z128.astype(np.complex256)

    def time_long_double_complex(self):
        abs_sq(self.z256)

    def time_double_complex(self):
        abs_sq(self.z128)

    def time_float_complex(self):
        abs_sq(self.z64)

    def time_long_double(self):
        abs_sq(self.d128)

    def time_double(self):
        abs_sq(self.d64)

    def time_float(self):
        abs_sq(self.d32)
