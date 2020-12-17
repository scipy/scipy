from itertools import product

import numpy as np
from .common import Benchmark, safe_import

with safe_import():
    import scipy.signal as signal


class Resample(Benchmark):

    # Some slow (prime), some fast (in radix)
    param_names = ['N', 'num']
    params = [[977, 9973, 2 ** 14, 2 ** 16]] * 2

    def setup(self, N, num):
        x = np.linspace(0, 10, N, endpoint=False)
        self.y = np.cos(-x**2/6.0)

    def time_complex(self, N, num):
        signal.resample(self.y + 0j, num)

    def time_real(self, N, num):
        signal.resample(self.y, num)


class CalculateWindowedFFT(Benchmark):

    def setup(self):
        np.random.seed(5678)
        # Create some long arrays for computation
        x = np.random.randn(2**20)
        y = np.random.randn(2**20)
        self.x = x
        self.y = y

    def time_welch(self):
        signal.welch(self.x)

    def time_csd(self):
        signal.csd(self.x, self.y)

    def time_periodogram(self):
        signal.periodogram(self.x)

    def time_spectrogram(self):
        signal.spectrogram(self.x)

    def time_coherence(self):
        signal.coherence(self.x, self.y)


class Convolve2D(Benchmark):
    param_names = ['mode', 'boundary']
    params = [
        ['full', 'valid', 'same'],
        ['fill', 'wrap', 'symm']
    ]

    def setup(self, mode, boundary):
        np.random.seed(1234)
        # sample a bunch of pairs of 2d arrays
        pairs = []
        for ma, na, mb, nb in product((8, 13, 30, 36), repeat=4):
            a = np.random.randn(ma, na)
            b = np.random.randn(mb, nb)
            pairs.append((a, b))
        self.pairs = pairs

    def time_convolve2d(self, mode, boundary):
        for a, b in self.pairs:
            if mode == 'valid':
                if b.shape[0] > a.shape[0] or b.shape[1] > a.shape[1]:
                    continue
            signal.convolve2d(a, b, mode=mode, boundary=boundary)

    def time_correlate2d(self, mode, boundary):
        for a, b in self.pairs:
            if mode == 'valid':
                if b.shape[0] > a.shape[0] or b.shape[1] > a.shape[1]:
                    continue
            signal.correlate2d(a, b, mode=mode, boundary=boundary)


class FFTConvolve(Benchmark):
    param_names = ['mode', 'size']
    params = [
        ['full', 'valid', 'same'],
        [(a,b) for a,b in product((1, 2, 8, 36, 60, 150, 200, 500), repeat=2)
         if b <= a]
    ]

    def setup(self, mode, size):
        np.random.seed(1234)
        # sample a bunch of pairs of 2d arrays
        pairs = []
        self.a = np.random.randn(size[0])
        self.b = np.random.randn(size[1])

    def time_convolve2d(self, mode, size):
        signal.fftconvolve(self.a, self.b, mode=mode)


class OAConvolve(Benchmark):
    param_names = ['mode', 'size']
    params = [
        ['full', 'valid', 'same'],
        [(a, b) for a, b in product((40, 200, 3000), repeat=2)
         if b < a]
    ]

    def setup(self, mode, size):
        np.random.seed(1234)
        # sample a bunch of pairs of 2d arrays
        pairs = []
        self.a = np.random.randn(size[0])
        self.b = np.random.randn(size[1])

    def time_convolve2d(self, mode, size):
        signal.oaconvolve(self.a, self.b, mode=mode)


class Convolve(Benchmark):
    param_names = ['mode']
    params = [
        ['full', 'valid', 'same']
    ]

    def setup(self, mode):
        np.random.seed(1234)
        # sample a bunch of pairs of 2d arrays
        pairs = {'1d': [], '2d': []}
        for ma, nb in product((1, 2, 8, 13, 30, 36, 50, 75), repeat=2):
            a = np.random.randn(ma)
            b = np.random.randn(nb)
            pairs['1d'].append((a, b))

        for n_image in [256, 512, 1024]:
            for n_kernel in [3, 5, 7]:
                x = np.random.randn(n_image, n_image)
                h = np.random.randn(n_kernel, n_kernel)
                pairs['2d'].append((x, h))
        self.pairs = pairs

    def time_convolve(self, mode):
        for a, b in self.pairs['1d']:
            if b.shape[0] > a.shape[0]:
                continue
            signal.convolve(a, b, mode=mode)

    def time_convolve2d(self, mode):
        for a, b in self.pairs['2d']:
            if mode == 'valid':
                if b.shape[0] > a.shape[0] or b.shape[1] > a.shape[1]:
                    continue
            signal.convolve(a, b, mode=mode)

    def time_correlate(self, mode):
        for a, b in self.pairs['1d']:
            if b.shape[0] > a.shape[0]:
                continue
            signal.correlate(a, b, mode=mode)

    def time_correlate2d(self, mode):
        for a, b in self.pairs['2d']:
            if mode == 'valid':
                if b.shape[0] > a.shape[0] or b.shape[1] > a.shape[1]:
                    continue
            signal.correlate(a, b, mode=mode)


class LTI(Benchmark):

    def setup(self):
        self.system = signal.lti(1.0, [1, 0, 1])
        self.t = np.arange(0, 100, 0.5)
        self.u = np.sin(2 * self.t)

    def time_lsim(self):
        signal.lsim(self.system, self.u, self.t)

    def time_lsim2(self):
        signal.lsim2(self.system, self.u, self.t)

    def time_step(self):
        signal.step(self.system, T=self.t)

    def time_impulse(self):
        signal.impulse(self.system, T=self.t)

    def time_bode(self):
        signal.bode(self.system)


class Upfirdn1D(Benchmark):
    param_names = ['up', 'down']
    params = [
        [1, 4],
        [1, 4]
    ]

    def setup(self, up, down):
        np.random.seed(1234)
        # sample a bunch of pairs of 2d arrays
        pairs = []
        for nfilt in [8, ]:
            for n in [32, 128, 512, 2048]:
                h = np.random.randn(nfilt)
                x = np.random.randn(n)
                pairs.append((h, x))

        self.pairs = pairs

    def time_upfirdn1d(self, up, down):
        for h, x in self.pairs:
            signal.upfirdn(h, x, up=up, down=down)


class Upfirdn2D(Benchmark):
    param_names = ['up', 'down', 'axis']
    params = [
        [1, 4],
        [1, 4],
        [0, -1],
    ]

    def setup(self, up, down, axis):
        np.random.seed(1234)
        # sample a bunch of pairs of 2d arrays
        pairs = []
        for nfilt in [8, ]:
            for n in [32, 128, 512]:
                h = np.random.randn(nfilt)
                x = np.random.randn(n, n)
                pairs.append((h, x))

        self.pairs = pairs

    def time_upfirdn2d(self, up, down, axis):
        for h, x in self.pairs:
            signal.upfirdn(h, x, up=up, down=down, axis=axis)


class FIRLS(Benchmark):
    param_names = ['n', 'edges']
    params = [
        [21, 101, 1001, 2001],
        [(0.1, 0.9), (0.01, 0.99)],
        ]

    def time_firls(self, n, edges):
        signal.firls(n, (0,) + edges + (1,), [1, 1, 0, 0])
