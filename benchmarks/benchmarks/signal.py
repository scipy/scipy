from __future__ import division, absolute_import, print_function

from itertools import product

import numpy as np

try:
    import scipy.signal as signal
except ImportError:
    pass

from .common import Benchmark


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
        for ma, na, mb, nb in product((1, 2, 8, 13, 30), repeat=4):
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
    param_names = ['mode']    
    params = [
        ['full', 'valid', 'same']
    ]

    def setup(self, mode):
        np.random.seed(1234)
        # sample a bunch of pairs of 2d arrays
        pairs = []
        for ma, nb in product((1, 2, 8, 13, 30, 36, 50, 75), repeat=2):
            a = np.random.randn(ma)
            b = np.random.randn(nb)
            pairs.append((a, b))
        self.pairs = pairs

    def time_convolve2d(self, mode):
        for a, b in self.pairs:
            if b.shape[0] > a.shape[0]:
                continue
            signal.fftconvolve(a, b, mode=mode)


class Convolve(Benchmark):
    param_names = ['mode']    
    params = [
        ['full', 'valid', 'same']
    ]

    def setup(self, mode):
        np.random.seed(1234)
        # sample a bunch of pairs of 2d arrays
        pairs = []
        for ma, nb in product((1, 2, 8, 13, 30, 36, 50, 75), repeat=2):
            a = np.random.randn(ma)
            b = np.random.randn(nb)
            pairs.append((a, b))
        self.pairs = pairs

    def time_convolve(self, mode):
        for a, b in self.pairs:
            if b.shape[0] > a.shape[0]:
                continue
            signal.convolve(a, b, mode=mode)

    def time_correlate(self, mode):
        for a, b in self.pairs:
            if b.shape[0] > a.shape[0]:
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
