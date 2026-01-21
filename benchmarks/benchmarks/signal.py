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
        rng = np.random.default_rng(5678)
        # Create some long arrays for computation
        x = rng.standard_normal(2**20)
        y = rng.standard_normal(2**20)
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
        rng = np.random.default_rng(1234)
        # sample a bunch of pairs of 2d arrays
        pairs = []
        for ma, na, mb, nb in product((8, 13, 30, 36), repeat=4):
            a = rng.standard_normal((ma, na))
            b = rng.standard_normal((mb, nb))
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
        rng = np.random.default_rng(1234)
        self.a = rng.standard_normal(size[0])
        self.b = rng.standard_normal(size[1])

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
        rng = np.random.default_rng(1234)
        self.a = rng.standard_normal(size[0])
        self.b = rng.standard_normal(size[1])

    def time_convolve2d(self, mode, size):
        signal.oaconvolve(self.a, self.b, mode=mode)


class Convolve(Benchmark):
    param_names = ['mode']
    params = [
        ['full', 'valid', 'same']
    ]

    def setup(self, mode):
        rng = np.random.default_rng(1234)
        # sample a bunch of pairs of 2d arrays
        pairs = {'1d': [], '2d': []}
        for ma, nb in product((1, 2, 8, 13, 30, 36, 50, 75), repeat=2):
            a = rng.standard_normal(ma)
            b = rng.standard_normal(nb)
            pairs['1d'].append((a, b))

        for n_image in [256, 512, 1024]:
            for n_kernel in [3, 5, 7]:
                x = rng.standard_normal((n_image, n_image))
                h = rng.standard_normal((n_kernel, n_kernel))
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
        rng = np.random.default_rng(1234)
        # sample a bunch of pairs of 2d arrays
        pairs = []
        for nfilt in [8, ]:
            for n in [32, 128, 512, 2048]:
                h = rng.standard_normal(nfilt)
                x = rng.standard_normal(n)
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
        rng = np.random.default_rng(1234)
        # sample a bunch of pairs of 2d arrays
        pairs = []
        for nfilt in [8, ]:
            for n in [32, 128, 512]:
                h = rng.standard_normal(nfilt)
                x = rng.standard_normal((n, n))
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


class IdentifyRidgeLines(Benchmark):
    """Benchmark _identify_ridge_lines performance directly"""

    param_names = ['matrix_type', 'size']
    params = [
        ['synthetic_sparse', 'synthetic_dense', 'cwt_real'],
        ['small', 'medium', 'large']
    ]

    def setup(self, matrix_type, size):
        from scipy.signal._wavelets import _cwt, _ricker

        # Setup test matrices based on parameters
        sizes = {
            'small': (50, 500),
            'medium': (100, 1000),
            'large': (200, 5000)
        }
        n_rows, n_cols = sizes[size]

        if matrix_type == 'synthetic_sparse':
            rng = np.random.RandomState(42)
            self.matr = rng.random((n_rows, n_cols))
            mask = rng.random((n_rows, n_cols)) > 0.01
            self.matr[mask] = 0
        elif matrix_type == 'synthetic_dense':
            rng = np.random.RandomState(42)
            self.matr = rng.random((n_rows, n_cols))
            mask = rng.random((n_rows, n_cols)) > 0.1
            self.matr[mask] = 0
        else:  # cwt_real
            # Generate CWT matrix
            signal_length = n_cols
            n_widths = n_rows
            # Create signal with Gaussian peaks
            rng = np.random.RandomState(42)
            sig = np.sin(np.linspace(0, 10*np.pi, signal_length))
            sig += rng.normal(0, 0.1, signal_length)
            widths = np.arange(1, n_widths + 1)
            self.matr = _cwt(sig, _ricker, widths)

        self.max_distances = np.full(n_rows, 3)
        self.gap_thresh = 2

    def time_identify_ridge_lines(self, matrix_type, size):
        from scipy.signal._peak_finding import _identify_ridge_lines
        _identify_ridge_lines(self.matr, self.max_distances, self.gap_thresh)

    def peakmem_identify_ridge_lines(self, matrix_type, size):
        from scipy.signal._peak_finding import _identify_ridge_lines
        _identify_ridge_lines(self.matr, self.max_distances, self.gap_thresh)


class FindPeaksCWT(Benchmark):
    """Benchmark find_peaks_cwt performance (end-to-end with _identify_ridge_lines)"""

    param_names = ['signal_length', 'n_widths']
    params = [
        [500, 1000, 2000, 5000],  # Signal lengths
        [30, 50, 100]  # Number of widths
    ]

    def setup(self, signal_length, n_widths):
        # Create realistic signal with multiple Gaussian peaks
        rng = np.random.RandomState(42)
        t = np.linspace(0, 1, signal_length)

        # Add several Gaussian peaks at different positions
        self.signal = np.zeros(signal_length)
        peak_positions = [0.2, 0.4, 0.6, 0.8]
        for pos in peak_positions:
            self.signal += np.exp(-((t - pos) / 0.05) ** 2)

        # Add noise
        self.signal += rng.normal(0, 0.1, signal_length)

        # Widths for CWT
        self.widths = np.arange(1, n_widths + 1)

    def time_find_peaks_cwt(self, signal_length, n_widths):
        signal.find_peaks_cwt(self.signal, self.widths)

    def peakmem_find_peaks_cwt(self, signal_length, n_widths):
        signal.find_peaks_cwt(self.signal, self.widths)
