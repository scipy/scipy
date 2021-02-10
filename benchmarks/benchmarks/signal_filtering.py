import numpy as np
from numpy.testing import assert_allclose
import timeit
from concurrent.futures import ThreadPoolExecutor, as_completed, wait

from .common import Benchmark, safe_import

with safe_import():
    from scipy.signal import (lfilter, firwin, decimate, butter, sosfilt,
                              medfilt2d)


class Decimate(Benchmark):
    param_names = ['q', 'ftype', 'zero_phase']
    params = [
        [2, 10, 30],
        ['iir', 'fir'],
        [True, False]
    ]

    def setup(self, q, ftype, zero_phase):
        np.random.seed(123456)
        sample_rate = 10000.
        t = np.arange(int(1e6), dtype=np.float64) / sample_rate
        self.sig = np.sin(2*np.pi*500*t) + 0.3 * np.sin(2*np.pi*4e3*t)

    def time_decimate(self, q, ftype, zero_phase):
        decimate(self.sig, q, ftype=ftype, zero_phase=zero_phase)


class Lfilter(Benchmark):
    param_names = ['n_samples', 'numtaps']
    params = [
        [1e3, 50e3, 1e6],
        [9, 23, 51]
    ]

    def setup(self, n_samples, numtaps):
        np.random.seed(125678)
        sample_rate = 25000.
        t = np.arange(n_samples, dtype=np.float64) / sample_rate
        nyq_rate = sample_rate / 2.
        cutoff_hz = 3000.0
        self.sig = np.sin(2*np.pi*500*t) + 0.3 * np.sin(2*np.pi*11e3*t)
        self.coeff = firwin(numtaps, cutoff_hz/nyq_rate)

    def time_lfilter(self, n_samples, numtaps):
        lfilter(self.coeff, 1.0, self.sig)

class ParallelSosfilt(Benchmark):
    timeout = 100
    timer = timeit.default_timer

    param_names = ['n_samples', 'threads']
    params = [
        [1e3, 10e3],
        [1, 2, 4]
    ]

    def setup(self, n_samples, threads):
        self.filt = butter(8, 8e-6, "lowpass", output="sos")
        self.data = np.arange(int(n_samples) * 3000).reshape(int(n_samples), 3000)
        self.chunks = np.array_split(self.data, threads)

    def time_sosfilt(self, n_samples, threads):
        with ThreadPoolExecutor(max_workers=threads) as pool:
            futures = []
            for i in range(threads):
                futures.append(pool.submit(sosfilt, self.filt, self.chunks[i]))

            wait(futures)


class Sosfilt(Benchmark):
    param_names = ['n_samples', 'order']
    params = [
        [1000, 1000000],
        [6, 20]
    ]

    def setup(self, n_samples, order):
        self.sos = butter(order, [0.1575, 0.1625], 'band', output='sos')
        self.y = np.random.RandomState(0).randn(n_samples)

    def time_sosfilt_basic(self, n_samples, order):
        sosfilt(self.sos, self.y)


class MedFilt2D(Benchmark):
    def setup_cache(self):
        np.random.seed(8176)
        data = np.random.randn(2000, 4000)
        expected = medfilt2d(data, 5)
        np.save("medfilt2d_data.npy", data)
        np.save("medfilt2d_expected.npy", expected)

    def setup(self):
        self.data = np.load("medfilt2d_data.npy")
        self.expected = np.load("medfilt2d_expected.npy")

    def teardown(self):
        assert_allclose(self.output, self.expected)

    def _multithreaded(self):
        # Lets take four overlapping chunks over the last axis. First pair of
        # values is range to select from input, second pair range of the result
        # to use, last pair the range of the output array to store it in.
        chunks = (
            (0, 1002, None, -2, 0, 1000),
            (998, 2002, 2, -2, 1000, 2000),
            (1998, 3002, 2, -2, 2000, 3000),
            (2998, 4000, 2, None, 3000, 4000)
        )

        self.output = np.empty(self.data.shape, dtype=self.data.dtype)

        def apply(chunk):
            return medfilt2d(self.data[:, chunk[0]:chunk[1]], 5)

        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = {executor.submit(apply, chk): chk for chk in chunks}
            for future in as_completed(futures):
                chk = futures[future]
                result = future.result()
                self.output[:, chk[4]:chk[5]] = result[:, chk[2]:chk[3]]

    def time_singlethreaded(self):
        self.output = medfilt2d(self.data, 5)

    def time_multithreaded(self):
        self._multithreaded()

    def peakmem_singlethreaded(self):
        self.output = medfilt2d(self.data, 5)

    def peakmem_multithreaded(self):
        self._multithreaded()
