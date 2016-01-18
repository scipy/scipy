from __future__ import division, absolute_import, print_function

import numpy as np

try:
    from scipy.signal import lfilter, firwin
except ImportError:
    pass

from .common import Benchmark


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
