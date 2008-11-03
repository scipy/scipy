#this program corresponds to special.py

import numpy as np
from numpy.testing import *

import scipy.signal as signal

class TestChirp(TestCase):
    def test_log_chirp_at_zero(self):
        assert_almost_equal(signal.waveforms.chirp(t=0, method='log'),
                            1.0)

if __name__ == "__main__":
    run_module_suite()
