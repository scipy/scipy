import warnings

import numpy as np
from numpy.testing import TestCase, assert_array_almost_equal

from scipy.signal import tf2zpk, bessel, BadCoefficients, kaiserord, firwin, freqz


class TestTf2zpk(TestCase):
    def test_simple(self):
        z_r = np.array([0.5, -0.5])
        p_r = np.array([1.j / np.sqrt(2), -1.j / np.sqrt(2)])
        # Sort the zeros/poles so that we don't fail the test if the order
        # changes
        z_r.sort()
        p_r.sort()
        b = np.poly(z_r)
        a = np.poly(p_r)

        z, p, k = tf2zpk(b, a)
        z.sort()
        p.sort()
        assert_array_almost_equal(z, z_r)
        assert_array_almost_equal(p, p_r)

    def test_bad_filter(self):
        """Regression test for #651: better handling of badly conditioned
        filter coefficients."""
        warnings.simplefilter("error", BadCoefficients)
        try:
            try:
                b, a = bessel(20, 0.1)
                z, p, k = tf2zpk(b, a)
                raise AssertionError("tf2zpk did not warn about bad "\
                                     "coefficients")
            except BadCoefficients:
                pass
        finally:
            warnings.simplefilter("always", BadCoefficients)


class TestFirWin(TestCase):
    
    def test_lowpass(self):
        width = 0.04
        ntaps, beta = kaiserord(120, width)
        taps = firwin(ntaps, cutoff=0.5, window=('kaiser', beta))
        freq_samples = np.array([0.0, 0.25, 0.5-width/2, 0.5+width/2, 0.75, 1.0])
        freqs, response = freqz(taps, worN=np.pi*freq_samples)
        assert_array_almost_equal(np.abs(response),
                                    [1.0, 1.0, 1.0, 0.0, 0.0, 0.0], decimal=5)
