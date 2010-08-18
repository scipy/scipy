import warnings

import numpy as np
from numpy.testing import TestCase, assert_array_almost_equal, assert_

from scipy.signal import tf2zpk, bessel, BadCoefficients, kaiserord, firwin, freqz, remez


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

class TestRemez(TestCase):

    def test_hilbert(self):
        N = 11 # number of taps in the filter
        a = 0.1 # width of the transition band

        # design an unity gain hilbert bandpass filter from w to 0.5-w
        h = remez(11, [ a, 0.5-a ], [ 1 ], type='hilbert')

        # make sure the filter has correct # of taps
        assert_(len(h) == N, "Number of Taps")

        # make sure it is type III (anti-symmtric tap coefficients)
        assert_array_almost_equal(h[:(N-1)/2], -h[:-(N-1)/2-1:-1])

        # Since the requested response is symmetric, all even coeffcients
        # should be zero (or in this case really small)
        assert_((abs(h[1::2]) < 1e-15).all(), "Even Coefficients Equal Zero")

        # now check the frequency response
        w, H = freqz(h, 1)
        f = w/2/np.pi
        Hmag = abs(H)

        # should have a zero at 0 and pi (in this case close to zero)
        assert_((Hmag[ [0,-1] ] < 0.02).all(), "Zero at zero and pi")

        # check that the pass band is close to unity
        idx = (f > a) * (f < 0.5-a)
        assert_((abs(Hmag[idx] - 1) < 0.015).all(), "Pass Band Close To Unity")
