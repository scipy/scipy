
import numpy as np
from numpy.testing import TestCase, run_module_suite, assert_raises, \
        assert_almost_equal, assert_array_almost_equal, assert_equal, \
        assert_
from scipy.special import sinc

from scipy.signal import kaiser_beta, kaiser_atten, kaiserord, \
        firwin, firwin2, freqz, remez


def test_kaiser_beta():
    b = kaiser_beta(58.7)
    assert_almost_equal(b, 0.1102 * 50.0)
    b = kaiser_beta(22.0)
    assert_almost_equal(b, 0.5842 + 0.07886)
    b = kaiser_beta(21.0)
    assert_equal(b, 0.0)
    b = kaiser_beta(10.0)
    assert_equal(b, 0.0)


def test_kaiser_atten():
    a = kaiser_atten(1, 1.0)
    assert_equal(a, 7.95)
    a = kaiser_atten(2, 1/np.pi)
    assert_equal(a, 2.285 + 7.95)


def test_kaiserord():
    assert_raises(ValueError, kaiserord, 1.0, 1.0)
    numtaps, beta = kaiserord(2.285 + 7.95 - 0.001, 1/np.pi)
    assert_equal((numtaps, beta), (2, 0.0))


class TestFirwin(TestCase):

    def check_response(self, h, expected_response, tol=.05):
        N = len(h)
        alpha = 0.5 * (N-1)
        m = np.arange(0,N) - alpha   # time indices of taps
        for freq, expected in expected_response:
            actual = abs(np.sum(h*np.exp(-1.j*np.pi*m*freq)))
            mse = abs(actual-expected)**2
            self.assertTrue(mse < tol, 'response not as expected, mse=%g > %g'\
               %(mse, tol))

    def test_response(self):
        N = 51
        f = .5
        # increase length just to try even/odd
        h = firwin(N, f) # low-pass from 0 to f
        self.check_response(h, [(.25,1), (.75,0)])

        h = firwin(N+1, f, window='nuttall') # specific window
        self.check_response(h, [(.25,1), (.75,0)])

        h = firwin(N+2, f, pass_zero=False) # stop from 0 to f --> high-pass
        self.check_response(h, [(.25,0), (.75,1)])

        f1, f2, f3, f4 = .2, .4, .6, .8
        h = firwin(N+3, [f1, f2], pass_zero=False) # band-pass filter
        self.check_response(h, [(.1,0), (.3,1), (.5,0)])

        h = firwin(N+4, [f1, f2]) # band-stop filter
        self.check_response(h, [(.1,1), (.3,0), (.5,1)])

        h = firwin(N+5, [f1, f2, f3, f4], pass_zero=False, scale=False)
        self.check_response(h, [(.1,0), (.3,1), (.5,0), (.7,1), (.9,0)])

        h = firwin(N+6, [f1, f2, f3, f4])  # multiband filter
        self.check_response(h, [(.1,1), (.3,0), (.5,1), (.7,0), (.9,1)])

        h = firwin(N+7, 0.1, width=.03) # low-pass
        self.check_response(h, [(.05,1), (.75,0)])

        h = firwin(N+8, 0.1, pass_zero=False) # high-pass
        self.check_response(h, [(.05,0), (.75,1)])

    def mse(self, h, bands):
        """Compute mean squared error versus ideal response across frequency
        band.
          h -- coefficients
          bands -- list of (left, right) tuples relative to 1==Nyquist of
            passbands
        """
        w, H = freqz(h, worN=1024)
        f = w/np.pi
        passIndicator = np.zeros(len(w), bool)
        for left, right in bands:
            passIndicator |= (f>=left) & (f<right)
        Hideal = np.where(passIndicator, 1, 0)
        mse = np.mean(abs(abs(H)-Hideal)**2)
        return mse

    def test_scaling(self):
        """
        For one lowpass, bandpass, and highpass example filter, this test
        checks two things:
          - the mean squared error over the frequency domain of the unscaled
            filter is smaller than the scaled filter (true for rectangular
            window)
          - the response of the scaled filter is exactly unity at the center
            of the first passband
        """
        N = 11
        cases = [
            ([.5],      True,   (0, 1)),
            ([0.2, .6], False,  (.4, 1)),
            ([.5],      False,  (1, 1)),
        ]
        for cutoff, pass_zero, expected_response in cases:
            h = firwin(N, cutoff, scale=False, pass_zero=pass_zero, window='ones')
            hs = firwin(N, cutoff, scale=True, pass_zero=pass_zero, window='ones')
            if len(cutoff) == 1:
                if pass_zero:
                    cutoff = [0] + cutoff
                else:
                    cutoff = cutoff + [1]
            self.assertTrue(self.mse(h, [cutoff]) < self.mse(hs, [cutoff]),
                'least squares violation')
            self.check_response(hs, [expected_response], 1e-12)


class TestFirWinMore(TestCase):
    """Different author, different style, different tests..."""

    def test_lowpass(self):
        width = 0.04
        ntaps, beta = kaiserord(120, width)
        taps = firwin(ntaps, cutoff=0.5, window=('kaiser', beta), scale=False)

        # Check the symmetry of taps.
        assert_array_almost_equal(taps[:ntaps//2], taps[ntaps:ntaps-ntaps//2-1:-1])

        # Check the gain at a few samples where we know it should be approximately 0 or 1.
        freq_samples = np.array([0.0, 0.25, 0.5-width/2, 0.5+width/2, 0.75, 1.0])
        freqs, response = freqz(taps, worN=np.pi*freq_samples)
        assert_array_almost_equal(np.abs(response),
                                    [1.0, 1.0, 1.0, 0.0, 0.0, 0.0], decimal=5)

    def test_highpass(self):
        width = 0.04
        ntaps, beta = kaiserord(120, width)

        # Ensure that ntaps is odd.
        ntaps |= 1

        taps = firwin(ntaps, cutoff=0.5, window=('kaiser', beta),
                        pass_zero=False, scale=False)

        # Check the symmetry of taps.
        assert_array_almost_equal(taps[:ntaps//2], taps[ntaps:ntaps-ntaps//2-1:-1])

        # Check the gain at a few samples where we know it should be approximately 0 or 1.
        freq_samples = np.array([0.0, 0.25, 0.5-width/2, 0.5+width/2, 0.75, 1.0])
        freqs, response = freqz(taps, worN=np.pi*freq_samples)
        assert_array_almost_equal(np.abs(response),
                                    [0.0, 0.0, 0.0, 1.0, 1.0, 1.0], decimal=5)

    def test_bandpass(self):
        width = 0.04
        ntaps, beta = kaiserord(120, width)
        taps = firwin(ntaps, cutoff=[0.3, 0.7], window=('kaiser', beta),
                        pass_zero=False, scale=False)

        # Check the symmetry of taps.
        assert_array_almost_equal(taps[:ntaps//2], taps[ntaps:ntaps-ntaps//2-1:-1])

        # Check the gain at a few samples where we know it should be approximately 0 or 1.
        freq_samples = np.array([0.0, 0.2, 0.3-width/2, 0.3+width/2, 0.5,
                                0.7-width/2, 0.7+width/2, 0.8, 1.0])
        freqs, response = freqz(taps, worN=np.pi*freq_samples)
        assert_array_almost_equal(np.abs(response),
                [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], decimal=5)

    def test_multi(self):
        width = 0.04
        ntaps, beta = kaiserord(120, width)
        taps = firwin(ntaps, cutoff=[0.2, 0.5, 0.8], window=('kaiser', beta),
                        pass_zero=True, scale=False)

        # Check the symmetry of taps.
        assert_array_almost_equal(taps[:ntaps//2], taps[ntaps:ntaps-ntaps//2-1:-1])

        # Check the gain at a few samples where we know it should be approximately 0 or 1.
        freq_samples = np.array([0.0, 0.1, 0.2-width/2, 0.2+width/2, 0.35,
                                0.5-width/2, 0.5+width/2, 0.65,
                                0.8-width/2, 0.8+width/2, 0.9, 1.0])
        freqs, response = freqz(taps, worN=np.pi*freq_samples)
        assert_array_almost_equal(np.abs(response),
                [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
                decimal=5)

    def test_nyq(self):
        """Test the nyq keyword."""
        nyquist = 1000
        width = 40.0
        relative_width = width/nyquist
        ntaps, beta = kaiserord(120, relative_width)
        taps = firwin(ntaps, cutoff=[300, 700], window=('kaiser', beta),
                        pass_zero=False, scale=False, nyq=nyquist)

        # Check the symmetry of taps.
        assert_array_almost_equal(taps[:ntaps//2], taps[ntaps:ntaps-ntaps//2-1:-1])

        # Check the gain at a few samples where we know it should be approximately 0 or 1.
        freq_samples = np.array([0.0, 200, 300-width/2, 300+width/2, 500,
                                700-width/2, 700+width/2, 800, 1000])
        freqs, response = freqz(taps, worN=np.pi*freq_samples/nyquist)
        assert_array_almost_equal(np.abs(response),
                [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], decimal=5)

    def test_bad_cutoff(self):
        """Test that invalid cutoff argument raises ValueError."""
        # cutoff values must be greater than 0 and less than 1.
        assert_raises(ValueError, firwin, 99, -0.5)
        assert_raises(ValueError, firwin, 99, 1.5)
        # Don't allow 0 or 1 in cutoff.
        assert_raises(ValueError, firwin, 99, [0, 0.5])
        assert_raises(ValueError, firwin, 99, [0.5, 1])
        # cutoff values must be strictly increasing.
        assert_raises(ValueError, firwin, 99, [0.1, 0.5, 0.2])
        assert_raises(ValueError, firwin, 99, [0.1, 0.5, 0.5])
        # Must have at least one cutoff value.
        assert_raises(ValueError, firwin, 99, [])
        # 2D array not allowed.
        assert_raises(ValueError, firwin, 99, [[0.1, 0.2],[0.3, 0.4]])
        # cutoff values must be less than nyq.
        assert_raises(ValueError, firwin, 99, 50.0, nyq=40)
        assert_raises(ValueError, firwin, 99, [10, 20, 30], nyq=25)

    def test_even_highpass_raises_value_error(self):
        """Test that attempt to create a highpass filter with an even number
        of taps raises a ValueError exception."""
        assert_raises(ValueError, firwin, 40, 0.5, pass_zero=False)
        assert_raises(ValueError, firwin, 40, [.25, 0.5])




class TestFirwin2(TestCase):

    def test_invalid_args(self):
        # `freq` and `gain` have different lengths.
        assert_raises(ValueError, firwin2, 50, [0, 0.5, 1], [0.0, 1.0])
        # `nfreqs` is less than `ntaps`.
        assert_raises(ValueError, firwin2, 50, [0, 0.5, 1], [0.0, 1.0, 1.0], nfreqs=33)
        # Decreasing value in `freq`
        assert_raises(ValueError, firwin2, 50, [0, 0.5, 0.4, 1.0], [0, .25, .5, 1.0])
        # Value in `freq` repeated more than once.
        assert_raises(ValueError, firwin2, 50, [  0,   .1,   .1,   .1, 1.0],
                                               [0.0,  0.5, 0.75,  1.0, 1.0])
        # `freq` does not start at 0.0.
        assert_raises(ValueError, firwin2, 50, [0.5, 1.0], [0.0, 1.0])

        # Type II filter, but the gain at nyquist rate is not zero.
        assert_raises(ValueError, firwin2, 16, [0.0, 0.5, 1.0], [0.0, 1.0, 1.0])

        # Type III filter, but the gains at nyquist and zero rate are not zero.
        assert_raises(ValueError, firwin2, 17, [0.0, 0.5, 1.0], [0.0, 1.0, 1.0],
                      antisymmetric=True)
        assert_raises(ValueError, firwin2, 17, [0.0, 0.5, 1.0], [1.0, 1.0, 0.0],
                      antisymmetric=True)
        assert_raises(ValueError, firwin2, 17, [0.0, 0.5, 1.0], [1.0, 1.0, 1.0],
                      antisymmetric=True)

        # Type VI filter, but the gain at zero rate is not zero.
        assert_raises(ValueError, firwin2, 16, [0.0, 0.5, 1.0], [1.0, 1.0, 0.0],
                      antisymmetric=True)

    def test01(self):
        width = 0.04
        beta = 12.0
        ntaps = 400
        # Filter is 1 from w=0 to w=0.5, then decreases linearly from 1 to 0 as w
        # increases from w=0.5 to w=1  (w=1 is the Nyquist frequency).
        freq = [0.0, 0.5, 1.0]
        gain = [1.0, 1.0, 0.0]
        taps = firwin2(ntaps, freq, gain, window=('kaiser', beta))
        freq_samples = np.array([0.0, 0.25, 0.5-width/2, 0.5+width/2,
                                                        0.75, 1.0-width/2])
        freqs, response = freqz(taps, worN=np.pi*freq_samples)
        assert_array_almost_equal(np.abs(response),
                        [1.0, 1.0, 1.0, 1.0-width, 0.5, width], decimal=5)

    def test02(self):
        width = 0.04
        beta = 12.0
        # ntaps must be odd for positive gain at Nyquist.
        ntaps = 401
        # An ideal highpass filter.
        freq = [0.0, 0.5, 0.5, 1.0]
        gain = [0.0, 0.0, 1.0, 1.0]
        taps = firwin2(ntaps, freq, gain, window=('kaiser', beta))
        freq_samples = np.array([0.0, 0.25, 0.5-width, 0.5+width, 0.75, 1.0])
        freqs, response = freqz(taps, worN=np.pi*freq_samples)
        assert_array_almost_equal(np.abs(response),
                                [0.0, 0.0, 0.0, 1.0, 1.0, 1.0], decimal=5)

    def test03(self):
        width = 0.02
        ntaps, beta = kaiserord(120, width)
        # ntaps must be odd for positive gain at Nyquist.
        ntaps = int(ntaps) | 1
        freq = [0.0, 0.4, 0.4, 0.5, 0.5, 1.0]
        gain = [1.0, 1.0, 0.0, 0.0, 1.0, 1.0]
        taps = firwin2(ntaps, freq, gain, window=('kaiser', beta))
        freq_samples = np.array([0.0, 0.4-width, 0.4+width, 0.45,
                                    0.5-width, 0.5+width, 0.75, 1.0])
        freqs, response = freqz(taps, worN=np.pi*freq_samples)
        assert_array_almost_equal(np.abs(response),
                    [1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0], decimal=5)

    def test04(self):
        """Test firwin2 when window=None."""
        ntaps = 5
        # Ideal lowpass: gain is 1 on [0,0.5], and 0 on [0.5, 1.0]
        freq = [0.0, 0.5, 0.5, 1.0]
        gain = [1.0, 1.0, 0.0, 0.0]
        taps = firwin2(ntaps, freq, gain, window=None, nfreqs=8193)
        alpha = 0.5 * (ntaps - 1)
        m = np.arange(0, ntaps) - alpha
        h = 0.5 * sinc(0.5 * m)
        assert_array_almost_equal(h, taps)

    def test05(self):
        """Test firwin2 for calculating Type IV filters"""
        ntaps = 1500

        freq = [0.0, 1.0]
        gain = [0.0, 1.0]
        taps = firwin2(ntaps, freq, gain, window=None, antisymmetric=True)
        assert_array_almost_equal(taps[: ntaps // 2], -taps[ntaps // 2:][::-1])


        freqs, response = freqz(taps, worN=2048)
        assert_array_almost_equal(abs(response), freqs / np.pi, decimal=4)

    def test06(self):
        """Test firwin2 for calculating Type III filters"""
        ntaps = 1501

        freq = [0.0, 0.5, 0.55,  1.0]
        gain = [0.0, 0.5, 0.0,   0.0]
        taps = firwin2(ntaps, freq, gain, window=None, antisymmetric=True)
        assert_equal(taps[ntaps // 2], 0.0)
        assert_array_almost_equal(taps[: ntaps // 2], -taps[ntaps // 2 + 1:][::-1])

        freqs, response1 = freqz(taps, worN=2048)
        response2        = np.interp(freqs / np.pi, freq, gain)
        assert_array_almost_equal(abs(response1), response2, decimal=3)


    def test_nyq(self):
        taps1 = firwin2(80, [0.0, 0.5, 1.0], [1.0, 1.0, 0.0])
        taps2 = firwin2(80, [0.0, 30.0, 60.0], [1.0, 1.0, 0.0], nyq=60.0)
        assert_array_almost_equal(taps1, taps2)


class TestRemez(TestCase):

    def test_bad_args(self):
        assert_raises(ValueError, remez, 11, [0.1, 0.4], [1], type='pooka')

    def test_hilbert(self):
        N = 11 # number of taps in the filter
        a = 0.1 # width of the transition band

        # design an unity gain hilbert bandpass filter from w to 0.5-w
        h = remez(11, [ a, 0.5-a ], [ 1 ], type='hilbert')

        # make sure the filter has correct # of taps
        assert_(len(h) == N, "Number of Taps")

        # make sure it is type III (anti-symmetric tap coefficients)
        assert_array_almost_equal(h[:(N-1)//2], -h[:-(N-1)//2-1:-1])

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


if __name__ == "__main__":
    run_module_suite()
