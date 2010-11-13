
import numpy as np
from numpy.testing import TestCase, run_module_suite, assert_raises, \
        assert_array_almost_equal

from scipy.signal import firwin, kaiserord, freqz


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
        assert_array_almost_equal(taps[:ntaps/2], taps[ntaps:ntaps-ntaps/2-1:-1])

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
        assert_array_almost_equal(taps[:ntaps/2], taps[ntaps:ntaps-ntaps/2-1:-1])

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
        assert_array_almost_equal(taps[:ntaps/2], taps[ntaps:ntaps-ntaps/2-1:-1])

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
        assert_array_almost_equal(taps[:ntaps/2], taps[ntaps:ntaps-ntaps/2-1:-1])

        # Check the gain at a few samples where we know it should be approximately 0 or 1.
        freq_samples = np.array([0.0, 0.1, 0.2-width/2, 0.2+width/2, 0.35,
                                0.5-width/2, 0.5+width/2, 0.65,
                                0.8-width/2, 0.8+width/2, 0.9, 1.0])
        freqs, response = freqz(taps, worN=np.pi*freq_samples)
        assert_array_almost_equal(np.abs(response),
                [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
                decimal=5)

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


    def test_even_highpass_raises_value_error(self):
        """Test that attempt to create a highpass filter with an even number
        of taps raises a ValueError exception."""
        assert_raises(ValueError, firwin, 40, 0.5, pass_zero=False)
        assert_raises(ValueError, firwin, 40, [.25, 0.5])


if __name__ == "__main__":
    run_module_suite()
