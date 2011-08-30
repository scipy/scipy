"""Unit tests for Lomb Scargle routines.
"""

import numpy as np
from numpy.testing import dec, assert_raises, assert_equal, \
                          assert_almost_equal, assert_array_equal, \
                          assert_array_almost_equal, assert_approx_equal, \
                          assert_, run_module_suite

from scipy.signal.spectral import lombscargle


class TestLombscargle:
    def test_frequency(self):
        """Test if frequency location of peak corresponds to frequency of
        generated input signal.
        """

        # Input parameters
        ampl = 2.
        w = 1.
        phi = 0.5 * np.pi
        nin = 100
        nout = 1000
        p = 0.7 # Fraction of points to select

        # Randomly select a fraction of an array with timesteps
        np.random.seed(2353425)
        r = np.random.rand(nin)
        t = np.linspace(0.01*np.pi, 10.*np.pi, nin)[r >= p]

        # Plot a sine wave for the selected times
        x = ampl * np.sin(w*t + phi)

        # Define the array of frequencies for which to compute the periodogram
        f = np.linspace(0.01, 10., nout)

        # Calculate Lomb-Scargle periodogram
        P = lombscargle(t, x, f)

        # Check if difference between found frequency maximum and input
        # frequency is less than accuracy
        delta = f[1] - f[0]
        assert_(w - f[np.argmax(P)] < (delta/2.))

    def test_amplitude(self):
        """Test if height of peak in normalized Lomb-Scargle periodogram
        corresponds to amplitude of the generated input signal.
        """

        # Input parameters
        ampl = 2.
        w = 1.
        phi = 0.5 * np.pi
        nin = 100
        nout = 1000
        p = 0.7 # Fraction of points to select

        # Randomly select a fraction of an array with timesteps
        np.random.seed(2353425)
        r = np.random.rand(nin)
        t = np.linspace(0.01*np.pi, 10.*np.pi, nin)[r >= p]

        # Plot a sine wave for the selected times
        x = ampl * np.sin(w*t + phi)

        # Define the array of frequencies for which to compute the periodogram
        f = np.linspace(0.01, 10., nout)

        # Calculate Lomb-Scargle periodogram
        pgram = lombscargle(t, x, f)

        # Normalize
        pgram = np.sqrt(4 * pgram / t.shape[0])

        # Check if difference between found frequency maximum and input
        # frequency is less than accuracy
        assert_approx_equal(np.max(pgram), ampl, significant=2)

    def test_wrong_shape(self):
        t = np.linspace(0, 1, 1)
        x = np.linspace(0, 1, 2)
        f = np.linspace(0, 1, 3)
        assert_raises(ValueError, lombscargle, t, x, f)

    def test_zero_division(self):
        t = np.zeros(1)
        x = np.zeros(1)
        f = np.zeros(1)
        assert_raises(ZeroDivisionError, lombscargle, t, x, f)


if __name__ == "__main__":
    run_module_suite()
