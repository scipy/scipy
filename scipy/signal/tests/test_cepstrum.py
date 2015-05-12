import numpy as np

from scipy.signal import (complex_cepstrum, inverse_complex_cepstrum, 
                          real_cepstrum, minimum_phase, sawtooth)

from numpy.testing import (assert_array_almost_equal)


def test_complex_cepstrum():
    """The period of a periodic harmonic will show up as a peak in a 
    complex cepstrum.
    """
    duration = 5.0
    fs = 8000.0
    samples = int(fs*duration)
    t = np.arange(samples) / fs
    fundamental = 100.0
    signal = sawtooth(2.*np.pi*fundamental*t)
    ceps, _ = complex_cepstrum(signal)
    assert(fundamental == 1.0/t[ceps.argmax()])

def test_inverse_complex_cepstrum():
    """Applying the complex cepstrum and then the inverse complex cepstrum 
    should result in the original sequence.
    """
    x = np.arange(10)
    ceps, ndelay = complex_cepstrum(x)
    y = inverse_complex_cepstrum(ceps, ndelay)
    assert_array_almost_equal(x, y)
