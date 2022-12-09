import numpy as np
from numpy.testing import assert_allclose
from scipy.signal import gausspulse

from scipy.signal import wigner_ville


class TestWigner:
    def test_wigner_ville_shapes(self):
        t, f, wv = wigner_ville(np.sin(np.linspace(0, 100, 1024)),
                                       resolution=4, window_size=128)
        assert len(t) == 256
        assert len(f) == 128
        assert wv.shape == (128, 256)

    def test_wigner_ville_tc_fc(self):
        fs = 100
        fc = 4.2
        tc = 512 / fs
        t = np.arange(0, 1024) / fs
        x_re, x_im = gausspulse(t - tc, fc, bw=0.2, retquad=True)
        t, f, wv = wigner_ville(x_re + 1j * x_im, fs)
        max_idx = np.unravel_index(np.argmax(wv), wv.shape)
        assert_allclose(fc, f[max_idx[0]], atol=f[1] - f[0])
        assert_allclose(tc, t[max_idx[1]], atol=t[1] - t[0])
