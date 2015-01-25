from __future__ import division, print_function, absolute_import

import warnings
import numpy as np
from numpy import array
from numpy.testing import (assert_array_almost_equal, assert_array_equal,
                           run_module_suite, assert_raises)
from scipy import signal


window_funcs = [
    ('boxcar', ()),
    ('triang', ()),
    ('parzen', ()),
    ('bohman', ()),
    ('blackman', ()),
    ('nuttall', ()),
    ('blackmanharris', ()),
    ('flattop', ()),
    ('bartlett', ()),
    ('hanning', ()),
    ('barthann', ()),
    ('hamming', ()),
    ('kaiser', (1,)),
    ('gaussian', (0.5,)),
    ('general_gaussian', (1.5, 2)),
    ('chebwin', (1,)),
    ('slepian', (2,)),
    ('cosine', ()),
    ('hann', ()),
    ]


cheb_odd_true = array([0.200938, 0.107729, 0.134941, 0.165348,
                       0.198891, 0.235450, 0.274846, 0.316836,
                       0.361119, 0.407338, 0.455079, 0.503883,
                       0.553248, 0.602637, 0.651489, 0.699227,
                       0.745266, 0.789028, 0.829947, 0.867485,
                       0.901138, 0.930448, 0.955010, 0.974482,
                       0.988591, 0.997138, 1.000000, 0.997138,
                       0.988591, 0.974482, 0.955010, 0.930448,
                       0.901138, 0.867485, 0.829947, 0.789028,
                       0.745266, 0.699227, 0.651489, 0.602637,
                       0.553248, 0.503883, 0.455079, 0.407338,
                       0.361119, 0.316836, 0.274846, 0.235450,
                       0.198891, 0.165348, 0.134941, 0.107729,
                       0.200938])

cheb_even_true = array([0.203894, 0.107279, 0.133904,
                        0.163608, 0.196338, 0.231986,
                        0.270385, 0.311313, 0.354493,
                        0.399594, 0.446233, 0.493983,
                        0.542378, 0.590916, 0.639071,
                        0.686302, 0.732055, 0.775783,
                        0.816944, 0.855021, 0.889525,
                        0.920006, 0.946060, 0.967339,
                        0.983557, 0.994494, 1.000000,
                        1.000000, 0.994494, 0.983557,
                        0.967339, 0.946060, 0.920006,
                        0.889525, 0.855021, 0.816944,
                        0.775783, 0.732055, 0.686302,
                        0.639071, 0.590916, 0.542378,
                        0.493983, 0.446233, 0.399594,
                        0.354493, 0.311313, 0.270385,
                        0.231986, 0.196338, 0.163608,
                        0.133904, 0.107279, 0.203894])


class TestChebWin(object):

    def test_cheb_odd_high_attenuation(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            cheb_odd = signal.chebwin(53, at=-40)
        assert_array_almost_equal(cheb_odd, cheb_odd_true, decimal=4)

    def test_cheb_even_high_attenuation(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            cheb_even = signal.chebwin(54, at=-40)
        assert_array_almost_equal(cheb_even, cheb_even_true, decimal=4)

    def test_cheb_odd_low_attenuation(self):
        cheb_odd_low_at_true = array([1.000000, 0.519052, 0.586405,
                                      0.610151, 0.586405, 0.519052,
                                      1.000000])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            cheb_odd = signal.chebwin(7, at=-10)
        assert_array_almost_equal(cheb_odd, cheb_odd_low_at_true, decimal=4)

    def test_cheb_even_low_attenuation(self):
        cheb_even_low_at_true = array([1.000000, 0.451924, 0.51027,
                                       0.541338, 0.541338, 0.51027,
                                       0.451924, 1.000000])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            cheb_even = signal.chebwin(8, at=-10)
        assert_array_almost_equal(cheb_even, cheb_even_low_at_true, decimal=4)


class TestGetWindow(object):

    def test_boxcar(self):
        w = signal.get_window('boxcar', 12)
        assert_array_equal(w, np.ones_like(w))

    def test_cheb_odd(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            w = signal.get_window(('chebwin', -40), 53, fftbins=False)
        assert_array_almost_equal(w, cheb_odd_true, decimal=4)

    def test_cheb_even(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            w = signal.get_window(('chebwin', -40), 54, fftbins=False)
        assert_array_almost_equal(w, cheb_even_true, decimal=4)

    def test_array_as_window(self):
        # github issue 3603
        osfactor = 128
        sig = np.arange(128)

        win = signal.get_window(('kaiser', 8.0), osfactor // 2)
        assert_raises(ValueError, signal.resample, (sig, len(sig) * osfactor), {'window': win})


def test_windowfunc_basics():
    for window_name, params in window_funcs:
        window = getattr(signal, window_name)
        w1 = window(7, *params, sym=True)
        w2 = window(7, *params, sym=False)
        assert_array_almost_equal(w1, w2)
        # just check the below runs
        window(6, *params, sym=True)
        window(6, *params, sym=False)


if __name__ == "__main__":
    run_module_suite()
