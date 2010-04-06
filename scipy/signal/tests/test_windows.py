
from numpy import array, ones_like
from numpy.testing import assert_array_almost_equal, assert_array_equal
from scipy import signal


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

    def test_cheb_odd(self):
        cheb_odd = signal.chebwin(53, at=-40)
        assert_array_almost_equal(cheb_odd, cheb_odd_true, decimal=4)

    def test_cheb_even(self):
        cheb_even = signal.chebwin(54, at=-40)
        assert_array_almost_equal(cheb_even, cheb_even_true, decimal=4)


class TestGetWindow(object):

    def test_boxcar(self):
        w = signal.get_window('boxcar', 12)
        assert_array_equal(w, ones_like(w))

    def test_cheb_odd(self):
        w = signal.get_window(('chebwin', -40), 53, fftbins=False)
        assert_array_almost_equal(w, cheb_odd_true, decimal=4)

    def test_cheb_even(self):
        w = signal.get_window(('chebwin', -40), 54, fftbins=False)
        assert_array_almost_equal(w, cheb_even_true, decimal=4)
