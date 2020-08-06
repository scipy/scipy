import numpy as np
from numpy.testing import assert_equal, assert_allclose
from scipy.special import iv, logiv, logive
from scipy.special._log_iv import _log_iv_asym


def assert_iv(v, z, expected_iv):
    # Underlying asymptotic implementation
    assert_allclose(_log_iv_asym(v, z), expected_iv, rtol=1e-10)

    # Wrapper which chooses between log(iv()) and log_iv_asym()
    assert_allclose(logiv(v, z), expected_iv, rtol=3e-14)
    assert_allclose(logive(v, z), expected_iv - abs(z), 3e-14)

class TestLogIv(object):
    def test_int(self):
        assert_iv(55, 27, -22.01164919937816)

    def test_float(self):
        assert_iv(55.0, 27, -22.01164919937816)
        assert_iv(55, 27.0, -22.01164919937816)
        assert_iv(55.0, 27.0, -22.01164919937816)

    def test_nan(self):
        assert_iv(5, np.inf, np.nan)
        assert_iv(5, np.nan, np.nan)
        assert_iv(np.nan, 5, np.nan)

    def test_inf(self):
        assert_iv(np.inf, 5, -np.inf)

    def test_large(self):
        assert_iv(50, 1, np.log(iv(50, 1)))
        assert_iv(100, 1, -433.0516183940659)
        assert_iv(400, 1, -2277.758946766306)
        assert_iv(800, 1, -5106.468163036196)

        assert_iv(100, 10, np.log(iv(100, 10)))
        assert_iv(100, 40, np.log(iv(100, 40)))
        assert_iv(100, 80, np.log(iv(100, 80)))

        sqrt2 = np.sqrt(2)
        assert_iv(100, 10*sqrt2, np.log(iv(100, 10*sqrt2)))
        assert_iv(100*sqrt2, 10, np.log(iv(100*sqrt2, 10)))

        assert_iv(1000, 10, -4302.665291340331)
        assert_iv(1000, 10*sqrt2, -3956.066726969071)
        assert_iv(1000, 100, -1997.610772811001)
        assert_iv(1000, 100*sqrt2, -1648.548946015068)
        assert_iv(1000, 1000, 528.2938870936566)
        assert_iv(1000, 1000*sqrt2, 1068.924421927351)

        assert_iv(60, 100, 79.19174702982741)
        assert_iv(60*sqrt2, 1000*sqrt2, 1407.121827193677)
