from __future__ import division, print_function, absolute_import
from scipy.special._precompute import utils

try:
    import mpmath as mp
except ImportError:
    try:
        import sympy.mpmath as mp
    except ImportError:
        pass


class TestInversion():
    @utils.skip()
    def test_log(self):
        with mp.workdps(30):
            logcoeffs = mp.taylor(lambda x: mp.log(1 + x), 0, 10)
            expcoeffs = mp.taylor(lambda x: mp.exp(x) - 1, 0, 10)
            invlogcoeffs = utils.lagrange_inversion(logcoeffs)
            utils.mpf_assert_allclose(invlogcoeffs, expcoeffs)

    @utils.skip()
    def test_sin(self):
        with mp.workdps(30):
            sincoeffs = mp.taylor(mp.sin, 0, 10)
            asincoeffs = mp.taylor(mp.asin, 0, 10)
            invsincoeffs = utils.lagrange_inversion(sincoeffs)
            utils.mpf_assert_allclose(invsincoeffs, asincoeffs, atol=1e-30)
