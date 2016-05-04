from __future__ import division, print_function, absolute_import
from scipy._lib.six import with_metaclass
from numpy.testing import dec

from scipy.special._testutils import MissingModule, check_version, SystematicMeta
from scipy.special._precompute import utils

try:
    import sympy
    from sympy import mpmath as mp
except ImportError:
    sympy = MissingModule('sympy')


class TestInversion(with_metaclass(SystematicMeta, object)):
    decodict = {dec.slow: (), check_version: (sympy, '1.0')}

    def test_log(self):
        with mp.workdps(30):
            logcoeffs = mp.taylor(lambda x: mp.log(1 + x), 0, 10)
            expcoeffs = mp.taylor(lambda x: mp.exp(x) - 1, 0, 10)
            invlogcoeffs = utils.lagrange_inversion(logcoeffs)
            utils.mpf_assert_allclose(invlogcoeffs, expcoeffs)

    def test_sin(self):
        with mp.workdps(30):
            sincoeffs = mp.taylor(mp.sin, 0, 10)
            asincoeffs = mp.taylor(mp.asin, 0, 10)
            invsincoeffs = utils.lagrange_inversion(sincoeffs)
            utils.mpf_assert_allclose(invsincoeffs, asincoeffs, atol=1e-30)
