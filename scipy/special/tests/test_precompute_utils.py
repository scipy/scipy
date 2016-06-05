from __future__ import division, print_function, absolute_import
from scipy._lib.six import with_metaclass
from numpy.testing import dec

from scipy.special._testutils import MissingModule, check_version, DecoratorMeta
from scipy.special._mptestutils import mp_assert_allclose
from scipy.special._precompute.utils import lagrange_inversion

try:
    import sympy
except ImportError:
    sympy = MissingModule('sympy')

try:
    import mpmath as mp
except ImportError:
    try:
        from sympy import mpmath as mp
    except ImportError:
        mp = MissingModule('mpmath')


class TestInversion(with_metaclass(DecoratorMeta, object)):
    decorators = [(dec.slow, None),
                  (check_version, (sympy, '0.7')),
                  (check_version, (mp, '0.19'))]

    def test_log(self):
        with mp.workdps(30):
            logcoeffs = mp.taylor(lambda x: mp.log(1 + x), 0, 10)
            expcoeffs = mp.taylor(lambda x: mp.exp(x) - 1, 0, 10)
            invlogcoeffs = lagrange_inversion(logcoeffs)
            mp_assert_allclose(invlogcoeffs, expcoeffs)

    def test_sin(self):
        with mp.workdps(30):
            sincoeffs = mp.taylor(mp.sin, 0, 10)
            asincoeffs = mp.taylor(mp.asin, 0, 10)
            invsincoeffs = lagrange_inversion(sincoeffs)
            mp_assert_allclose(invsincoeffs, asincoeffs, atol=1e-30)
