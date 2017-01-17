from __future__ import division, print_function, absolute_import
from scipy._lib.six import with_metaclass
import numpy as np
from numpy.testing import dec, run_module_suite

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
    mp = MissingModule('mpmath')


_is_32bit_platform = np.intp(0).itemsize < 8


class TestInversion(with_metaclass(DecoratorMeta, object)):
    decorators = [(dec.slow, None),
                  (check_version, (sympy, '0.7')),
                  (check_version, (mp, '0.19'))]

    @dec.knownfailureif(_is_32bit_platform, "rtol only 2e-9, see gh-6938")
    def test_log(self):
        with mp.workdps(30):
            logcoeffs = mp.taylor(lambda x: mp.log(1 + x), 0, 10)
            expcoeffs = mp.taylor(lambda x: mp.exp(x) - 1, 0, 10)
            invlogcoeffs = lagrange_inversion(logcoeffs)
            mp_assert_allclose(invlogcoeffs, expcoeffs)

    @dec.knownfailureif(_is_32bit_platform, "rtol only 1e-15, see gh-6938")
    def test_sin(self):
        with mp.workdps(30):
            sincoeffs = mp.taylor(mp.sin, 0, 10)
            asincoeffs = mp.taylor(mp.asin, 0, 10)
            invsincoeffs = lagrange_inversion(sincoeffs)
            mp_assert_allclose(invsincoeffs, asincoeffs, atol=1e-30)


if __name__ == "__main__":
    run_module_suite()
