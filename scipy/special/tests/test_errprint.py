from __future__ import division, print_function, absolute_import

import warnings

from numpy.testing import assert_

from scipy import special
from scipy.special import cython_special


def _check_errprint(mod):
    flag = mod.errprint(True)
    try:
        assert_(isinstance(flag, bool))
        with warnings.catch_warnings(record=True) as w:
            mod.loggamma(0)
            assert_(w[-1].category is special.SpecialFunctionWarning)
    finally:
        mod.errprint(flag)


def test_special_errprint():
    _check_errprint(special)


def test_cython_special_errprint():
    _check_errprint(cython_special)
