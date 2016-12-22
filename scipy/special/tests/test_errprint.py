from __future__ import division, print_function, absolute_import

import warnings

from numpy.testing import assert_

import scipy. special as sc


def test_errprint():
    flag = sc.errprint(True)
    try:
        assert_(isinstance(flag, bool))
        with warnings.catch_warnings(record=True) as w:
            sc.loggamma(0)
            assert_(w[-1].category is sc.SpecialFunctionWarning)
    finally:
        sc.errprint(flag)
