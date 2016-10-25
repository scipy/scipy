from __future__ import division, print_function, absolute_import

import warnings

from numpy.testing import dec, assert_

from scipy import special
from scipy.special import _test_sf_error
from scipy.special.cython_special import have_openmp


def test_errprint():
    flag = special.errprint(True)
    try:
        assert_(isinstance(flag, bool))
        with warnings.catch_warnings(record=True) as w:
            special.loggamma(0)
            assert_(w[-1].category is special.SpecialFunctionWarning)
    finally:
        special.errprint(flag)


def test_cython_special_error_cfunc():
    _test_sf_error.test_cython_special_error_cfunc()


def test_cython_special_error_cyfunc():
    _test_sf_error.test_cython_special_error_cyfunc()


def test_cython_special_error_cppfunc():
    _test_sf_error.test_cython_special_error_cppfunc()


@dec.skipif(not have_openmp())
def test_cython_special_error_parallel():
    _test_sf_error.test_cython_special_error_parallel()


def test_cython_special_error_serial():
    _test_sf_error.test_cython_special_error_serial()
