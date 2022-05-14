
import pytest
import numpy as np
from numpy.testing import assert_allclose
from scipy.stats import _boost


rtols = {np.float32: 32*np.finfo(np.float32).eps,
         np.float64: 32*np.finfo(np.float64).eps,
         np.longdouble: 32*np.finfo(np.longdouble).eps}


# Each item in this list is
#   (func, args, expected_value)
# All the values can be represented exactly, even with np.float32.
#
# This is not an exhaustive test data set of all the functions!
# It is a spot check of several functions, primarily for
# checking that the different data types are handled correctly.
test_data = [
    (_boost._beta_cdf, (0.5, 2, 3), 0.6875),
    (_boost._beta_ppf, (0.6875, 2, 3), 0.5),
    (_boost._beta_pdf, (0.5, 2, 3), 1.5),
    (_boost._beta_sf, (0.5, 2, 1), 0.75),
    (_boost._beta_isf, (0.75, 2, 1), 0.5),
    (_boost._binom_cdf, (1, 3, 0.5), 0.5),
    (_boost._binom_pdf, (1, 4, 0.5), 0.25),
    (_boost._hypergeom_cdf, (2, 3, 5, 6), 0.5),
    (_boost._nbinom_cdf, (1, 4, 0.25), 0.015625),
    (_boost._ncf_mean, (10, 12, 2.5), 1.5),
]


@pytest.mark.parametrize('func, args, expected', test_data)
@pytest.mark.parametrize('typ, rtol', rtols.items())
def test_stats_boost_ufunc(func, args, expected, typ, rtol):
    args = (typ(arg) for arg in args)
    value = func(*args)
    assert_allclose(value, expected, rtol=rtol)
