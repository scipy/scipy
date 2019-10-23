import numpy as np
from numpy import isnan

import pytest
from numpy.testing import (assert_almost_equal, assert_array_equal)

from scipy import special
import scipy.special._ufuncs as cephes


class TestCephes(object):
    # old tests with int args
    def test_pdtr(self):
        val = cephes.pdtr(0, 1)
        assert_almost_equal(val, np.exp(-1))
        # Edge case: m = 0.
        val = cephes.pdtr([0, 1, 2], 0)
        assert_array_equal(val, [1, 1, 1])
        # test the rounding from double -> int
        double_val = cephes.pdtr([0.1, 1.1, 2.1], 1.0)
        int_val = cephes.pdtr([0, 1, 2], 1.0)
        assert_array_equal(double_val, int_val)

    def test_pdtrc(self):
        val = cephes.pdtrc(0, 1)
        assert_almost_equal(val, 1 - np.exp(-1))
        # Edge case: m = 0.
        val = cephes.pdtrc([0, 1, 2], 0.0)
        assert_array_equal(val, [0, 0, 0])
        # test the rounding from double -> int
        double_val = cephes.pdtrc([0.1, 1.1, 2.1], 1.0)
        int_val = cephes.pdtrc([0, 1, 2], 1.0)
        assert_array_equal(double_val, int_val)

    #test np.inf arguments
    def test_pdtr_np_inf(self):
        val = cephes.pdtr(np.inf, 1.0)
        assert_almost_equal(val, 1.0)

    def test_pdtrc_np_inf(self):
        val = cephes.pdtrc(np.inf, 1.0)
        assert_almost_equal(val, 0.0)

    #test -1.1 is np.nan, negative values are not valid domain
    def test_pdtr_np_inf(self):
        val = cephes.pdtr(-1.1, 1.0)
        assert isnan(val)

    def test_pdtrc_np_inf(self):
        val = cephes.pdtrc(-1.1, 1.0)
        assert isnan(val)
