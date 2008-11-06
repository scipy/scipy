import warnings

import numpy as np
from numpy.testing import TestCase, assert_array_almost_equal

from scipy.signal import tf2zpk, bessel, BadCoefficients

class TestTf2zpk(TestCase):
    def test_simple(self):
        z_r = np.array([0.5, -0.5])
        p_r = np.array([1.j / np.sqrt(2), -1.j / np.sqrt(2)])
        # Sort the zeros/poles so that we don't fail the test if the order
        # changes
        z_r.sort()
        p_r.sort()
        b = np.poly(z_r)
        a = np.poly(p_r)

        z, p, k = tf2zpk(b, a)
        z.sort()
        p.sort()
        assert_array_almost_equal(z, z_r)
        assert_array_almost_equal(p, p_r)

    def test_bad_filter(self):
        """Regression test for #651: better handling of badly conditionned
        filter coefficients."""
        b, a = bessel(20, 0.1)
        warnings.simplefilter("error", BadCoefficients)
        try:
            try:
                z, p, k = tf2zpk(b, a)
                raise AssertionError("tf2zpk did not warn about bad "\
                                     "coefficients")
            except BadCoefficients:
                pass
        finally:
            warnings.simplefilter("always", BadCoefficients)
