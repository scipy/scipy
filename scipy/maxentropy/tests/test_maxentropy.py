#!/usr/bin/env python

""" Test functions for maximum entropy module.

Author: Ed Schofield, 2003-2005
Copyright: Ed Schofield, 2003-2005
"""

from numpy.testing import *
from numpy import arange, log, exp, ones
from scipy.maxentropy.maxentropy import *

class TestMaxentropy(TestCase):
    """Test whether logsumexp() function correctly handles large
    inputs.
    """
    def test_logsumexp(self):
        a = arange(200)
        desired = log(sum(exp(a)))
        assert_almost_equal(logsumexp(a), desired)

        # Now test with large numbers
        b = [1000,1000]
        desired = 1000.0 + log(2.0)
        assert_almost_equal(logsumexp(b), desired)

        n = 1000
        b = ones(n)*10000
        desired = 10000.0 + log(n)
        assert_almost_equal(logsumexp(b), desired)

    def test_simple(self):
        # Write me!
        pass


if __name__ == "__main__":
    run_module_suite()
