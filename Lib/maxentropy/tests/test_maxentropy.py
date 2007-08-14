#!/usr/bin/env python

""" Test functions for maximum entropy module.

Author: Ed Schofield, 2003-2005
Copyright: Ed Schofield, 2003-2005
"""

import sys
from numpy.testing import *
from numpy import arange, add, array, dot, zeros, identity, log, exp, ones
set_package_path()
from scipy.maxentropy.maxentropy import *
restore_path()

import unittest


class test_maxentropy(NumpyTestCase):
    """Test whether logsumexp() function correctly handles large
    inputs.
    """
    def check_logsumexp(self, level=1):
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

    def check_simple(self, level=1):
        # Write me!
        pass


if __name__ == "__main__":
    NumpyTest().run()
