#!/usr/bin/env python

""" Test functions for maximum entropy module.

Author: Ed Schofield, 2003-2005
Copyright: Ed Schofield, 2003-2005

"""

__usage__ = """
Build linalg:
  python setup.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.maxent.test(<level>)'
Run tests if maxent is not installed:
  python tests/test_maxent.py [<level>]
"""

import scipy.base as Numeric
from scipy.base import arange, add, array, dot, zeros, identity

import sys
from scipy.test.testing import *
set_package_path()
from scipy.base import *
from maxent import *
restore_path()

import unittest


class test_logsumexp(ScipyTestCase):
    """Test whether logsumexp() function correctly and handles large
    inputs.
    """
    def check_simple(self):
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


class test_maxent(ScipyTestCase):
    def check_simple(self):
        # Write me!
        pass


if __name__ == "__main__":
    ScipyTest().run()
