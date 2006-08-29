#!/usr/bin/env python

""" Test functions for generic discrete sampler 'intsampler'

Author: Ed Schofield, 2003-2006
Copyright: Ed Schofield, 2003-2006

"""

__usage__ = """
Build intsampler:
  python setup.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.intsampler.test(<level>)'
Run tests if intsampler is not installed:
  python tests/test_intsampler.py [<level>]
"""

from numpy import arange, add, array, dot, zeros, identity

import sys
from numpy.testing import *
set_package_path()
from numpy import *
#from scipy.montecarlo import *
from scipy.sandbox.montecarlo import *
from scipy import stats
restore_path()

import unittest


class test_intsampler(ScipyTestCase):
    def check_simple(self):
        # Sample from a Poisson distribution, P(lambda = 10.0)
        lam = 10.0
        n = 35
        numsamples = 10000
        a = array([stats.poisson.pmf(i, lam) for i in range(n)])
        sampler = intsampler(a)
        x = sampler.sample(numsamples)
        m = x.mean()
        s = x.std()
        # Use a normal approximation for confidence intervals for the mean
        z = 2.5758   # = norminv(0.995), for a 1% confidence interval
        assert abs(m - lam) < z * lam/sqrt(numsamples)

    def check_sanity(self):
        # Sample from this pmf:
        #      x        0       1       2       3       4
        #      p(x)     0.5     0.1     0.15    0       0.25

        # The true mean and variance are:
        truemean = 1.4
        truevar = 2.74

        numsamples = 10000
        a = array([0.5, 0.1, 0.15, 0., 0.25])
        sampler = intsampler(a)
        x = sampler.sample(numsamples)
        m = x.mean()
        v = x.var()
        assert x.max() == 4
        assert x.min() == 0
        assert sum(x==3,axis=0) == 0
        assert 0.08 < average(x==1,axis=0) < 0.12
        # Use a normal approx for confidence intervals for the mean
        z = 2.5758   # = norminv(0.995), for a 1% confidence interval
        assert abs(m - truemean) < z * sqrt(truevar/numsamples)



if __name__ == "__main__":
    ScipyTest().run()
