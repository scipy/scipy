#!/usr/bin/env python

""" Test functions for the 'dictsampler' derived class.

Author: Ed Schofield, 2003-2006
Copyright: Ed Schofield, 2003-2006

"""

from numpy import arange, add, array, dot, zeros, identity

import sys
from scipy.testing import *

from numpy import sum, array, average, sqrt
#from scipy.montecarlo import *
from scipy.sandbox.montecarlo import *
from scipy import stats


class test_dict_sampler(TestCase):
    def test_simple(self):
        """
        # Sample from this discrete distribution:
        #    x       'a'       'b'       'c'
        #    p(x)    10/180    150/180   20/180
        """
        table = {}
        table['a'] = 10
        table['b'] = 150
        table['c'] = 20
        sampler = dictsampler(table)
        n = 1000
        #import pdb
        #pdb.set_trace()
        s = sampler.sample(n)
        assert sum((s[i]=='b' for i in range(n)),axis=0)*1./n > 0.75

        #lam = 10.0
        #n = 35
        #numsamples = 10000
        #a = array([stats.poisson.pmf(i, lam) for i in range(n)])
        #sampler = intsampler(a)
        #x = sampler.sample(numsamples)
        #m = x.mean()
        #s = x.std()
        ## Use a normal approximation for confidence intervals for the mean
        #z = 2.5758   # = norminv(0.995), for a 1% confidence interval
        #assert abs(m - lam) < z * lam/sqrt(numsamples)

    def test_sanity(self):
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
    nose.run(argv=['', __file__])
