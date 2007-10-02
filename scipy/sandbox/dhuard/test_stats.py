"""
Test statistical functions.
"""


from numpy.testing import *
import stats
import numpy as np

N = 100
np.random.seed(2)
r = np.random.randn(N)

class TestEmpiricalCDF(NumpyTestCase):
        def check_hazen(self):

                f = stats.empiricalcdf(r)
                assert_equal(len(f), len(r))
                assert_array_equal(np.argsort(r), np.argsort(f))
                assert_array_equal(np.sort(f), (np.arange(N)+.5)/N)

        def check_weibull(self):
                f = stats.empiricalcdf(r, 'weibull')
                assert_array_equal(np.sort(f), (np.arange(N)+1.)/(N+1.))

        def check_california(self):
                f = stats.empiricalcdf(r, 'california')
                assert_array_equal(np.sort(f), (np.arange(N))/float(N))

class TestScoreAtPercentile(NumpyTestCase):
        def check_simple(self):
                r = np.random.randn(1000)
                s = stats.scoreatpercentile(r, [15.9,50,84.1])
                assert_array_almost_equal(s, [-1,0,1], 1)

class TestPercentileOfScore(NumpyTestCase):
        def check_simple(self):
                r = np.random.randn(3000)
                p = stats.percentileofscore(r, [-1,0,1])
                assert_array_almost_equal(p, [15.9, 50, 84.1], 0)



if __name__ == '__main__':
        NumpyTest().run()
