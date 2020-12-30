from __future__ import division, print_function, absolute_import

import numpy as np
from itertools import product
from numpy.testing import (assert_, assert_equal, assert_allclose,
                           assert_almost_equal)  # avoid new uses

import pytest
from pytest import raises as assert_raises
from scipy.stats._hypotests import (epps_singleton_2samp, cramervonmises,
                                    _cdf_cvm)
from scipy.stats._mannwhitneyu import mannwhitneyu2, _mwu_state
from scipy.stats import distributions
from .common_tests import check_named_results


class TestEppsSingleton(object):
    def test_statistic_1(self):
        # first example in Goerg & Kaiser, also in original paper of
        # Epps & Singleton. Note: values do not match exactly, the
        # value of the interquartile range varies depending on how
        # quantiles are computed
        x = np.array([-0.35, 2.55, 1.73, 0.73, 0.35,
                      2.69, 0.46, -0.94, -0.37, 12.07])
        y = np.array([-1.15, -0.15, 2.48, 3.25, 3.71,
                      4.29, 5.00, 7.74, 8.38, 8.60])
        w, p = epps_singleton_2samp(x, y)
        assert_almost_equal(w, 15.14, decimal=1)
        assert_almost_equal(p, 0.00442, decimal=3)

    def test_statistic_2(self):
        # second example in Goerg & Kaiser, again not a perfect match
        x = np.array((0, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 5, 5, 5, 5, 6, 10,
                      10, 10, 10))
        y = np.array((10, 4, 0, 5, 10, 10, 0, 5, 6, 7, 10, 3, 1, 7, 0, 8, 1,
                      5, 8, 10))
        w, p = epps_singleton_2samp(x, y)
        assert_allclose(w, 8.900, atol=0.001)
        assert_almost_equal(p, 0.06364, decimal=3)

    def test_epps_singleton_array_like(self):
        np.random.seed(1234)
        x, y = np.arange(30), np.arange(28)

        w1, p1 = epps_singleton_2samp(list(x), list(y))
        w2, p2 = epps_singleton_2samp(tuple(x), tuple(y))
        w3, p3 = epps_singleton_2samp(x, y)

        assert_(w1 == w2 == w3)
        assert_(p1 == p2 == p3)

    def test_epps_singleton_size(self):
        # raise error if less than 5 elements
        x, y = (1, 2, 3, 4), np.arange(10)
        assert_raises(ValueError, epps_singleton_2samp, x, y)

    def test_epps_singleton_nonfinite(self):
        # raise error if there are non-finite values
        x, y = (1, 2, 3, 4, 5, np.inf), np.arange(10)
        assert_raises(ValueError, epps_singleton_2samp, x, y)
        x, y = np.arange(10), (1, 2, 3, 4, 5, np.nan)
        assert_raises(ValueError, epps_singleton_2samp, x, y)

    def test_epps_singleton_1d_input(self):
        x = np.arange(100).reshape(-1, 1)
        assert_raises(ValueError, epps_singleton_2samp, x, x)

    def test_names(self):
        x, y = np.arange(20), np.arange(30)
        res = epps_singleton_2samp(x, y)
        attributes = ('statistic', 'pvalue')
        check_named_results(res, attributes)


class TestCvm(object):
    # the expected values of the cdfs are taken from Table 1 in
    # Csorgo / Faraway: The Exact and Asymptotic Distribution of
    # Cramér-von Mises Statistics, 1996.
    def test_cdf_4(self):
        assert_allclose(
                _cdf_cvm([0.02983, 0.04111, 0.12331, 0.94251], 4),
                [0.01, 0.05, 0.5, 0.999],
                atol=1e-4)

    def test_cdf_10(self):
        assert_allclose(
                _cdf_cvm([0.02657, 0.03830, 0.12068, 0.56643], 10),
                [0.01, 0.05, 0.5, 0.975],
                atol=1e-4)

    def test_cdf_1000(self):
        assert_allclose(
                _cdf_cvm([0.02481, 0.03658, 0.11889, 1.16120], 1000),
                [0.01, 0.05, 0.5, 0.999],
                atol=1e-4)

    def test_cdf_inf(self):
        assert_allclose(
                _cdf_cvm([0.02480, 0.03656, 0.11888, 1.16204]),
                [0.01, 0.05, 0.5, 0.999],
                atol=1e-4)

    def test_cdf_support(self):
        # cdf has support on [1/(12*n), n/3]
        assert_equal(_cdf_cvm([1/(12*533), 533/3], 533), [0, 1])
        assert_equal(_cdf_cvm([1/(12*(27 + 1)), (27 + 1)/3], 27), [0, 1])

    def test_cdf_large_n(self):
        # test that asymptotic cdf and cdf for large samples are close
        assert_allclose(
                _cdf_cvm([0.02480, 0.03656, 0.11888, 1.16204, 100], 10000),
                _cdf_cvm([0.02480, 0.03656, 0.11888, 1.16204, 100]),
                atol=1e-4)

    def test_large_x(self):
        # for large values of x and n, the series used to compute the cdf
        # converges slowly.
        # this leads to bug in R package goftest and MAPLE code that is
        # the basis of the implemenation in scipy
        # note: cdf = 1 for x >= 1000/3 and n = 1000
        assert_(0.99999 < _cdf_cvm(333.3, 1000) < 1.0)
        assert_(0.99999 < _cdf_cvm(333.3) < 1.0)

    def test_low_p(self):
        # _cdf_cvm can return values larger than 1. In that case, we just
        # return a p-value of zero.
        n = 12
        res = cramervonmises(np.ones(n)*0.8, 'norm')
        assert_(_cdf_cvm(res.statistic, n) > 1.0)
        assert_equal(res.pvalue, 0)

    def test_invalid_input(self):
        x = np.arange(10).reshape((2, 5))
        assert_raises(ValueError, cramervonmises, x, "norm")
        assert_raises(ValueError, cramervonmises, [1.5], "norm")
        assert_raises(ValueError, cramervonmises, (), "norm")

    def test_values_R(self):
        # compared against R package goftest, version 1.1.1
        # goftest::cvm.test(c(-1.7, 2, 0, 1.3, 4, 0.1, 0.6), "pnorm")
        res = cramervonmises([-1.7, 2, 0, 1.3, 4, 0.1, 0.6], "norm")
        assert_allclose(res.statistic, 0.288156, atol=1e-6)
        assert_allclose(res.pvalue, 0.1453465, atol=1e-6)

        # goftest::cvm.test(c(-1.7, 2, 0, 1.3, 4, 0.1, 0.6),
        #                   "pnorm", mean = 3, sd = 1.5)
        res = cramervonmises([-1.7, 2, 0, 1.3, 4, 0.1, 0.6], "norm", (3, 1.5))
        assert_allclose(res.statistic, 0.9426685, atol=1e-6)
        assert_allclose(res.pvalue, 0.002026417, atol=1e-6)

        # goftest::cvm.test(c(1, 2, 5, 1.4, 0.14, 11, 13, 0.9, 7.5), "pexp")
        res = cramervonmises([1, 2, 5, 1.4, 0.14, 11, 13, 0.9, 7.5], "expon")
        assert_allclose(res.statistic, 0.8421854, atol=1e-6)
        assert_allclose(res.pvalue, 0.004433406, atol=1e-6)

    def test_callable_cdf(self):
        x, args = np.arange(5), (1.4, 0.7)
        r1 = cramervonmises(x, distributions.expon.cdf)
        r2 = cramervonmises(x, "expon")
        assert_equal((r1.statistic, r1.pvalue), (r2.statistic, r2.pvalue))

        r1 = cramervonmises(x, distributions.beta.cdf, args)
        r2 = cramervonmises(x, "beta", args)
        assert_equal((r1.statistic, r1.pvalue), (r2.statistic, r2.pvalue))


class TestMannWhitneyU():

    # These are tabulated values of the CDF of the exact distribution of
    # the test statistic from pg 52 of reference [1] (Mann-Whitney Original)
    pn3 = {1: [0.25, 0.5, 0.75], 2: [0.1, 0.2, 0.4, 0.6],
           3: [0.05, .1, 0.2, 0.35, 0.5, 0.65]}
    pn4 = {1: [0.2, 0.4, 0.6], 2: [0.067, 0.133, 0.267, 0.4, 0.6],
           3: [0.028, 0.057, 0.114, 0.2, .314, 0.429, 0.571],
           4: [0.014, 0.029, 0.057, 0.1, 0.171, 0.243, 0.343, 0.443, 0.557]}
    pm5 = {1: [0.167, 0.333, 0.5, 0.667],
           2: [0.047, 0.095, 0.19, 0.286, 0.429, 0.571],
           3: [0.018, 0.036, 0.071, 0.125, 0.196, 0.286, 0.393, 0.5, 0.607],
           4: [0.008, 0.016, 0.032, 0.056, 0.095, 0.143,
               0.206, 0.278, 0.365, 0.452, 0.548],
           5: [0.004, 0.008, 0.016, 0.028, 0.048, 0.075, 0.111,
               0.155, 0.21, 0.274, 0.345, .421, 0.5, 0.579]}
    pm6 = {1: [0.143, 0.286, 0.428, 0.571],
           2: [0.036, 0.071, 0.143, 0.214, 0.321, 0.429, 0.571],
           3: [0.012, 0.024, 0.048, 0.083, 0.131,
               0.19, 0.274, 0.357, 0.452, 0.548],
           4: [0.005, 0.01, 0.019, 0.033, 0.057, 0.086, 0.129,
               0.176, 0.238, 0.305, 0.381, 0.457, 0.543],  # the last element
              # of the previous list, 0.543, has been modified from 0.545;
              # I assume it was a typo
           5: [0.002, 0.004, 0.009, 0.015, 0.026, 0.041, 0.063, 0.089,
               0.123, 0.165, 0.214, 0.268, 0.331, 0.396, 0.465, 0.535],
           6: [0.001, 0.002, 0.004, 0.008, 0.013, 0.021, 0.032, 0.047,
               0.066, 0.09, 0.12, 0.155, 0.197, 0.242, 0.294, 0.350,
               0.409, 0.469, 0.531]}

    def test_exact_distribution(self):
        # I considered parametrize. I decided against it.
        p_tables = {3: self.pn3, 4: self.pn4, 5: self.pm5, 6: self.pm6}
        for n, table in p_tables.items():
            for m, p in table.items():
                # check p-value against table
                u = np.arange(0, len(p))
                assert_allclose(_mwu_state.cdf(k=u, m=m, n=n), p, atol=1e-3)

                # check identity CDF + SF - PMF = 1
                # ( In this implementation, SF(U) includes PMF(U) )
                u2 = np.arange(0, m*n+1)
                assert_allclose(_mwu_state.cdf(k=u2, m=m, n=n)
                                + _mwu_state.sf(k=u2, m=m, n=n)
                                - _mwu_state.pmf(k=u2, m=m, n=n), 1)

                # check symmetry about mean of U, i.e. pmf(U) = pmf(m*n-U)
                pmf = _mwu_state.pmf(k=u2, m=m, n=n)
                assert_allclose(pmf, pmf[::-1])

                # check symmetry w.r.t. interchange of m, n
                pmf2 = _mwu_state.pmf(k=u2, m=n, n=m)
                assert_allclose(pmf, pmf2)

    def test_input_validation(self):
        x = np.array([1, 2]) # generic, valid inputs
        y = np.array([3, 4])
        with assert_raises(ValueError, match="`x` and `y` must be of nonzero"):
            mannwhitneyu2([], y)
        with assert_raises(ValueError, match="`x` and `y` must be of nonzero"):
            mannwhitneyu2(x, [])
        with assert_raises(ValueError, match="`x` and `y` must not contain"):
            mannwhitneyu2([np.nan, 2], y)
        with assert_raises(ValueError, match="`use_continuity` must be one"):
            mannwhitneyu2(x, y, use_continuity='ekki')
        with assert_raises(ValueError, match="`alternative` must be one of"):
            mannwhitneyu2(x, y, alternative='ekki')
        with assert_raises(ValueError, match="`axis` must be an integer"):
            mannwhitneyu2(x, y, axis=1.5)
        with assert_raises(ValueError, match="`method` must be one of"):
            mannwhitneyu2(x, y, method='ekki')

    def test_corner_cases(self):
        # tests behavior for cases that required special handling
        pass
        # mannwhitneyu2([1], [1])
        # mannwhitneyu2([1, 2], [1, 2])

    def test_unusual_cases(self):
        # test behavior for unusual cases that did not require special handling
        mannwhitneyu2(1, 2, alternative='two-sided')
        mannwhitneyu2(1, 2, alternative='greater')
        mannwhitneyu2(1, 2, alternative='less')

    @pytest.mark.parametrize("method", ["asymptotic", "exact"])
    def test_gh_12837_11113(self, method):
        # test that behavior for broadcastable nd arrays is appropriate:
        # output shape is correct and all values are equal to when the test
        # is performed on one pair of samples at a time.
        # tests that gh-12837 and gh-11113 (requests for n-d input)
        # are resolved
        np.random.seed(0)

        # arrays are broadcastable except for axis = -3
        axis = -3
        m, n = 7, 10  # sample sizes
        x = np.random.rand(m, 3, 8)
        y = np.random.rand(6, n, 1, 8) + 0.1
        res = mannwhitneyu2(x, y, method=method, axis=axis)

        shape = (6, 3, 8) # appropriate shape of outputs, given inputs
        assert(res.pvalue.shape == shape)
        assert(res.statistic.shape == shape)

        # move axis of test to end for simplicity
        x, y = np.moveaxis(x, axis, -1), np.moveaxis(y, axis, -1)

        x = x[None, ...] # give x a zeroth dimension
        assert(x.ndim == y.ndim)

        x = np.broadcast_to(x, shape + (m,))
        y = np.broadcast_to(y, shape + (n,))
        assert(x.shape[:-1] == shape)
        assert(y.shape[:-1] == shape)

        # loop over pairs of samples
        statistics = np.zeros(shape)
        pvalues = np.zeros(shape)
        for indices in product(*[range(i) for i in shape]):
            xi = x[indices]
            yi = y[indices]
            temp = mannwhitneyu2(xi, yi, method=method)
            statistics[indices] = temp.statistic
            pvalues[indices] = temp.pvalue

        np.testing.assert_equal(res.pvalue, pvalues)
        np.testing.assert_equal(res.statistic, statistics)

    def test_gh_11355(self):
        x = [1, 2, 3, 4]
        y = [3, 6, 7, 8, 9, 3, 2, 1, 4, 4, 5]

        res1 = mannwhitneyu2(x, y)

        y[4] = np.nan
        with assert_raises(ValueError, match="`x` and `y` must not contain"):
            mannwhitneyu2(x, y)

        y[4] = np.inf
        res2 = mannwhitneyu2([1, 2, 3, 4], [3, 6, 7, 8, np.inf, 3, 2, 1, 4, 4, 5])

        assert_equal(res1.statistic, res2.statistic)
        assert_equal(res1.pvalue, res1.pvalue)

    settings = [[True, "less", "asymptotic", 0.900775348204],
                [True, "greater", "asymptotic", 0.1223118025635],
                [True, "two-sided", "asymptotic", 0.244623605127],
                [False, "less", "asymptotic", 0.8896643190401],
                [False, "greater", "asymptotic", 0.1103356809599],
                [False, "two-sided", "asymptotic", 0.2206713619198],
                [True, "less", "exact", 0.8967698967699],
                [True, "greater", "exact", 0.1272061272061],
                [True, "two-sided", "exact", 0.2544122544123]]
    @pytest.mark.parametrize(("use_continuity", "alternative",
                              "method", "pvalue_exp"), settings)
    def test_gh_9184(self, use_continuity, alternative, method, pvalue_exp):
        # gh-9184 might be considered a doc-only bug. Please see the
        # documentation to confirm that mannwhitneyu2 correctly notes
        # that the output statistic is that of the first sample (x). In any
        # case, check the case provided there against output from R.
        # R code:
        # options(digits=16)
        # x <- c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
        # y <- c(1.15, 0.88, 0.90, 0.74, 1.21)
        # wilcox.test(x, y, alternative = "less", exact = FALSE)
        # wilcox.test(x, y, alternative = "greater", exact = FALSE)
        # wilcox.test(x, y, alternative = "two.sided", exact = FALSE)
        # wilcox.test(x, y, alternative = "less", exact = FALSE,
        #             correct=FALSE)
        # wilcox.test(x, y, alternative = "greater", exact = FALSE,
        #             correct=FALSE)
        # wilcox.test(x, y, alternative = "two.sided", exact = FALSE,
        #             correct=FALSE)
        # wilcox.test(x, y, alternative = "less", exact = TRUE)
        # wilcox.test(x, y, alternative = "greater", exact = TRUE)
        # wilcox.test(x, y, alternative = "two.sided", exact = TRUE)
        statistic_exp = 35
        x = (0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
        y = (1.15, 0.88, 0.90, 0.74, 1.21)
        res = mannwhitneyu2(x, y, use_continuity=use_continuity,
                            alternative=alternative, method=method)
        assert_equal(res.statistic, statistic_exp)
        assert_allclose(res.pvalue, pvalue_exp)

    def test_gh_6897(self):
        with assert_raises(ValueError, match="`x` and `y` must be of nonzero"):
            mannwhitneyu2([], [])

    def test_gh_4067(self):
        a=np.array([np.nan,np.nan,np.nan,np.nan,np.nan])
        b=np.array([np.nan,np.nan,np.nan,np.nan,np.nan])
        with assert_raises(ValueError, match="`x` and `y` must not contain"):
            mannwhitneyu2(a,b)
