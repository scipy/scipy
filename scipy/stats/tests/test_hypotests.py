from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (assert_, assert_equal, assert_allclose,
                           assert_almost_equal)  # avoid new uses

import pytest
from pytest import raises as assert_raises
from scipy.stats._hypotests import (epps_singleton_2samp, cramervonmises,
                                    _cdf_cvm, cramervonmises_2samp,
                                    _pval_cvm_2samp_exact, barnard_exact)
import scipy.stats as stats
from scipy.stats import distributions
from .common_tests import check_named_results


class TestEppsSingleton:
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


class TestCvm:
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


class TestSomersD:

    def test_like_kendalltau(self):
        # All tests correspond with one in test_stats.py `test_kendalltau`

        # case without ties, con-dis equal zero
        x = [5, 2, 1, 3, 6, 4, 7, 8]
        y = [5, 2, 6, 3, 1, 8, 7, 4]
        # Cross-check with result from SAS FREQ:
        expected = (0.000000000000000, 1.000000000000000)
        res = stats.somersd(x, y)
        assert_allclose(res.statistic, expected[0], atol=1e-15)
        assert_allclose(res.pvalue, expected[1], atol=1e-15)

        # case without ties, con-dis equal zero
        x = [0, 5, 2, 1, 3, 6, 4, 7, 8]
        y = [5, 2, 0, 6, 3, 1, 8, 7, 4]
        # Cross-check with result from SAS FREQ:
        expected = (0.000000000000000, 1.000000000000000)
        res = stats.somersd(x, y)
        assert_allclose(res.statistic, expected[0], atol=1e-15)
        assert_allclose(res.pvalue, expected[1], atol=1e-15)

        # case without ties, con-dis close to zero
        x = [5, 2, 1, 3, 6, 4, 7]
        y = [5, 2, 6, 3, 1, 7, 4]
        # Cross-check with result from SAS FREQ:
        expected = (-0.142857142857140, 0.630326953157670)
        res = stats.somersd(x, y)
        assert_allclose(res.statistic, expected[0], atol=1e-15)
        assert_allclose(res.pvalue, expected[1], atol=1e-15)

        # simple case without ties
        x = np.arange(10)
        y = np.arange(10)
        # Cross-check with result from SAS FREQ:
        # SAS p value is not provided.
        expected = (1.000000000000000, 0)
        res = stats.somersd(x, y)
        assert_allclose(res.statistic, expected[0], atol=1e-15)
        assert_allclose(res.pvalue, expected[1], atol=1e-15)

        # swap a couple values and a couple more
        x = np.arange(10)
        y = np.array([0, 2, 1, 3, 4, 6, 5, 7, 8, 9])
        # Cross-check with result from SAS FREQ:
        expected = (0.911111111111110, 0.000000000000000)
        res = stats.somersd(x, y)
        assert_allclose(res.statistic, expected[0], atol=1e-15)
        assert_allclose(res.pvalue, expected[1], atol=1e-15)

        # same in opposite direction
        x = np.arange(10)
        y = np.arange(10)[::-1]
        # Cross-check with result from SAS FREQ:
        # SAS p value is not provided.
        expected = (-1.000000000000000, 0)
        res = stats.somersd(x, y)
        assert_allclose(res.statistic, expected[0], atol=1e-15)
        assert_allclose(res.pvalue, expected[1], atol=1e-15)

        # swap a couple values and a couple more
        x = np.arange(10)
        y = np.array([9, 7, 8, 6, 5, 3, 4, 2, 1, 0])
        # Cross-check with result from SAS FREQ:
        expected = (-0.9111111111111111, 0.000000000000000)
        res = stats.somersd(x, y)
        assert_allclose(res.statistic, expected[0], atol=1e-15)
        assert_allclose(res.pvalue, expected[1], atol=1e-15)

        # with some ties
        x1 = [12, 2, 1, 12, 2]
        x2 = [1, 4, 7, 1, 0]
        # Cross-check with result from SAS FREQ:
        expected = (-0.500000000000000, 0.304901788178780)
        res = stats.somersd(x1, x2)
        assert_allclose(res.statistic, expected[0], atol=1e-15)
        assert_allclose(res.pvalue, expected[1], atol=1e-15)

        # with only ties in one or both inputs
        # SAS will not produce an output for these:
        # NOTE: No statistics are computed for x * y because x has fewer
        # than 2 nonmissing levels.
        # WARNING: No OUTPUT data set is produced for this table because a
        # row or column variable has fewer than 2 nonmissing levels and no
        # statistics are computed.

        res = stats.somersd([2, 2, 2], [2, 2, 2])
        assert_allclose(res.statistic, np.nan)
        assert_allclose(res.pvalue, np.nan)

        res = stats.somersd([2, 0, 2], [2, 2, 2])
        assert_allclose(res.statistic, np.nan)
        assert_allclose(res.pvalue, np.nan)

        res = stats.somersd([2, 2, 2], [2, 0, 2])
        assert_allclose(res.statistic, np.nan)
        assert_allclose(res.pvalue, np.nan)

        res = stats.somersd([0], [0])
        assert_allclose(res.statistic, np.nan)
        assert_allclose(res.pvalue, np.nan)

        # empty arrays provided as input
        res = stats.somersd([], [])
        assert_allclose(res.statistic, np.nan)
        assert_allclose(res.pvalue, np.nan)

        # test unequal length inputs
        x = np.arange(10.)
        y = np.arange(20.)
        assert_raises(ValueError, stats.somersd, x, y)

    def test_asymmetry(self):
        # test that somersd is asymmetric w.r.t. input order and that
        # convention is as described: first input is row variable & independent
        # data is from Wikipedia:
        # https://en.wikipedia.org/wiki/Somers%27_D
        # but currently that example contradicts itself - it says X is
        # independent yet take D_XY

        x = [1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 1, 2,
             2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3]
        y = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2,
             2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        # Cross-check with result from SAS FREQ:
        d_cr = 0.272727272727270
        d_rc = 0.342857142857140
        p = 0.092891940883700  # same p-value for either direction
        res = stats.somersd(x, y)
        assert_allclose(res.statistic, d_cr, atol=1e-15)
        assert_allclose(res.pvalue, p, atol=1e-4)
        assert_equal(res.table.shape, (3, 2))
        res = stats.somersd(y, x)
        assert_allclose(res.statistic, d_rc, atol=1e-15)
        assert_allclose(res.pvalue, p, atol=1e-15)
        assert_equal(res.table.shape, (2, 3))

    def test_somers_original(self):
        # test against Somers' original paper [1]

        # Table 5A
        # Somers' convention was column IV
        table = np.array([[8, 2], [6, 5], [3, 4], [1, 3], [2, 3]])
        # Our convention (and that of SAS FREQ) is row IV
        table = table.T
        dyx = 129/340
        assert_allclose(stats.somersd(table).statistic, dyx)

        # table 7A - d_yx = 1
        table = np.array([[25, 0], [85, 0], [0, 30]])
        dxy, dyx = 3300/5425, 3300/3300
        assert_allclose(stats.somersd(table).statistic, dxy)
        assert_allclose(stats.somersd(table.T).statistic, dyx)

        # table 7B - d_yx < 0
        table = np.array([[25, 0], [0, 30], [85, 0]])
        dyx = -1800/3300
        assert_allclose(stats.somersd(table.T).statistic, dyx)

    def test_contingency_table_with_zero_rows_cols(self):
        # test that zero rows/cols in contingency table don't affect result

        N = 100
        shape = 4, 6
        size = np.prod(shape)

        np.random.seed(0)
        s = stats.multinomial.rvs(N, p=np.ones(size)/size).reshape(shape)
        res = stats.somersd(s)

        s2 = np.insert(s, 2, np.zeros(shape[1]), axis=0)
        res2 = stats.somersd(s2)

        s3 = np.insert(s, 2, np.zeros(shape[0]), axis=1)
        res3 = stats.somersd(s3)

        s4 = np.insert(s2, 2, np.zeros(shape[0]+1), axis=1)
        res4 = stats.somersd(s4)

        # Cross-check with result from SAS FREQ:
        assert_allclose(res.statistic, -0.116981132075470, atol=1e-15)
        assert_allclose(res.statistic, res2.statistic)
        assert_allclose(res.statistic, res3.statistic)
        assert_allclose(res.statistic, res4.statistic)

        assert_allclose(res.pvalue, 0.156376448188150, atol=1e-15)
        assert_allclose(res.pvalue, res2.pvalue)
        assert_allclose(res.pvalue, res3.pvalue)
        assert_allclose(res.pvalue, res4.pvalue)

    def test_invalid_contingency_tables(self):
        N = 100
        shape = 4, 6
        size = np.prod(shape)

        np.random.seed(0)
        # start with a valid contingency table
        s = stats.multinomial.rvs(N, p=np.ones(size)/size).reshape(shape)

        s5 = s - 2
        message = "All elements of the contingency table must be non-negative"
        with assert_raises(ValueError, match=message):
            stats.somersd(s5)

        s6 = s + 0.01
        message = "All elements of the contingency table must be integer"
        with assert_raises(ValueError, match=message):
            stats.somersd(s6)

        message = ("At least two elements of the contingency "
                   "table must be nonzero.")
        with assert_raises(ValueError, match=message):
            stats.somersd([[]])

        with assert_raises(ValueError, match=message):
            stats.somersd([[1]])

        s7 = np.zeros((3, 3))
        with assert_raises(ValueError, match=message):
            stats.somersd(s7)

        s7[0, 1] = 1
        with assert_raises(ValueError, match=message):
            stats.somersd(s7)

    def test_only_ranks_matter(self):
        # only ranks of input data should matter
        x = [1, 2, 3]
        x2 = [-1, 2.1, np.inf]
        y = [3, 2, 1]
        y2 = [0, -0.5, -np.inf]
        res = stats.somersd(x, y)
        res2 = stats.somersd(x2, y2)
        assert_equal(res.statistic, res2.statistic)
        assert_equal(res.pvalue, res2.pvalue)

    def test_contingency_table_return(self):
        # check that contingency table is returned
        x = np.arange(10)
        y = np.arange(10)
        res = stats.somersd(x, y)
        assert_equal(res.table, np.eye(10))


class TestBarnardExact:
    """Some tests to show that barnard_exact() works correctly."""

    @pytest.mark.parametrize(
        "input_sample,expected",
        [
            ([[43, 40], [10, 39]], (3.555406779643, 0.000362832367)),
            ([[100, 2], [1000, 5]], (-1.776382925679, 0.135126970878)),
            ([[2, 7], [8, 2]], (-2.518474945157, 0.019210815430)),
            ([[5, 1], [10, 10]], (1.449486150679, 0.156277546306)),
            ([[5, 15], [20, 20]], (-1.851640199545, 0.066363501421)),
            ([[5, 16], [20, 25]], (-1.609639949352, 0.116984852192)),
            ([[10, 5], [10, 1]], (-1.449486150679, 0.177536588915)),
            ([[5, 0], [1, 4]], (2.581988897472, 0.013671875000)),
            ([[0, 1], [3, 2]], (-1.095445115010, 0.509667991877)),
            ([[0, 2], [6, 4]], (-1.549193338483, 0.197019618792)),
            ([[2, 7], [8, 2]], (-2.518474945157, 0.019210815430)),
        ],
    )
    def test_precise(self, input_sample, expected):
        """The expected values have been generated by R, using a resolution
        for the nuisance parameter of 1e-6 :
        ```R
        library(Barnard)
        options(digits=10)
        barnard.test(43, 40, 10, 39, dp=1e-6, pooled=TRUE)
        ```
        """
        res = barnard_exact(input_sample)
        statistic, pvalue = res.statistic, res.pvalue
        assert_allclose([statistic, pvalue], expected)

    @pytest.mark.parametrize(
        "input_sample,expected",
        [
            ([[43, 40], [10, 39]], (3.920362887717, 0.000289470662)),
            ([[100, 2], [1000, 5]], (-1.139432816087, 0.950272080594)),
            ([[2, 7], [8, 2]], (-3.079373904042, 0.020172119141)),
            ([[5, 1], [10, 10]], (1.622375939458, 0.150599922226)),
            ([[5, 15], [20, 20]], (-1.974771239528, 0.063038448651)),
            ([[5, 16], [20, 25]], (-1.722122973346, 0.133329494287)),
            ([[10, 5], [10, 1]], (-1.765469659009, 0.250566655215)),
            ([[5, 0], [1, 4]], (5.477225575052, 0.007812500000)),
            ([[0, 1], [3, 2]], (-1.224744871392, 0.509667991877)),
            ([[0, 2], [6, 4]], (-1.732050807569, 0.197019618792)),
            ([[2, 7], [8, 2]], (-3.079373904042, 0.020172119141)),
        ],
    )
    def test_pooled_param(self, input_sample, expected):
        """The expected values have been generated by R, using a resolution
        for the nuisance parameter of 1e-6 :
        ```R
        library(Barnard)
        options(digits=10)
        barnard.test(43, 40, 10, 39, dp=1e-6, pooled=FALSE)
        ```
        """
        res = barnard_exact(input_sample, pooled=False)
        statistic, pvalue = res.statistic, res.pvalue
        assert_allclose([statistic, pvalue], expected)

    def test_raises(self):
        # test we raise an error for wrong input number of nuisances.
        error_msg = (
            "Number of points `n` must be strictly positive, found 0"
        )
        with assert_raises(ValueError, match=error_msg):
            barnard_exact([[1, 2], [3, 4]], n=0)

        # test we raise an error for wrong shape of input.
        error_msg = "The input `table` must be of shape \\(2, 2\\)."
        with assert_raises(ValueError, match=error_msg):
            barnard_exact(np.arange(6).reshape(2, 3))

        # Test all values must be positives
        error_msg = "All values in `table` must be nonnegative."
        with assert_raises(ValueError, match=error_msg):
            barnard_exact([[-1, 2], [3, 4]])

        # Test value error on wrong alternative param
        error_msg = (
            "`alternative` should be one of {'two-sided', 'less', 'greater'},"
            " found .*"
        )
        with assert_raises(ValueError, match=error_msg):
            barnard_exact([[1, 2], [3, 4]], "not-correct")

    @pytest.mark.parametrize(
        "input_sample,expected",
        [
            ([[0, 0], [4, 3]], (1.0, 0)),
        ],
    )
    def test_edge_cases(self, input_sample, expected):
        res = barnard_exact(input_sample)
        statistic, pvalue = res.statistic, res.pvalue
        assert_equal(pvalue, expected[0])
        assert_equal(statistic, expected[1])

    @pytest.mark.parametrize(
        "input_sample,expected",
        [
            ([[0, 5], [0, 10]], (1.0, np.nan)),
            ([[5, 0], [10, 0]], (1.0, np.nan)),
        ],
    )
    def test_row_or_col_zero(self, input_sample, expected):
        res = barnard_exact(input_sample)
        statistic, pvalue = res.statistic, res.pvalue
        assert_equal(pvalue, expected[0])
        assert_equal(statistic, expected[1])

    @pytest.mark.parametrize(
        "input_sample,expected",
        [
            ([[2, 7], [8, 2]], (-2.518474945157, 0.009886140845)),
            ([[7, 200], [300, 8]], (-21.320036698460, 0.0)),
            ([[21, 28], [1957, 6]], (-30.489638143953, 0.0)),
        ],
    )
    @pytest.mark.parametrize("alternative", ["greater", "less"])
    def test_less_greater(self, input_sample, expected, alternative):
        """
        "The expected values have been generated by R, using a resolution
        for the nuisance parameter of 1e-6 :
        ```R
        library(Barnard)
        options(digits=10)
        a = barnard.test(2, 7, 8, 2, dp=1e-6, pooled=TRUE)
        a$p.value[1]
        ```
        In this test, we are using the "one-sided" return value `a$p.value[1]`
        to test our pvalue.
        """
        expected_stat, less_pvalue_expect = expected

        if alternative == "greater":
            input_sample = np.array(input_sample)[:, ::-1]
            expected_stat = -expected_stat

        res = barnard_exact(input_sample, alternative=alternative)
        statistic, pvalue = res.statistic, res.pvalue
        assert_allclose(
            [statistic, pvalue], [expected_stat, less_pvalue_expect], atol=1e-7
        )


class TestCvm_2samp:
    def test_invalid_input(self):
        x = np.arange(10).reshape((2, 5))
        y = np.arange(5)
        msg = 'The samples must be one-dimensional'
        with pytest.raises(ValueError, match=msg):
            cramervonmises_2samp(x, y)
        with pytest.raises(ValueError, match=msg):
            cramervonmises_2samp(y, x)
        msg = 'x and y must contain at least two observations.'
        with pytest.raises(ValueError, match=msg):
            cramervonmises_2samp([], y)
        with pytest.raises(ValueError, match=msg):
            cramervonmises_2samp(y, [1])
        msg = 'method must be either auto, exact or asymptotic'
        with pytest.raises(ValueError, match=msg):
            cramervonmises_2samp(y, y, 'xyz')

    def test_list_input(self):
        x = [2, 3, 4, 7, 6]
        y = [0.2, 0.7, 12, 18]
        r1 = cramervonmises_2samp(x, y)
        r2 = cramervonmises_2samp(np.array(x), np.array(y))
        assert_equal((r1.statistic, r1.pvalue), (r2.statistic, r2.pvalue))

    def test_example_conover(self):
        # Example 2 in Section 6.2 of W.J. Conover: Practical Nonparametric
        # Statistics, 1971.
        x = [7.6, 8.4, 8.6, 8.7, 9.3, 9.9, 10.1, 10.6, 11.2]
        y = [5.2, 5.7, 5.9, 6.5, 6.8, 8.2, 9.1, 9.8, 10.8, 11.3, 11.5, 12.3,
             12.5, 13.4, 14.6]
        r = cramervonmises_2samp(x, y)
        assert_allclose(r.statistic, 0.262, atol=1e-3)
        assert_allclose(r.pvalue, 0.18, atol=1e-2)

    @pytest.mark.parametrize('statistic, m, n, pval',
                             [(710, 5, 6, 48./462),
                              (1897, 7, 7, 117./1716),
                              (576, 4, 6, 2./210),
                              (1764, 6, 7, 2./1716)])
    def test_exact_pvalue(self, statistic, m, n, pval):
        # the exact values are taken from Anderson: On the distribution of the
        # two-sample Cramer-von-Mises criterion, 1962.
        # The values are taken from Table 2, 3, 4 and 5
        assert_equal(_pval_cvm_2samp_exact(statistic, m, n), pval)

    def test_large_sample(self):
        # for large samples, the statistic U gets very large
        # do a sanity check that p-value is not 0, 1 or nan
        np.random.seed(4367)
        x = distributions.norm.rvs(size=1000000)
        y = distributions.norm.rvs(size=900000)
        r = cramervonmises_2samp(x, y)
        assert_(0 < r.pvalue < 1)
        r = cramervonmises_2samp(x, y+0.1)
        assert_(0 < r.pvalue < 1)

    def test_exact_vs_asymptotic(self):
        np.random.seed(0)
        x = np.random.rand(7)
        y = np.random.rand(8)
        r1 = cramervonmises_2samp(x, y, method='exact')
        r2 = cramervonmises_2samp(x, y, method='asymptotic')
        assert_equal(r1.statistic, r2.statistic)
        assert_allclose(r1.pvalue, r2.pvalue, atol=1e-2)

    #@pytest.mark.slow
    def test_method_auto(self):
        x = np.arange(10)
        y = [0.5, 4.7, 13.1]
        r1 = cramervonmises_2samp(x, y, method='exact')
        r2 = cramervonmises_2samp(x, y, method='auto')
        assert_equal(r1.pvalue, r2.pvalue)
        # switch to asymptotic if one sample has more than 10 observations
        x = np.arange(11)
        r1 = cramervonmises_2samp(x, y, method='asymptotic')
        r2 = cramervonmises_2samp(x, y, method='auto')
        assert_equal(r1.pvalue, r2.pvalue)

    def test_same_input(self):
        # make sure trivial edge case can be handled
        # note that _cdf_cvm_inf(0) = nan. implementation avoids nan by
        # returning pvalue=1 for very small values of the statistic
        x = np.arange(15)
        res = cramervonmises_2samp(x, x)
        assert_equal((res.statistic, res.pvalue), (0.0, 1.0))
        # check exact p-value
        res = cramervonmises_2samp(x[:4], x[:4])
        assert_equal((res.statistic, res.pvalue), (0.0, 1.0))
