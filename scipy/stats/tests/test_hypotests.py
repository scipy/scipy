from numpy.testing import (assert_, assert_equal,
                           assert_almost_equal, assert_array_almost_equal,
                           assert_array_equal, assert_approx_equal,
                           assert_allclose, assert_warns, suppress_warnings)
from pytest import raises as assert_raises
import numpy as np
import scipy.stats as stats


def assert_raises_with_match(exception_type, match, function, *args, **kwargs):
    with assert_raises(exception_type, match=match):
        function(*args, **kwargs)


class TestSomersD(object):

    def test_like_kendalltau(self):
        # All tests correspond with one in test_stats.py `test_kendalltau

        # case without ties, con-dis equal zero
        x = [5, 2, 1, 3, 6, 4, 7, 8]
        y = [5, 2, 6, 3, 1, 8, 7, 4]
        # Cross-check with result from SAS FREQ:
        expected = (0.0, 1.0)
        res = stats.somersd(x, y)
        assert_allclose(res.statistic, expected[0])
        assert_allclose(res.pvalue, expected[1])

        # case without ties, con-dis equal zero
        x = [0, 5, 2, 1, 3, 6, 4, 7, 8]
        y = [5, 2, 0, 6, 3, 1, 8, 7, 4]
        # Cross-check with result from SAS FREQ:
        expected = (0.0, 1.0)
        res = stats.somersd(x, y)
        assert_allclose(res.statistic, expected[0])
        assert_allclose(res.pvalue, expected[1])

        # case without ties, con-dis close to zero
        x = [5, 2, 1, 3, 6, 4, 7]
        y = [5, 2, 6, 3, 1, 7, 4]
        # Cross-check with result from SAS FREQ:
        expected = (-0.1429, 0.6303)
        res = stats.somersd(x, y)
        assert_allclose(res.statistic, expected[0], atol=1e-4)
        assert_allclose(res.pvalue, expected[1], atol=1e-4)

        # simple case without ties
        x = np.arange(10)
        y = np.arange(10)
        # Cross-check with result from SAS FREQ:
        expected = (1.0, 0)
        res = stats.somersd(x, y)
        assert_allclose(res.statistic, expected[0], atol=1e-4)
        assert_allclose(res.pvalue, expected[1], atol=1e-4)

        # swap a couple values and a couple more
        x = np.arange(10)
        y = np.array([0, 2, 1, 3, 4, 6, 5, 7, 8, 9])
        # Cross-check with result from SAS FREQ:
        expected = (0.9111, 0)
        res = stats.somersd(x, y)
        assert_allclose(res.statistic, expected[0], atol=1e-4)
        assert_allclose(res.pvalue, expected[1], atol=1e-4)

        # same in opposite direction
        x = np.arange(10)
        y = np.arange(10)[::-1]
        # Cross-check with result from SAS FREQ:
        expected = (-1.0, 5.511463844797e-07)
        res = stats.somersd(x, y)
        assert_allclose(res.statistic, expected[0], atol=1e-4)
        assert_allclose(res.pvalue, expected[1], atol=1e-4)

        # swap a couple values and a couple more
        x = np.arange(10)
        y = np.array([9, 7, 8, 6, 5, 3, 4, 2, 1, 0])
        expected = (-0.9111111111111111, 2.976190476190e-05)
        res = stats.somersd(x, y)
        assert_allclose(res.statistic, expected[0], atol=1e-4)
        assert_allclose(res.pvalue, expected[1], atol=1e-4)

        # with some ties
        # Cross-check with result from SAS FREQ:
        x1 = [12, 2, 1, 12, 2]
        x2 = [1, 4, 7, 1, 0]
        expected = (-0.5, 0.3049017881787882)
        res = stats.somersd(x1, x2)
        assert_allclose(res.statistic, expected[0])
        assert_allclose(res.pvalue, expected[1])

        # with only ties in one or both inputs
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
        d_cr = 0.2727
        d_rc = 0.3429
        p = 0.0929  # same p-value for either direction
        res = stats.somersd(x, y)
        assert_allclose(res.statistic, d_cr, atol=1e-4)
        assert_allclose(res.pvalue, p, atol=1e-4)
        assert_equal(res.table.shape, (3, 2))
        res = stats.somersd(y, x)
        assert_allclose(res.statistic, d_rc, atol=1e-4)
        assert_allclose(res.pvalue, p, atol=1e-4)
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
        assert_allclose(res.statistic, -0.1169811320754717)
        assert_allclose(res.statistic, res2.statistic)
        assert_allclose(res.statistic, res3.statistic)
        assert_allclose(res.statistic, res4.statistic)

        assert_allclose(res.pvalue, 0.15637644818814952)
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
        assert_raises_with_match(ValueError, "All elements of the contingency "
                                 "table must be non-negative", stats.somersd,
                                 s5)

        s6 = s + 0.01
        assert_raises_with_match(ValueError, "All elements of the contingency "
                                 "table must be integer", stats.somersd, s6)

        assert_raises_with_match(ValueError, "At least two elements of the "
                                 "contingency table must be nonzero.",
                                 stats.somersd, [[]])

        assert_raises_with_match(ValueError, "At least two elements of the "
                                 "contingency table must be nonzero.",
                                 stats.somersd, [[1]])

        s7 = np.zeros((3, 3))
        assert_raises_with_match(ValueError, "At least two elements of the "
                                 "contingency table must be nonzero.",
                                 stats.somersd, s7)

        s7[0, 1] = 1
        assert_raises_with_match(ValueError, "At least two elements of the "
                                 "contingency table must be nonzero.",
                                 stats.somersd, s7)

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
