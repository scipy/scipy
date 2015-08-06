from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_array_almost_equal, run_module_suite
from scipy.stats import (binned_statistic, binned_statistic_2d,
                         binned_statistic_dd)

from scipy._lib.six import u
from common_tests import check_named_results


class TestBinnedStatistic(object):

    @classmethod
    def setup_class(cls):
        np.random.seed(9865)
        cls.x = np.random.random(100)
        cls.y = np.random.random(100)
        cls.v = np.random.random(100)
        cls.X = np.random.random((100, 3))

    def test_1d_count(self):
        x = self.x
        v = self.v

        count1, edges1, bc = binned_statistic(x, v, 'count', bins=10)
        count2, edges2 = np.histogram(x, bins=10)

        assert_array_almost_equal(count1, count2)
        assert_array_almost_equal(edges1, edges2)

    def test_1d_result_attributes(self):
        x = self.x
        v = self.v

        res = binned_statistic(x, v, 'count', bins=10)
        attributes = ('statistic', 'bin_edges', 'binnumber')
        check_named_results(res, attributes)

    def test_1d_sum(self):
        x = self.x
        v = self.v

        sum1, edges1, bc = binned_statistic(x, v, 'sum', bins=10)
        sum2, edges2 = np.histogram(x, bins=10, weights=v)

        assert_array_almost_equal(sum1, sum2)
        assert_array_almost_equal(edges1, edges2)

    def test_1d_mean(self):
        x = self.x
        v = self.v

        stat1, edges1, bc = binned_statistic(x, v, 'mean', bins=10)
        stat2, edges2, bc = binned_statistic(x, v, np.mean, bins=10)

        assert_array_almost_equal(stat1, stat2)
        assert_array_almost_equal(edges1, edges2)

    def test_1d_std(self):
        x = self.x
        v = self.v

        stat1, edges1, bc = binned_statistic(x, v, 'std', bins=10)
        stat2, edges2, bc = binned_statistic(x, v, np.std, bins=10)

        assert_array_almost_equal(stat1, stat2)
        assert_array_almost_equal(edges1, edges2)

    def test_1d_median(self):
        x = self.x
        v = self.v

        stat1, edges1, bc = binned_statistic(x, v, 'median', bins=10)
        stat2, edges2, bc = binned_statistic(x, v, np.median, bins=10)

        assert_array_almost_equal(stat1, stat2)
        assert_array_almost_equal(edges1, edges2)

    def test_1d_max(self):
        x = self.x
        v = self.v

        stat1, edges1, bc = binned_statistic(x, v, 'max', bins=10)
        stat2, edges2, bc = binned_statistic(x, v, np.max, bins=10)

        assert_array_almost_equal(stat1, stat2)
        assert_array_almost_equal(edges1, edges2)

    def test_1d_min(self):
        x = self.x
        v = self.v

        stat1, edges1, bc = binned_statistic(x, v, 'min', bins=10)
        stat2, edges2, bc = binned_statistic(x, v, np.min, bins=10)

        assert_array_almost_equal(stat1, stat2)
        assert_array_almost_equal(edges1, edges2)

    def test_1d_bincode(self):
        x = self.x[:20]
        v = self.v[:20]

        count1, edges1, bc = binned_statistic(x, v, 'count', bins=3)
        bc2 = np.array([3, 2, 1, 3, 2, 3, 3, 3, 3, 1, 1, 3, 3, 1, 2, 3, 1,
                        1, 2, 1])

        bcount = [(bc == i).sum() for i in np.unique(bc)]

        assert_array_almost_equal(bc, bc2)
        assert_array_almost_equal(bcount, count1)

    def test_1d_range_keyword(self):
        # Regression test for gh-3063, range can be (min, max) or [(min, max)]
        np.random.seed(9865)
        x = np.arange(30)
        data = np.random.random(30)

        mean, bins, _ = binned_statistic(x[:15], data[:15])
        mean_range, bins_range, _ = binned_statistic(x, data, range=[(0, 14)])
        mean_range2, bins_range2, _ = binned_statistic(x, data, range=(0, 14))

        assert_array_almost_equal(mean, mean_range)
        assert_array_almost_equal(bins, bins_range)
        assert_array_almost_equal(mean, mean_range2)
        assert_array_almost_equal(bins, bins_range2)

    def test_2d_count(self):
        x = self.x
        y = self.y
        v = self.v

        count1, binx1, biny1, bc = binned_statistic_2d(x, y, v, 'count', bins=5)
        count2, binx2, biny2 = np.histogram2d(x, y, bins=5)

        assert_array_almost_equal(count1, count2)
        assert_array_almost_equal(binx1, binx2)
        assert_array_almost_equal(biny1, biny2)

    def test_2d_result_attributes(self):
        x = self.x
        y = self.y
        v = self.v

        res = binned_statistic_2d(x, y, v, 'count', bins=5)
        attributes = ('statistic', 'x_edge', 'y_edge', 'binnumber')
        check_named_results(res, attributes)

    def test_2d_sum(self):
        x = self.x
        y = self.y
        v = self.v

        sum1, binx1, biny1, bc = binned_statistic_2d(x, y, v, 'sum', bins=5)
        sum2, binx2, biny2 = np.histogram2d(x, y, bins=5, weights=v)

        assert_array_almost_equal(sum1, sum2)
        assert_array_almost_equal(binx1, binx2)
        assert_array_almost_equal(biny1, biny2)

    def test_2d_mean(self):
        x = self.x
        y = self.y
        v = self.v

        stat1, binx1, biny1, bc = binned_statistic_2d(x, y, v, 'mean', bins=5)
        stat2, binx2, biny2, bc = binned_statistic_2d(x, y, v, np.mean, bins=5)

        assert_array_almost_equal(stat1, stat2)
        assert_array_almost_equal(binx1, binx2)
        assert_array_almost_equal(biny1, biny2)

    def test_2d_mean_unicode(self):
        x = self.x
        y = self.y
        v = self.v
        stat1, binx1, biny1, bc = binned_statistic_2d(x, y, v, u('mean'), bins=5)
        stat2, binx2, biny2, bc = binned_statistic_2d(x, y, v, np.mean, bins=5)
        assert_array_almost_equal(stat1, stat2)
        assert_array_almost_equal(binx1, binx2)
        assert_array_almost_equal(biny1, biny2)

    def test_2d_std(self):
        x = self.x
        y = self.y
        v = self.v

        stat1, binx1, biny1, bc = binned_statistic_2d(x, y, v, 'std', bins=5)
        stat2, binx2, biny2, bc = binned_statistic_2d(x, y, v, np.std, bins=5)

        assert_array_almost_equal(stat1, stat2)
        assert_array_almost_equal(binx1, binx2)
        assert_array_almost_equal(biny1, biny2)

    def test_2d_median(self):
        x = self.x
        y = self.y
        v = self.v

        stat1, binx1, biny1, bc = binned_statistic_2d(x, y, v, 'median', bins=5)
        stat2, binx2, biny2, bc = binned_statistic_2d(x, y, v, np.median, bins=5)

        assert_array_almost_equal(stat1, stat2)
        assert_array_almost_equal(binx1, binx2)
        assert_array_almost_equal(biny1, biny2)

    def test_2d_max(self):
        x = self.x
        y = self.y
        v = self.v

        stat1, binx1, biny1, bc = binned_statistic_2d(x, y, v, 'max', bins=5)
        stat2, binx2, biny2, bc = binned_statistic_2d(x, y, v, np.max, bins=5)

        assert_array_almost_equal(stat1, stat2)
        assert_array_almost_equal(binx1, binx2)
        assert_array_almost_equal(biny1, biny2)

    def test_2d_min(self):
        x = self.x
        y = self.y
        v = self.v

        stat1, binx1, biny1, bc = binned_statistic_2d(x, y, v, 'min', bins=5)
        stat2, binx2, biny2, bc = binned_statistic_2d(x, y, v, np.min, bins=5)

        assert_array_almost_equal(stat1, stat2)
        assert_array_almost_equal(binx1, binx2)
        assert_array_almost_equal(biny1, biny2)

    def test_2d_bincode(self):
        x = self.x[:20]
        y = self.y[:20]
        v = self.v[:20]

        count1, binx1, biny1, bc = binned_statistic_2d(x, y, v, 'count', bins=3)
        bc2 = np.array([17, 11, 6, 16, 11, 17, 18, 17, 17, 7, 6, 18, 16,
                        6, 11, 16, 6, 6, 11, 8])

        bcount = [(bc == i).sum() for i in np.unique(bc)]

        assert_array_almost_equal(bc, bc2)
        count1adj = count1[count1.nonzero()]
        assert_array_almost_equal(bcount, count1adj)

    def test_dd_count(self):
        X = self.X
        v = self.v

        count1, edges1, bc = binned_statistic_dd(X, v, 'count', bins=3)
        count2, edges2 = np.histogramdd(X, bins=3)

        assert_array_almost_equal(count1, count2)
        assert_array_almost_equal(edges1, edges2)

    def test_dd_result_attributes(self):
        X = self.X
        v = self.v

        res = binned_statistic_dd(X, v, 'count', bins=3)
        attributes = ('statistic', 'bin_edges', 'binnumber')
        check_named_results(res, attributes)

    def test_dd_sum(self):
        X = self.X
        v = self.v

        sum1, edges1, bc = binned_statistic_dd(X, v, 'sum', bins=3)
        sum2, edges2 = np.histogramdd(X, bins=3, weights=v)

        assert_array_almost_equal(sum1, sum2)
        assert_array_almost_equal(edges1, edges2)

    def test_dd_mean(self):
        X = self.X
        v = self.v

        stat1, edges1, bc = binned_statistic_dd(X, v, 'mean', bins=3)
        stat2, edges2, bc = binned_statistic_dd(X, v, np.mean, bins=3)

        assert_array_almost_equal(stat1, stat2)
        assert_array_almost_equal(edges1, edges2)

    def test_dd_std(self):
        X = self.X
        v = self.v

        stat1, edges1, bc = binned_statistic_dd(X, v, 'std', bins=3)
        stat2, edges2, bc = binned_statistic_dd(X, v, np.std, bins=3)

        assert_array_almost_equal(stat1, stat2)
        assert_array_almost_equal(edges1, edges2)

    def test_dd_median(self):
        X = self.X
        v = self.v

        stat1, edges1, bc = binned_statistic_dd(X, v, 'median', bins=3)
        stat2, edges2, bc = binned_statistic_dd(X, v, np.median, bins=3)

        assert_array_almost_equal(stat1, stat2)
        assert_array_almost_equal(edges1, edges2)

    def test_dd_max(self):
        X = self.X
        v = self.v

        stat1, edges1, bc = binned_statistic_dd(X, v, 'max', bins=3)
        stat2, edges2, bc = binned_statistic_dd(X, v, np.max, bins=3)

        assert_array_almost_equal(stat1, stat2)
        assert_array_almost_equal(edges1, edges2)

    def test_dd_min(self):
        X = self.X
        v = self.v

        stat1, edges1, bc = binned_statistic_dd(X, v, 'min', bins=3)
        stat2, edges2, bc = binned_statistic_dd(X, v, np.min, bins=3)

        assert_array_almost_equal(stat1, stat2)
        assert_array_almost_equal(edges1, edges2)

    def test_dd_bincode(self):
        X = self.X[:20]
        v = self.v[:20]

        count1, edges1, bc = binned_statistic_dd(X, v, 'count', bins=3)
        bc2 = np.array([63, 33, 86, 83, 88, 67, 57, 33, 42, 41, 82, 83, 92,
                        32, 36, 91, 43, 87, 81, 81])

        bcount = [(bc == i).sum() for i in np.unique(bc)]

        assert_array_almost_equal(bc, bc2)
        count1adj = count1[count1.nonzero()]
        assert_array_almost_equal(bcount, count1adj)


if __name__ == "__main__":
    run_module_suite()
