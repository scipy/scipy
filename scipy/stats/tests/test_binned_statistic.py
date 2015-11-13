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
        cls.w = np.random.random(100)

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
        attributes = ('statistic', 'bin_edges', 'binnumbers')
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

    def test_1d_multi_values(self):
        x = self.x
        v = self.v
        w = self.w

        stat1v, edges1v, bc1v = binned_statistic(x, v, 'mean', bins=10)
        stat1w, edges1w, bc1w = binned_statistic(x, w, 'mean', bins=10)
        stat2, edges2, bc2 = binned_statistic(x, [v, w], 'mean', bins=10)

        assert_array_almost_equal(stat2[0], stat1v)
        assert_array_almost_equal(stat2[1], stat1w)
        assert_array_almost_equal(edges1v, edges2)
        assert_array_almost_equal(bc1v, bc2)



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
        attributes = ('statistic', 'x_edges', 'y_edges', 'binnumbers')
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

    def test_2d_multi_values(self):
        x = self.x
        y = self.y
        v = self.v
        w = self.w

        stat1v, binx1v, biny1v, bc1v = binned_statistic_2d(
            x, y, v, 'mean', bins=8)
        stat1w, binx1w, biny1w, bc1w = binned_statistic_2d(
            x, y, w, 'mean', bins=8)
        stat2, binx2, biny2, bc2 = binned_statistic_2d(
            x, y, [v, w], 'mean', bins=8)

        assert_array_almost_equal(stat2[0], stat1v)
        assert_array_almost_equal(stat2[1], stat1w)
        assert_array_almost_equal(binx1v, binx2)
        assert_array_almost_equal(biny1w, biny2)
        assert_array_almost_equal(bc1v, bc2)

    def test_2d_binnumbers_unraveled(self):
        x = self.x
        y = self.y
        v = self.v

        stat, edgesx, bcx = binned_statistic(x, v, 'mean', bins=10)
        stat, edgesy, bcy = binned_statistic(y, v, 'mean', bins=10)

        stat2, edgesx2, edgesy2, bc2 = binned_statistic_2d(
            x, y, v, 'mean', bins=10, expand_binnumbers=True)

        bcx3 = np.searchsorted(edgesx, x, side='right')
        bcy3 = np.searchsorted(edgesy, y, side='right')

        # `numpy.searchsorted` is non-inclusive on right-edge, compensate
        bcx3[x == x.max()] -= 1
        bcy3[y == y.max()] -= 1

        bcx4 = np.digitize(x, edgesx, right=True)
        bcy4 = np.digitize(y, edgesy, right=True)

        # `numpy.digitize` is non-inclusive on left-edge, compensate
        bcx4[x == x.min()] += 1
        bcy4[y == y.min()] += 1

        assert_array_almost_equal(bcx, bc2[0])
        assert_array_almost_equal(bcy, bc2[1])
        assert_array_almost_equal(bcx4, bc2[0])
        assert_array_almost_equal(bcy4, bc2[1])
        assert_array_almost_equal(bcx3, bc2[0])
        assert_array_almost_equal(bcy3, bc2[1])


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
        attributes = ('statistic', 'bin_edges', 'binnumbers')
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

    def test_dd_multi_values(self):
        X = self.X
        v = self.v
        w = self.w

        stat1v, edges1v, bc1v = binned_statistic_dd(X, v, np.std, bins=8)
        stat1w, edges1w, bc1w = binned_statistic_dd(X, w, np.std, bins=8)
        stat2, edges2, bc2 = binned_statistic_dd(X, [v, w], np.std, bins=8)

        assert_array_almost_equal(stat2[0], stat1v)
        assert_array_almost_equal(stat2[1], stat1w)
        assert_array_almost_equal(edges1v, edges2)
        assert_array_almost_equal(edges1w, edges2)
        assert_array_almost_equal(bc1v, bc2)

    def test_dd_binnumbers_unraveled(self):
        X = self.X
        v = self.v

        stat, edgesx, bcx = binned_statistic(X[:, 0], v, 'mean', bins=10)
        stat, edgesy, bcy = binned_statistic(X[:, 1], v, 'mean', bins=10)
        stat, edgesz, bcz = binned_statistic(X[:, 2], v, 'mean', bins=10)

        stat2, edges2, bc2 = binned_statistic_dd(
            X, v, 'mean', bins=10, expand_binnumbers=True)

        assert_array_almost_equal(bcx, bc2[0])
        assert_array_almost_equal(bcy, bc2[1])
        assert_array_almost_equal(bcz, bc2[2])

if __name__ == "__main__":
    run_module_suite()
