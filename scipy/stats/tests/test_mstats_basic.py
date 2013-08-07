"""
Tests for the stats.mstats module (support for maskd arrays)
"""
from __future__ import division, print_function, absolute_import

import warnings

import numpy as np
from numpy import nan
import numpy.ma as ma
from numpy.ma import masked, nomask

import scipy.stats.mstats as mstats
from scipy import stats
from numpy.testing import TestCase, run_module_suite
from numpy.ma.testutils import (assert_equal, assert_almost_equal,
    assert_array_almost_equal, assert_array_almost_equal_nulp, assert_)


class TestMquantiles(TestCase):
    """Regression tests for mstats module."""
    def test_mquantiles_limit_keyword(self):
        """Ticket #867"""
        data = np.array([[6., 7., 1.],
                         [47., 15., 2.],
                         [49., 36., 3.],
                         [15., 39., 4.],
                         [42., 40., -999.],
                         [41., 41., -999.],
                         [7., -999., -999.],
                         [39., -999., -999.],
                         [43., -999., -999.],
                         [40., -999., -999.],
                         [36., -999., -999.]])
        desired = [[19.2, 14.6, 1.45],
                   [40.0, 37.5, 2.5],
                   [42.8, 40.05, 3.55]]
        quants = mstats.mquantiles(data, axis=0, limit=(0, 50))
        assert_almost_equal(quants, desired)


class TestGMean(TestCase):
    def test_1D(self):
        a = (1,2,3,4)
        actual = mstats.gmean(a)
        desired = np.power(1*2*3*4,1./4.)
        assert_almost_equal(actual, desired,decimal=14)

        desired1 = mstats.gmean(a,axis=-1)
        assert_almost_equal(actual, desired1, decimal=14)
        assert_(not isinstance(desired1, ma.MaskedArray))

        a = ma.array((1,2,3,4),mask=(0,0,0,1))
        actual = mstats.gmean(a)
        desired = np.power(1*2*3,1./3.)
        assert_almost_equal(actual, desired,decimal=14)

        desired1 = mstats.gmean(a,axis=-1)
        assert_almost_equal(actual, desired1, decimal=14)

    def test_2D(self):
        a = ma.array(((1,2,3,4),(1,2,3,4),(1,2,3,4)),
                     mask=((0,0,0,0),(1,0,0,1),(0,1,1,0)))
        actual = mstats.gmean(a)
        desired = np.array((1,2,3,4))
        assert_array_almost_equal(actual, desired, decimal=14)
        #
        desired1 = mstats.gmean(a,axis=0)
        assert_array_almost_equal(actual, desired1, decimal=14)
        #
        actual = mstats.gmean(a, -1)
        desired = ma.array((np.power(1*2*3*4,1./4.),
                            np.power(2*3,1./2.),
                            np.power(1*4,1./2.)))
        assert_array_almost_equal(actual, desired, decimal=14)


class TestHMean(TestCase):
    def test_1D(self):
        a = (1,2,3,4)
        actual = mstats.hmean(a)
        desired = 4. / (1./1 + 1./2 + 1./3 + 1./4)
        assert_almost_equal(actual, desired, decimal=14)
        desired1 = mstats.hmean(ma.array(a),axis=-1)
        assert_almost_equal(actual, desired1, decimal=14)
        #
        a = ma.array((1,2,3,4),mask=(0,0,0,1))
        actual = mstats.hmean(a)
        desired = 3. / (1./1 + 1./2 + 1./3)
        assert_almost_equal(actual, desired,decimal=14)
        desired1 = mstats.hmean(a,axis=-1)
        assert_almost_equal(actual, desired1, decimal=14)

    def test_2D(self):
        a = ma.array(((1,2,3,4),(1,2,3,4),(1,2,3,4)),
                     mask=((0,0,0,0),(1,0,0,1),(0,1,1,0)))
        actual = mstats.hmean(a)
        desired = ma.array((1,2,3,4))
        assert_array_almost_equal(actual, desired, decimal=14)
        #
        actual1 = mstats.hmean(a,axis=-1)
        desired = (4./(1/1.+1/2.+1/3.+1/4.),
                   2./(1/2.+1/3.),
                   2./(1/1.+1/4.)
                   )
        assert_array_almost_equal(actual1, desired, decimal=14)


class TestRanking(TestCase):

    def __init__(self, *args, **kwargs):
        TestCase.__init__(self, *args, **kwargs)

    def test_ranking(self):
        x = ma.array([0,1,1,1,2,3,4,5,5,6,])
        assert_almost_equal(mstats.rankdata(x),[1,3,3,3,5,6,7,8.5,8.5,10])
        x[[3,4]] = masked
        assert_almost_equal(mstats.rankdata(x),[1,2.5,2.5,0,0,4,5,6.5,6.5,8])
        assert_almost_equal(mstats.rankdata(x,use_missing=True),
                            [1,2.5,2.5,4.5,4.5,4,5,6.5,6.5,8])
        x = ma.array([0,1,5,1,2,4,3,5,1,6,])
        assert_almost_equal(mstats.rankdata(x),[1,3,8.5,3,5,7,6,8.5,3,10])
        x = ma.array([[0,1,1,1,2], [3,4,5,5,6,]])
        assert_almost_equal(mstats.rankdata(x),[[1,3,3,3,5],[6,7,8.5,8.5,10]])
        assert_almost_equal(mstats.rankdata(x,axis=1),[[1,3,3,3,5],[1,2,3.5,3.5,5]])
        assert_almost_equal(mstats.rankdata(x,axis=0),[[1,1,1,1,1],[2,2,2,2,2,]])


class TestCorr(TestCase):

    def test_pearsonr(self):
        # Tests some computations of Pearson's r
        x = ma.arange(10)
        with warnings.catch_warnings():
            # The tests in this context are edge cases, with perfect
            # correlation or anticorrelation, or totally masked data.
            # None of these should trigger a RuntimeWarning.
            warnings.simplefilter("error", RuntimeWarning)

            assert_almost_equal(mstats.pearsonr(x, x)[0], 1.0)
            assert_almost_equal(mstats.pearsonr(x, x[::-1])[0], -1.0)

            x = ma.array(x, mask=True)
            pr = mstats.pearsonr(x, x)
            assert_(pr[0] is masked)
            assert_(pr[1] is masked)

        x1 = ma.array([-1.0, 0.0, 1.0])
        y1 = ma.array([0, 0, 3])
        r, p = mstats.pearsonr(x1, y1)
        assert_almost_equal(r, np.sqrt(3)/2)
        assert_almost_equal(p, 1.0/3)

        # (x2, y2) have the same unmasked data as (x1, y1).
        mask = [False, False, False, True]
        x2 = ma.array([-1.0, 0.0, 1.0, 99.0], mask=mask)
        y2 = ma.array([0, 0, 3, -1], mask=mask)
        r, p = mstats.pearsonr(x2, y2)
        assert_almost_equal(r, np.sqrt(3)/2)
        assert_almost_equal(p, 1.0/3)

    def test_spearmanr(self):
        "Tests some computations of Spearman's rho"
        (x, y) = ([5.05,6.75,3.21,2.66],[1.65,2.64,2.64,6.95])
        assert_almost_equal(mstats.spearmanr(x,y)[0], -0.6324555)
        (x, y) = ([5.05,6.75,3.21,2.66,np.nan],[1.65,2.64,2.64,6.95,np.nan])
        (x, y) = (ma.fix_invalid(x), ma.fix_invalid(y))
        assert_almost_equal(mstats.spearmanr(x,y)[0], -0.6324555)
        #
        x = [2.0, 47.4, 42.0, 10.8, 60.1, 1.7, 64.0, 63.1,
              1.0, 1.4, 7.9, 0.3, 3.9, 0.3, 6.7]
        y = [22.6, 08.3, 44.4, 11.9, 24.6, 0.6, 5.7, 41.6,
              0.0, 0.6, 6.7, 3.8, 1.0, 1.2, 1.4]
        assert_almost_equal(mstats.spearmanr(x,y)[0], 0.6887299)
        x = [2.0, 47.4, 42.0, 10.8, 60.1, 1.7, 64.0, 63.1,
              1.0, 1.4, 7.9, 0.3, 3.9, 0.3, 6.7, np.nan]
        y = [22.6, 08.3, 44.4, 11.9, 24.6, 0.6, 5.7, 41.6,
              0.0, 0.6, 6.7, 3.8, 1.0, 1.2, 1.4, np.nan]
        (x, y) = (ma.fix_invalid(x), ma.fix_invalid(y))
        assert_almost_equal(mstats.spearmanr(x,y)[0], 0.6887299)

    def test_kendalltau(self):
        "Tests some computations of Kendall's tau"
        x = ma.fix_invalid([5.05, 6.75, 3.21, 2.66,np.nan])
        y = ma.fix_invalid([1.65, 26.5, -5.93, 7.96, np.nan])
        z = ma.fix_invalid([1.65, 2.64, 2.64, 6.95, np.nan])
        assert_almost_equal(np.asarray(mstats.kendalltau(x,y)),
                            [+0.3333333,0.4969059])
        assert_almost_equal(np.asarray(mstats.kendalltau(x,z)),
                            [-0.5477226,0.2785987])
        #
        x = ma.fix_invalid([0, 0, 0, 0,20,20, 0,60, 0,20,
                            10,10, 0,40, 0,20, 0, 0, 0, 0, 0, np.nan])
        y = ma.fix_invalid([0,80,80,80,10,33,60, 0,67,27,
                            25,80,80,80,80,80,80, 0,10,45, np.nan, 0])
        result = mstats.kendalltau(x,y)
        assert_almost_equal(np.asarray(result), [-0.1585188, 0.4128009])

    def test_kendalltau_seasonal(self):
        "Tests the seasonal Kendall tau."
        x = [[nan,nan, 4, 2, 16, 26, 5, 1, 5, 1, 2, 3, 1],
             [4, 3, 5, 3, 2, 7, 3, 1, 1, 2, 3, 5, 3],
             [3, 2, 5, 6, 18, 4, 9, 1, 1,nan, 1, 1,nan],
             [nan, 6, 11, 4, 17,nan, 6, 1, 1, 2, 5, 1, 1]]
        x = ma.fix_invalid(x).T
        output = mstats.kendalltau_seasonal(x)
        assert_almost_equal(output['global p-value (indep)'], 0.008, 3)
        assert_almost_equal(output['seasonal p-value'].round(2),
                            [0.18,0.53,0.20,0.04])

    def test_pointbiserial(self):
        "Tests point biserial"
        x = [1,0,1,1,1,1,0,1,0,0,0,1,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,
             0,0,0,0,1,-1]
        y = [14.8,13.8,12.4,10.1,7.1,6.1,5.8,4.6,4.3,3.5,3.3,3.2,3.0,
             2.8,2.8,2.5,2.4,2.3,2.1,1.7,1.7,1.5,1.3,1.3,1.2,1.2,1.1,
             0.8,0.7,0.6,0.5,0.2,0.2,0.1,np.nan]
        assert_almost_equal(mstats.pointbiserialr(x, y)[0], 0.36149, 5)


class TestTrimming(TestCase):

    def test_trim(self):
        "Tests trimming"
        a = ma.arange(10)
        assert_equal(mstats.trim(a), [0,1,2,3,4,5,6,7,8,9])
        a = ma.arange(10)
        assert_equal(mstats.trim(a,(2,8)), [None,None,2,3,4,5,6,7,8,None])
        a = ma.arange(10)
        assert_equal(mstats.trim(a,limits=(2,8),inclusive=(False,False)),
                     [None,None,None,3,4,5,6,7,None,None])
        a = ma.arange(10)
        assert_equal(mstats.trim(a,limits=(0.1,0.2),relative=True),
                     [None,1,2,3,4,5,6,7,None,None])
        #
        a = ma.arange(12)
        a[[0,-1]] = a[5] = masked
        assert_equal(mstats.trim(a,(2,8)),
                     [None,None,2,3,4,None,6,7,8,None,None,None])
        #
        x = ma.arange(100).reshape(10,10)
        trimx = mstats.trim(x,(0.1,0.2),relative=True,axis=None)
        assert_equal(trimx._mask.ravel(),[1]*10+[0]*70+[1]*20)
        trimx = mstats.trim(x,(0.1,0.2),relative=True,axis=0)
        assert_equal(trimx._mask.ravel(),[1]*10+[0]*70+[1]*20)
        trimx = mstats.trim(x,(0.1,0.2),relative=True,axis=-1)
        assert_equal(trimx._mask.T.ravel(),[1]*10+[0]*70+[1]*20)
        #
        x = ma.arange(110).reshape(11,10)
        x[1] = masked
        trimx = mstats.trim(x,(0.1,0.2),relative=True,axis=None)
        assert_equal(trimx._mask.ravel(),[1]*20+[0]*70+[1]*20)
        trimx = mstats.trim(x,(0.1,0.2),relative=True,axis=0)
        assert_equal(trimx._mask.ravel(),[1]*20+[0]*70+[1]*20)
        trimx = mstats.trim(x.T,(0.1,0.2),relative=True,axis=-1)
        assert_equal(trimx.T._mask.ravel(),[1]*20+[0]*70+[1]*20)

    def test_trim_old(self):
        "Tests trimming."
        x = ma.arange(100)
        assert_equal(mstats.trimboth(x).count(), 60)
        assert_equal(mstats.trimtail(x,tail='r').count(), 80)
        x[50:70] = masked
        trimx = mstats.trimboth(x)
        assert_equal(trimx.count(), 48)
        assert_equal(trimx._mask, [1]*16 + [0]*34 + [1]*20 + [0]*14 + [1]*16)
        x._mask = nomask
        x.shape = (10,10)
        assert_equal(mstats.trimboth(x).count(), 60)
        assert_equal(mstats.trimtail(x).count(), 80)

    def test_trimmedmean(self):
        "Tests the trimmed mean."
        data = ma.array([77, 87, 88,114,151,210,219,246,253,262,
                         296,299,306,376,428,515,666,1310,2611])
        assert_almost_equal(mstats.trimmed_mean(data,0.1), 343, 0)
        assert_almost_equal(mstats.trimmed_mean(data,(0.1,0.1)), 343, 0)
        assert_almost_equal(mstats.trimmed_mean(data,(0.2,0.2)), 283, 0)

    def test_trimmed_stde(self):
        "Tests the trimmed mean standard error."
        data = ma.array([77, 87, 88,114,151,210,219,246,253,262,
                         296,299,306,376,428,515,666,1310,2611])
        assert_almost_equal(mstats.trimmed_stde(data,(0.2,0.2)), 56.13193, 5)
        assert_almost_equal(mstats.trimmed_stde(data,0.2), 56.13193, 5)

    def test_winsorization(self):
        "Tests the Winsorization of the data."
        data = ma.array([77, 87, 88,114,151,210,219,246,253,262,
                         296,299,306,376,428,515,666,1310,2611])
        assert_almost_equal(mstats.winsorize(data,(0.2,0.2)).var(ddof=1),
                            21551.4, 1)
        data[5] = masked
        winsorized = mstats.winsorize(data)
        assert_equal(winsorized.mask, data.mask)


class TestMoments(TestCase):
    # Comparison numbers are found using R v.1.5.1
    # note that length(testcase) = 4
    # testmathworks comes from documentation for the
    # Statistics Toolbox for Matlab and can be found at both
    # http://www.mathworks.com/access/helpdesk/help/toolbox/stats/kurtosis.shtml
    # http://www.mathworks.com/access/helpdesk/help/toolbox/stats/skewness.shtml
    # Note that both test cases came from here.
    testcase = [1,2,3,4]
    testmathworks = ma.fix_invalid([1.165, 0.6268, 0.0751, 0.3516, -0.6965,
                                    np.nan])
    testcase_2d = ma.array(
    np.array([[0.05245846, 0.50344235, 0.86589117, 0.36936353, 0.46961149],
           [0.11574073, 0.31299969, 0.45925772, 0.72618805, 0.75194407],
           [0.67696689, 0.91878127, 0.09769044, 0.04645137, 0.37615733],
           [0.05903624, 0.29908861, 0.34088298, 0.66216337, 0.83160998],
           [0.64619526, 0.94894632, 0.27855892, 0.0706151, 0.39962917]]),
    mask=np.array([[True, False, False, True, False],
           [True, True, True, False, True],
           [False, False, False, False, False],
           [True, True, True, True, True],
           [False, False, True, False, False]], dtype=np.bool))

    def test_moment(self):
        # mean((testcase-mean(testcase))**power,axis=0),axis=0))**power))
        y = mstats.moment(self.testcase,1)
        assert_almost_equal(y,0.0,10)
        y = mstats.moment(self.testcase,2)
        assert_almost_equal(y,1.25)
        y = mstats.moment(self.testcase,3)
        assert_almost_equal(y,0.0)
        y = mstats.moment(self.testcase,4)
        assert_almost_equal(y,2.5625)

    def test_variation(self):
        #variation = samplestd/mean
##        y = stats.variation(self.shoes[0])
##        assert_almost_equal(y,21.8770668)
        y = mstats.variation(self.testcase)
        assert_almost_equal(y,0.44721359549996, 10)

    def test_skewness(self):
        # sum((testmathworks-mean(testmathworks,axis=0))**3,axis=0)/((sqrt(var(testmathworks)*4/5))**3)/5
        y = mstats.skew(self.testmathworks)
        assert_almost_equal(y,-0.29322304336607,10)
        y = mstats.skew(self.testmathworks,bias=0)
        assert_almost_equal(y,-0.437111105023940,10)
        y = mstats.skew(self.testcase)
        assert_almost_equal(y,0.0,10)

    def test_kurtosis(self):
        #    sum((testcase-mean(testcase,axis=0))**4,axis=0)/((sqrt(var(testcase)*3/4))**4)/4
        #    sum((test2-mean(testmathworks,axis=0))**4,axis=0)/((sqrt(var(testmathworks)*4/5))**4)/5
        #    Set flags for axis = 0 and
        #    fisher=0 (Pearson's definition of kurtosis for compatibility with Matlab)
        y = mstats.kurtosis(self.testmathworks,0,fisher=0,bias=1)
        assert_almost_equal(y, 2.1658856802973,10)
        # Note that MATLAB has confusing docs for the following case
        #  kurtosis(x,0) gives an unbiased estimate of Pearson's skewness
        #  kurtosis(x)  gives a biased estimate of Fisher's skewness (Pearson-3)
        #  The MATLAB docs imply that both should give Fisher's
        y = mstats.kurtosis(self.testmathworks,fisher=0, bias=0)
        assert_almost_equal(y, 3.663542721189047,10)
        y = mstats.kurtosis(self.testcase,0,0)
        assert_almost_equal(y,1.64)

        # test that kurtosis works on multidimensional masked arrays
        correct_2d = ma.array(np.array([-1.5, -3., -1.47247052385, 0.,
                                        -1.26979517952]),
                              mask=np.array([False, False, False, True,
                                             False], dtype=np.bool))
        assert_array_almost_equal(mstats.kurtosis(self.testcase_2d, 1),
                                  correct_2d)
        for i, row in enumerate(self.testcase_2d):
            assert_almost_equal(mstats.kurtosis(row), correct_2d[i])

        correct_2d_bias_corrected = ma.array(
            np.array([-1.5, -3., -1.88988209538, 0., -0.5234638463918877]),
            mask=np.array([False, False, False, True, False], dtype=np.bool))
        assert_array_almost_equal(mstats.kurtosis(self.testcase_2d, 1,
                                                  bias=False),
                                  correct_2d_bias_corrected)
        for i, row in enumerate(self.testcase_2d):
            assert_almost_equal(mstats.kurtosis(row, bias=False),
                                correct_2d_bias_corrected[i])

        # Check consistency between stats and mstats implementations
        assert_array_almost_equal_nulp(mstats.kurtosis(self.testcase_2d[2, :]),
                                       stats.kurtosis(self.testcase_2d[2, :]))

    def test_mode(self):
        a1 = [0,0,0,1,1,1,2,3,3,3,3,4,5,6,7]
        a2 = np.reshape(a1, (3,5))
        a3 = np.array([1,2,3,4,5,6])
        a4 = np.reshape(a3, (3,2))
        ma1 = ma.masked_where(ma.array(a1) > 2, a1)
        ma2 = ma.masked_where(a2 > 2, a2)
        ma3 = ma.masked_where(a3 < 2, a3)
        ma4 = ma.masked_where(ma.array(a4) < 2, a4)
        assert_equal(mstats.mode(a1, axis=None), (3,4))
        assert_equal(mstats.mode(a1, axis=0), (3,4))
        assert_equal(mstats.mode(ma1, axis=None), (0,3))
        assert_equal(mstats.mode(a2, axis=None), (3,4))
        assert_equal(mstats.mode(ma2, axis=None), (0,3))
        assert_equal(mstats.mode(a3, axis=None), (1,1))
        assert_equal(mstats.mode(ma3, axis=None), (2,1))
        assert_equal(mstats.mode(a2, axis=0), ([[0,0,0,1,1]], [[1,1,1,1,1]]))
        assert_equal(mstats.mode(ma2, axis=0), ([[0,0,0,1,1]], [[1,1,1,1,1]]))
        assert_equal(mstats.mode(a2, axis=-1), ([[0],[3],[3]], [[3],[3],[1]]))
        assert_equal(mstats.mode(ma2, axis=-1), ([[0],[1],[0]], [[3],[1],[0]]))
        assert_equal(mstats.mode(ma4, axis=0), ([[3,2]], [[1,1]]))
        assert_equal(mstats.mode(ma4, axis=-1), ([[2],[3],[5]], [[1],[1],[1]]))


class TestPercentile(TestCase):
    def setUp(self):
        self.a1 = [3,4,5,10,-3,-5,6]
        self.a2 = [3,-6,-2,8,7,4,2,1]
        self.a3 = [3.,4,5,10,-3,-5,-6,7.0]

    def test_percentile(self):
        x = np.arange(8) * 0.5
        assert_equal(mstats.scoreatpercentile(x, 0), 0.)
        assert_equal(mstats.scoreatpercentile(x, 100), 3.5)
        assert_equal(mstats.scoreatpercentile(x, 50), 1.75)

    def test_2D(self):
        x = ma.array([[1, 1, 1],
                      [1, 1, 1],
                      [4, 4, 3],
                      [1, 1, 1],
                      [1, 1, 1]])
        assert_equal(mstats.scoreatpercentile(x,50), [1,1,1])


class TestVariability(TestCase):
    """  Comparison numbers are found using R v.1.5.1
         note that length(testcase) = 4
    """
    testcase = ma.fix_invalid([1,2,3,4,np.nan])

    def test_signaltonoise(self):
        """
        this is not in R, so used
        mean(testcase,axis=0)/(sqrt(var(testcase)*3/4)) """
        #y = stats.signaltonoise(self.shoes[0])
        #assert_approx_equal(y,4.5709967)
        y = mstats.signaltonoise(self.testcase)
        assert_almost_equal(y,2.236067977)

    def test_sem(self):
        """
        this is not in R, so used
        sqrt(var(testcase)*3/4)/sqrt(3)
        """
        #y = stats.sem(self.shoes[0])
        #assert_approx_equal(y,0.775177399)
        y = mstats.sem(self.testcase)
        assert_almost_equal(y,0.6454972244)

    def test_zmap(self):
        """
        not in R, so tested by using
        (testcase[i]-mean(testcase,axis=0))/sqrt(var(testcase)*3/4)
        """
        y = mstats.zmap(self.testcase, self.testcase)
        desired_unmaskedvals = ([-1.3416407864999, -0.44721359549996,
                                 0.44721359549996, 1.3416407864999])
        assert_array_almost_equal(desired_unmaskedvals,
                                  y.data[y.mask == False], decimal=12)

    def test_zscore(self):
        """
        not in R, so tested by using
        (testcase[i]-mean(testcase,axis=0))/sqrt(var(testcase)*3/4)
        """
        y = mstats.zscore(self.testcase)
        desired = ma.fix_invalid([-1.3416407864999, -0.44721359549996,
                                  0.44721359549996, 1.3416407864999, np.nan])
        assert_almost_equal(desired, y, decimal=12)


class TestMisc(TestCase):

    def test_obrientransform(self):
        "Tests Obrien transform"
        args = [[5]*5+[6]*11+[7]*9+[8]*3+[9]*2+[10]*2,
                [6]+[7]*2+[8]*4+[9]*9+[10]*16]
        result = [5*[3.1828]+11*[0.5591]+9*[0.0344]+3*[1.6086]+2*[5.2817]+2*[11.0538],
                  [10.4352]+2*[4.8599]+4*[1.3836]+9*[0.0061]+16*[0.7277]]
        assert_almost_equal(np.round(mstats.obrientransform(*args).T,4),
                            result,4)

    def test_kstwosamp(self):
        "Tests the Kolmogorov-Smirnov 2 samples test"
        x = [[nan,nan, 4, 2, 16, 26, 5, 1, 5, 1, 2, 3, 1],
             [4, 3, 5, 3, 2, 7, 3, 1, 1, 2, 3, 5, 3],
             [3, 2, 5, 6, 18, 4, 9, 1, 1,nan, 1, 1,nan],
             [nan, 6, 11, 4, 17,nan, 6, 1, 1, 2, 5, 1, 1]]
        x = ma.fix_invalid(x).T
        (winter,spring,summer,fall) = x.T

        assert_almost_equal(np.round(mstats.ks_twosamp(winter,spring),4),
                            (0.1818,0.9892))
        assert_almost_equal(np.round(mstats.ks_twosamp(winter,spring,'g'),4),
                            (0.1469,0.7734))
        assert_almost_equal(np.round(mstats.ks_twosamp(winter,spring,'l'),4),
                            (0.1818,0.6744))

    def test_friedmanchisq(self):
        "Tests the Friedman Chi-square test"
        # No missing values
        args = ([9.0,9.5,5.0,7.5,9.5,7.5,8.0,7.0,8.5,6.0],
                [7.0,6.5,7.0,7.5,5.0,8.0,6.0,6.5,7.0,7.0],
                [6.0,8.0,4.0,6.0,7.0,6.5,6.0,4.0,6.5,3.0])
        result = mstats.friedmanchisquare(*args)
        assert_almost_equal(result[0], 10.4737, 4)
        assert_almost_equal(result[1], 0.005317, 6)
        # Missing values
        x = [[nan,nan, 4, 2, 16, 26, 5, 1, 5, 1, 2, 3, 1],
             [4, 3, 5, 3, 2, 7, 3, 1, 1, 2, 3, 5, 3],
             [3, 2, 5, 6, 18, 4, 9, 1, 1,nan, 1, 1,nan],
             [nan, 6, 11, 4, 17,nan, 6, 1, 1, 2, 5, 1, 1]]
        x = ma.fix_invalid(x)
        result = mstats.friedmanchisquare(*x)
        assert_almost_equal(result[0], 2.0156, 4)
        assert_almost_equal(result[1], 0.5692, 4)


def test_regress_simple():
    """Regress a line with sinusoidal noise. Test for #1273."""
    x = np.linspace(0, 100, 100)
    y = 0.2 * np.linspace(0, 100, 100) + 10
    y += np.sin(np.linspace(0, 20, 100))

    slope, intercept, r_value, p_value, sterr = mstats.linregress(x, y)
    assert_almost_equal(slope, 0.19644990055858422)
    assert_almost_equal(intercept, 10.211269918932341)


def test_plotting_positions():
    """Regression test for #1256"""
    pos = mstats.plotting_positions(np.arange(3), 0, 0)
    assert_array_almost_equal(pos.data, np.array([0.25, 0.5, 0.75]))


if __name__ == "__main__":
    run_module_suite()
