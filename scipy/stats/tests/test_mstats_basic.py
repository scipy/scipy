"""
Tests for the stats.mstats module (support for masked arrays)
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
    assert_array_almost_equal, assert_array_almost_equal_nulp, assert_,
    assert_allclose, assert_raises)

from nose.tools import nottest

class TestMquantiles(TestCase):
    def test_mquantiles_limit_keyword(self):
        # Regression test for Trac ticket #867
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

        desired1 = mstats.gmean(a,axis=0)
        assert_array_almost_equal(actual, desired1, decimal=14)

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
        # Tests some computations of Spearman's rho
        (x, y) = ([5.05,6.75,3.21,2.66],[1.65,2.64,2.64,6.95])
        assert_almost_equal(mstats.spearmanr(x,y)[0], -0.6324555)
        (x, y) = ([5.05,6.75,3.21,2.66,np.nan],[1.65,2.64,2.64,6.95,np.nan])
        (x, y) = (ma.fix_invalid(x), ma.fix_invalid(y))
        assert_almost_equal(mstats.spearmanr(x,y)[0], -0.6324555)

        x = [2.0, 47.4, 42.0, 10.8, 60.1, 1.7, 64.0, 63.1,
              1.0, 1.4, 7.9, 0.3, 3.9, 0.3, 6.7]
        y = [22.6, 8.3, 44.4, 11.9, 24.6, 0.6, 5.7, 41.6,
              0.0, 0.6, 6.7, 3.8, 1.0, 1.2, 1.4]
        assert_almost_equal(mstats.spearmanr(x,y)[0], 0.6887299)
        x = [2.0, 47.4, 42.0, 10.8, 60.1, 1.7, 64.0, 63.1,
              1.0, 1.4, 7.9, 0.3, 3.9, 0.3, 6.7, np.nan]
        y = [22.6, 8.3, 44.4, 11.9, 24.6, 0.6, 5.7, 41.6,
              0.0, 0.6, 6.7, 3.8, 1.0, 1.2, 1.4, np.nan]
        (x, y) = (ma.fix_invalid(x), ma.fix_invalid(y))
        assert_almost_equal(mstats.spearmanr(x,y)[0], 0.6887299)

    def test_kendalltau(self):
        # Tests some computations of Kendall's tau
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
        # Tests the seasonal Kendall tau.
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
        x = [1,0,1,1,1,1,0,1,0,0,0,1,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,
             0,0,0,0,1,-1]
        y = [14.8,13.8,12.4,10.1,7.1,6.1,5.8,4.6,4.3,3.5,3.3,3.2,3.0,
             2.8,2.8,2.5,2.4,2.3,2.1,1.7,1.7,1.5,1.3,1.3,1.2,1.2,1.1,
             0.8,0.7,0.6,0.5,0.2,0.2,0.1,np.nan]
        assert_almost_equal(mstats.pointbiserialr(x, y)[0], 0.36149, 5)


class TestTrimming(TestCase):

    def test_trim(self):
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

        a = ma.arange(12)
        a[[0,-1]] = a[5] = masked
        assert_equal(mstats.trim(a,(2,8)),
                     [None,None,2,3,4,None,6,7,8,None,None,None])

        x = ma.arange(100).reshape(10,10)
        trimx = mstats.trim(x,(0.1,0.2),relative=True,axis=None)
        assert_equal(trimx._mask.ravel(),[1]*10+[0]*70+[1]*20)
        trimx = mstats.trim(x,(0.1,0.2),relative=True,axis=0)
        assert_equal(trimx._mask.ravel(),[1]*10+[0]*70+[1]*20)
        trimx = mstats.trim(x,(0.1,0.2),relative=True,axis=-1)
        assert_equal(trimx._mask.T.ravel(),[1]*10+[0]*70+[1]*20)

        x = ma.arange(110).reshape(11,10)
        x[1] = masked
        trimx = mstats.trim(x,(0.1,0.2),relative=True,axis=None)
        assert_equal(trimx._mask.ravel(),[1]*20+[0]*70+[1]*20)
        trimx = mstats.trim(x,(0.1,0.2),relative=True,axis=0)
        assert_equal(trimx._mask.ravel(),[1]*20+[0]*70+[1]*20)
        trimx = mstats.trim(x.T,(0.1,0.2),relative=True,axis=-1)
        assert_equal(trimx.T._mask.ravel(),[1]*20+[0]*70+[1]*20)

    def test_trim_old(self):
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
        data = ma.array([77, 87, 88,114,151,210,219,246,253,262,
                         296,299,306,376,428,515,666,1310,2611])
        assert_almost_equal(mstats.trimmed_mean(data,0.1), 343, 0)
        assert_almost_equal(mstats.trimmed_mean(data,(0.1,0.1)), 343, 0)
        assert_almost_equal(mstats.trimmed_mean(data,(0.2,0.2)), 283, 0)

    def test_trimmed_stde(self):
        data = ma.array([77, 87, 88,114,151,210,219,246,253,262,
                         296,299,306,376,428,515,666,1310,2611])
        assert_almost_equal(mstats.trimmed_stde(data,(0.2,0.2)), 56.13193, 5)
        assert_almost_equal(mstats.trimmed_stde(data,0.2), 56.13193, 5)

    def test_winsorization(self):
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
        y = mstats.moment(self.testcase,1)
        assert_almost_equal(y,0.0,10)
        y = mstats.moment(self.testcase,2)
        assert_almost_equal(y,1.25)
        y = mstats.moment(self.testcase,3)
        assert_almost_equal(y,0.0)
        y = mstats.moment(self.testcase,4)
        assert_almost_equal(y,2.5625)

    def test_variation(self):
        y = mstats.variation(self.testcase)
        assert_almost_equal(y,0.44721359549996, 10)

    def test_skewness(self):
        y = mstats.skew(self.testmathworks)
        assert_almost_equal(y,-0.29322304336607,10)
        y = mstats.skew(self.testmathworks,bias=0)
        assert_almost_equal(y,-0.437111105023940,10)
        y = mstats.skew(self.testcase)
        assert_almost_equal(y,0.0,10)

    def test_kurtosis(self):
        # Set flags for axis = 0 and fisher=0 (Pearson's definition of kurtosis
        # for compatibility with Matlab)
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
        # This is not in R, so used:
        #     mean(testcase, axis=0) / (sqrt(var(testcase)*3/4))
        y = mstats.signaltonoise(self.testcase)
        assert_almost_equal(y,2.236067977)

    def test_sem(self):
        # This is not in R, so used: sqrt(var(testcase)*3/4) / sqrt(3)
        y = mstats.sem(self.testcase)
        assert_almost_equal(y,0.6454972244)

    def test_zmap(self):
        # This is not in R, so tested by using:
        #    (testcase[i]-mean(testcase,axis=0)) / sqrt(var(testcase)*3/4)
        y = mstats.zmap(self.testcase, self.testcase)
        desired_unmaskedvals = ([-1.3416407864999, -0.44721359549996,
                                 0.44721359549996, 1.3416407864999])
        assert_array_almost_equal(desired_unmaskedvals,
                                  y.data[y.mask == False], decimal=12)

    def test_zscore(self):
        # This is not in R, so tested by using:
        #     (testcase[i]-mean(testcase,axis=0)) / sqrt(var(testcase)*3/4)
        y = mstats.zscore(self.testcase)
        desired = ma.fix_invalid([-1.3416407864999, -0.44721359549996,
                                  0.44721359549996, 1.3416407864999, np.nan])
        assert_almost_equal(desired, y, decimal=12)


class TestMisc(TestCase):

    def test_obrientransform(self):
        args = [[5]*5+[6]*11+[7]*9+[8]*3+[9]*2+[10]*2,
                [6]+[7]*2+[8]*4+[9]*9+[10]*16]
        result = [5*[3.1828]+11*[0.5591]+9*[0.0344]+3*[1.6086]+2*[5.2817]+2*[11.0538],
                  [10.4352]+2*[4.8599]+4*[1.3836]+9*[0.0061]+16*[0.7277]]
        assert_almost_equal(np.round(mstats.obrientransform(*args).T,4),
                            result,4)

    def test_kstwosamp(self):
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
    # Regress a line with sinusoidal noise. Test for #1273.
    x = np.linspace(0, 100, 100)
    y = 0.2 * np.linspace(0, 100, 100) + 10
    y += np.sin(np.linspace(0, 20, 100))

    slope, intercept, r_value, p_value, sterr = mstats.linregress(x, y)
    assert_almost_equal(slope, 0.19644990055858422)
    assert_almost_equal(intercept, 10.211269918932341)


def test_plotting_positions():
    # Regression test for #1256
    pos = mstats.plotting_positions(np.arange(3), 0, 0)
    assert_array_almost_equal(pos.data, np.array([0.25, 0.5, 0.75]))


class TestNormalitytests():

    def test_vs_nonmasked(self):
        x = np.array((-2,-1,0,1,2,3)*4)**2
        assert_array_almost_equal(mstats.normaltest(x), stats.normaltest(x))
        assert_array_almost_equal(mstats.skewtest(x), stats.skewtest(x))
        assert_array_almost_equal(mstats.kurtosistest(x),
                                  stats.kurtosistest(x))

        funcs = [stats.normaltest, stats.skewtest, stats.kurtosistest]
        mfuncs = [mstats.normaltest, mstats.skewtest, mstats.kurtosistest]
        x = [1, 2, 3, 4]
        for func, mfunc in zip(funcs, mfuncs):
            assert_raises(ValueError, func, x)
            assert_raises(ValueError, mfunc, x)

    def test_axis_None(self):
        # Test axis=None (equal to axis=0 for 1-D input)
        x = np.array((-2,-1,0,1,2,3)*4)**2
        assert_allclose(mstats.normaltest(x, axis=None), mstats.normaltest(x))
        assert_allclose(mstats.skewtest(x, axis=None), mstats.skewtest(x))
        assert_allclose(mstats.kurtosistest(x, axis=None),
                        mstats.kurtosistest(x))

    def test_maskedarray_input(self):
        # Add some masked values, test result doesn't change
        x = np.array((-2,-1,0,1,2,3)*4)**2
        xm = np.ma.array(np.r_[np.inf, x, 10],
                         mask=np.r_[True, [False] * x.size, True])
        assert_allclose(mstats.normaltest(xm), stats.normaltest(x))
        assert_allclose(mstats.skewtest(xm), stats.skewtest(x))
        assert_allclose(mstats.kurtosistest(xm), stats.kurtosistest(x))

    def test_nd_input(self):
        x = np.array((-2,-1,0,1,2,3)*4)**2
        x_2d = np.vstack([x] * 2).T
        for func in [mstats.normaltest, mstats.skewtest, mstats.kurtosistest]:
            res_1d = func(x)
            res_2d = func(x_2d)
            assert_allclose(res_2d[0], [res_1d[0]] * 2)
            assert_allclose(res_2d[1], [res_1d[1]] * 2)


#TODO: for all ttest functions, add tests with masked array inputs
class TestTtest_rel():

    def test_vs_nonmasked(self):
        np.random.seed(1234567)
        outcome = np.random.randn(20, 4) + [0, 0, 1, 2]

        # 1-D inputs
        res1 = stats.ttest_rel(outcome[:, 0], outcome[:, 1])
        res2 = mstats.ttest_rel(outcome[:, 0], outcome[:, 1])
        assert_allclose(res1, res2)

        # 2-D inputs
        res1 = stats.ttest_rel(outcome[:, 0], outcome[:, 1], axis=None)
        res2 = mstats.ttest_rel(outcome[:, 0], outcome[:, 1], axis=None)
        assert_allclose(res1, res2)
        res1 = stats.ttest_rel(outcome[:, :2], outcome[:, 2:], axis=0)
        res2 = mstats.ttest_rel(outcome[:, :2], outcome[:, 2:], axis=0)
        assert_allclose(res1, res2)

        # Check default is axis=0
        res3 = mstats.ttest_rel(outcome[:, :2], outcome[:, 2:])
        assert_allclose(res2, res3)

    def test_invalid_input_size(self):
        assert_raises(ValueError, mstats.ttest_rel,
                      np.arange(10), np.arange(11))
        x = np.arange(24)
        assert_raises(ValueError, mstats.ttest_rel,
                      x.reshape(2, 3, 4), x.reshape(2, 4, 3), axis=1)
        assert_raises(ValueError, mstats.ttest_rel,
                      x.reshape(2, 3, 4), x.reshape(2, 4, 3), axis=2)

    def test_empty(self):
        res1 = mstats.ttest_rel([], [])
        assert_(np.all(np.isnan(res1)))


class TestTtest_ind():

    def test_vs_nonmasked(self):
        np.random.seed(1234567)
        outcome = np.random.randn(20, 4) + [0, 0, 1, 2]

        # 1-D inputs
        res1 = stats.ttest_ind(outcome[:, 0], outcome[:, 1])
        res2 = mstats.ttest_ind(outcome[:, 0], outcome[:, 1])
        assert_allclose(res1, res2)

        # 2-D inputs
        res1 = stats.ttest_ind(outcome[:, 0], outcome[:, 1], axis=None)
        res2 = mstats.ttest_ind(outcome[:, 0], outcome[:, 1], axis=None)
        assert_allclose(res1, res2)
        res1 = stats.ttest_ind(outcome[:, :2], outcome[:, 2:], axis=0)
        res2 = mstats.ttest_ind(outcome[:, :2], outcome[:, 2:], axis=0)
        assert_allclose(res1, res2)

        # Check default is axis=0
        res3 = mstats.ttest_ind(outcome[:, :2], outcome[:, 2:])
        assert_allclose(res2, res3)

    def test_empty(self):
        res1 = mstats.ttest_ind([], [])
        assert_(np.all(np.isnan(res1)))


class TestTtest_1samp():

    def test_vs_nonmasked(self):
        np.random.seed(1234567)
        outcome = np.random.randn(20, 4) + [0, 0, 1, 2]

        # 1-D inputs
        res1 = stats.ttest_1samp(outcome[:, 0], 1)
        res2 = mstats.ttest_1samp(outcome[:, 0], 1)
        assert_allclose(res1, res2)

        # 2-D inputs
        res1 = stats.ttest_1samp(outcome[:, 0], outcome[:, 1], axis=None)
        res2 = mstats.ttest_1samp(outcome[:, 0], outcome[:, 1], axis=None)
        assert_allclose(res1, res2)
        res1 = stats.ttest_1samp(outcome[:, :2], outcome[:, 2:], axis=0)
        res2 = mstats.ttest_1samp(outcome[:, :2], outcome[:, 2:], axis=0)
        assert_allclose(res1, res2)

        # Check default is axis=0
        res3 = mstats.ttest_1samp(outcome[:, :2], outcome[:, 2:])
        assert_allclose(res2, res3)

    def test_empty(self):
        res1 = mstats.ttest_1samp([], 1)
        assert_(np.all(np.isnan(res1)))




class TestCompareWithStats(TestCase):
    """
    Class to compare mstats results with stats results
    It is in general assumed that scipy.stats is at a more mature stage than stats.mstats. If a routine in mstats
    results in similar results like in scipy.stats, this is considered also as a proper validation of scipy.mstats routine

    Different sample sizes are used for testing, as some problems between stats and mstats are dependent on sample size

    Author: Alexander Loew

    NOTE that some tests fail. This might be caused by
    a) actual differences or bugs between stats and mstats
    b) numerical inaccuracies
    c) different definitions of routine interfaces

    These failures need to be checked. Current workaround is to have disabled these tests,
    but issuing reports on scipy-dev

    """


    def get_n(self):
        """returns list of sample sizes to be used for comparison"""
        return [1000,100,10,5]

    def generate_xy_sample(self,n):
        """
        generate some sample data
        This routine generates numpy arrays and corresponding masked arrays with the same data, but additional
        masked values
        """

        assert(isinstance(n,int))
        assert(n>3)

        x = np.random.randn(n); y = x + np.random.randn(n) #normal numpy arrays
        xm = np.ones(len(x)+5)*np.nan; ym = np.ones(len(y)+5)*np.nan #add here a few artificial samples that will be masked
        xm[0:len(x)] = x; ym[0:len(y)] = y
        xm = np.ma.array(xm,mask=np.isnan(xm))
        ym = np.ma.array(ym,mask=np.isnan(ym))

        return x,y,xm,ym #x,y are numpy arrays, while xm,ym are masked arrays

    def test_linregress(self):
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)
            slope, intercept, r_value, p_value, std_err      = stats.linregress(x,y)
            slopem, interceptm, r_valuem, p_valuem, std_errm = stats.mstats.linregress(xm,ym)

            assert_almost_equal(slope,slopem,10)
            assert_almost_equal(intercept,interceptm,10)
            assert_almost_equal(r_value,r_valuem,10)
            #~ assert_almost_equal(p_value,p_valuem,10) #ERROR invalid p-value and std_err todo
            #~ assert_almost_equal(std_err,std_errm,10) #todo

    def test_pearsonr(self):
        """ test for pearsonr """
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)
            r,p   = stats.pearsonr(x,y)
            rm,pm = stats.mstats.pearsonr(xm,ym)

            assert_almost_equal(r,rm,10)
            assert_almost_equal(p,pm,10)

    def test_spearmanr(self):
        """ test spearmanr """
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)
            r,p = stats.spearmanr(x,y)
            rm,pm = stats.mstats.spearmanr(xm,ym)
            assert_almost_equal(r,rm,10)
            assert_almost_equal(p,pm,10)

    def test_gmean(self):
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)

            r  = stats.gmean(abs(x))
            rm = stats.mstats.gmean(abs(xm))
            assert_almost_equal(r,rm,10)

            r  = stats.gmean(abs(y))
            rm = stats.mstats.gmean(abs(ym))
            assert_almost_equal(r,rm,10)

    def test_hmean(self):
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)

            r  = stats.hmean(abs(x))
            rm = stats.mstats.hmean(abs(xm))
            assert_almost_equal(r,rm,10)

            r  = stats.hmean(abs(y))
            rm = stats.mstats.hmean(abs(ym))
            assert_almost_equal(r,rm,10)

    def test_skew(self):
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)

            r = stats.skew(x)
            rm = stats.mstats.skew(xm)
            assert_almost_equal(r,rm,10)

            r = stats.skew(y)
            rm = stats.mstats.skew(ym)
            assert_almost_equal(r,rm,10)

    def test_moment(self):
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)

            r = stats.moment(x)
            rm = stats.mstats.moment(xm)
            assert_almost_equal(r,rm,10)

            r = stats.moment(y)
            rm = stats.mstats.moment(ym)
            assert_almost_equal(r,rm,10)

    def test_signaltonoise(self):
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)

            r = stats.signaltonoise(x)
            rm = stats.mstats.signaltonoise(xm)
            assert_almost_equal(r,rm,10)

            r = stats.signaltonoise(y)
            rm = stats.mstats.signaltonoise(ym)
            assert_almost_equal(r,rm,10)

    def test_betai(self):
        """ test incomplete beta function """
        for i in range(10):
            a = np.random.rand()*5.; b=np.random.rand()*200.

            assert_equal(stats.betai(a,b,0.),0.) #for x=0 result should be always 0.
            assert_equal(stats.betai(a,b,1.),1.) #for x=1 result should be always 1.
            assert_equal(stats.mstats.betai(a,b,0.),0.) #for x=0 result should be always 0.
            assert_equal(stats.mstats.betai(a,b,1.),1.) #for x=1 result should be always 1.

        #now some further random samples
        for i in range(10):
            a = np.random.rand()*5.; b=np.random.rand()*200.; x=np.random.rand()
            assert_equal(stats.betai(a,b,x),stats.mstats.betai(a,b,x))


    def test_zscore(self):
        """ test zscore """
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)

            #reference solution
            zx = (x-x.mean()) / x.std()
            zy = (y-y.mean()) / y.std()

            #validate stats
            assert(np.any(abs(stats.zscore(x)-zx) < 1.E-10))
            assert(np.any(abs(stats.zscore(y)-zy) < 1.E-10))

            #compare stats and mstats
            assert(np.any(abs(stats.zscore(x)-stats.mstats.zscore(xm[0:len(x)])) < 1.E-10))
            assert(np.any(abs(stats.zscore(y)-stats.mstats.zscore(ym[0:len(y)])) < 1.E-10))


    def test_kurtosis(self):
        """ test kurtosis """
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)

            r  = stats.kurtosis(x)
            rm = stats.mstats.kurtosis(xm)
            assert_almost_equal(r,rm,10)

            r  = stats.kurtosis(y)
            rm = stats.mstats.kurtosis(ym)
            assert_almost_equal(r,rm,10)

    @nottest
    def test_sem(self):
        #example from stats.sem doc
        a = np.arange(20).reshape(5,4)
        am = np.ma.array(a)
        r = stats.sem(a)
        rm = stats.mstats.sem(am)

        assert(np.all(abs(r - 2.82842712)<1.E-5 ))
        assert(np.all(abs(rm - 2.82842712)<1.E-5))

        #Find standard error across the whole array, using n degrees of freedom
        #~ assert_almost_equal(stats.mstats.sem(am,axis=None,ddof=0),stats.sem(a, axis=None, ddof=0),10)   todo #ERROR: mstats contains no ddof parameter

        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n) #ddof default is 0
            assert_almost_equal(stats.mstats.sem(xm,axis=None),stats.sem(x, axis=None, ddof=0),10) # todo ERROR: results are different at 4-5 decimal
            assert_almost_equal(stats.mstats.sem(ym,axis=None),stats.sem(y, axis=None, ddof=0),10) #todo

    def test_describe(self):
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)
            r  = stats.describe(x)
            rm = stats.mstats.describe(xm)

            assert_almost_equal(r[0],rm[0],10) #n
            assert_almost_equal(r[1][0],rm[1][0],10) #min
            assert_almost_equal(r[1][1],rm[1][1],10) #max
            assert_almost_equal(r[2],rm[2],10) #mean
            #assert_almost_equal(r[3],rm[3],10) #unbiased variance ERROR: throws an assertion error! todo
            assert_almost_equal(r[4],rm[4],10) #biased skewness
            assert_almost_equal(r[5],rm[5],10) #biased kurtosis

    def test_rankdata(self):
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)
            r  = stats.rankdata(x)
            rm = stats.mstats.rankdata(x)

            assert(np.all( (r-rm) == 0. ))

    def test_tmean(self):
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)

            assert_almost_equal(stats.tmean(x),stats.mstats.tmean(xm),10)
            assert_almost_equal(stats.tmean(y),stats.mstats.tmean(ym),10)

            #assert_almost_equal(stats.tmean(x,limits=(1.,3.)),stats.mstats.tmean(xm,limits=(1.,0.7)),10)
            #assert_almost_equal(stats.tmean(y,limits=(0.5,0.7)),stats.mstats.tmean(ym,limits=(0.5,0.7)),10)


    def test_tmax(self):
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)
            assert_almost_equal(stats.tmax(x,2.),stats.mstats.tmax(xm,2.),10)
            assert_almost_equal(stats.tmax(y,2.),stats.mstats.tmax(ym,2.),10)

    def test_tmin(self):
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)
            #~ todo assert_almost_equal(stats.tmin(x),stats.mstats.tmin(xm),10) #ERROR: causes trouble without keyword lowerlimit in mstats
            #~ assert_almost_equal(stats.tmin(y),stats.mstats.tmin(ym),10) #todo

            assert_almost_equal(stats.tmin(x,lowerlimit=-1.),stats.mstats.tmin(xm,lowerlimit=-1.),10)
            assert_almost_equal(stats.tmin(y,lowerlimit=-1.),stats.mstats.tmin(ym,lowerlimit=-1.),10)

    def test_zmap(self):
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)
            z  = stats.zmap(x,y)
            zm = stats.mstats.zmap(xm,ym)
            assert(np.all( abs(z - zm[0:len(z)]) < 1.E-10 ))

    def test_variation(self):
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)
            assert_almost_equal(stats.variation(x),stats.mstats.variation(xm),10)
            assert_almost_equal(stats.variation(y),stats.mstats.variation(ym),10)

    @nottest
    def test_tvar(self):
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)
            assert_almost_equal(stats.tvar(x),stats.mstats.tvar(xm),10) #ERROR: throws an assertion error
            assert_almost_equal(stats.tvar(y),stats.mstats.tvar(ym),10)

    def test_trimboth(self):
        a = np.arange(20)
        b = stats.trimboth(a,0.1)
        bm = stats.mstats.trimboth(a,0.1)

        assert(np.all(b == bm.data[~bm.mask]))

    @nottest
    def test_tsem(self):
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)
            assert_almost_equal(stats.tsem(x),stats.mstats.tsem(xm),10)
            assert_almost_equal(stats.tsem(y),stats.mstats.tsem(ym),10) #error
            assert_almost_equal(stats.tsem(x,limits=(-2.,2.)),stats.mstats.tsem(xm,limits=(-2.,2.)),10) #ERROR: causes problems with limits!!! (at 4th digit)

    @nottest
    def test_skewtest(self):
        for n in self.get_n():
            if n>8:
                x,y,xm,ym = self.generate_xy_sample(n)
                r = stats.skewtest(x)
                rm = stats.mstats.skewtest(xm)
                assert_almost_equal(r[0],rm[0],10)
                assert_almost_equal(r[1],rm[1],10) #<<< seems to be an inconsistency between modules!

    def test_normaltest(self):
        for n in self.get_n():
            if n > 8:
                x,y,xm,ym = self.generate_xy_sample(n)
                r = stats.normaltest(x)
                rm = stats.mstats.normaltest(xm)
                assert_almost_equal(r[0],rm[0],10)
                assert_almost_equal(r[1],rm[1],10)

    def test_find_repeats(self):
        x = np.asarray([1,1,2,2,3,3,3,4,4,4,4]).astype('float')
        tmp = np.asarray([1,1,2,2,3,3,3,4,4,4,4,5,5,5,5]).astype('float')
        xm = np.ma.array(tmp,mask=tmp==5.)

        r  = stats.find_repeats(x)
        rm = stats.mstats.find_repeats(xm)

        assert_equal(r,rm)

    def test_kendalltau(self):
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)
            r  = stats.kendalltau(x,y)
            rm = stats.mstats.kendalltau(xm,ym)
            assert_almost_equal(r[0],rm[0],10)
            assert_almost_equal(r[1],rm[1],7)

    @nottest
    def test_obrientransform(self):    #todo causes error!
        for n in self.get_n():
            x,y,xm,ym = self.generate_xy_sample(n)
            r = stats.obrientransform(x)
            rm = stats.mstats.obrientransform(xm)
            assert_almost_equal(r,rm[0:len(x)]) #ERROR: returned array is transposed in mstats compared to stats


if __name__ == "__main__":
    run_module_suite()
