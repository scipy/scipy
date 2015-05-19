""" Test functions for stats module

    WRITTEN BY LOUIS LUANGKESORN <lluang@yahoo.com> FOR THE STATS MODULE
    BASED ON WILKINSON'S STATISTICS QUIZ
    http://www.stanford.edu/~clint/bench/wilk.txt

    Additional tests by a host of SciPy developers.
"""
from __future__ import division, print_function, absolute_import

import sys
import warnings
from collections import namedtuple

from numpy.testing import (TestCase, assert_, assert_equal,
                           assert_almost_equal, assert_array_almost_equal,
                           assert_array_equal, assert_approx_equal,
                           assert_raises, run_module_suite, assert_allclose,
                           dec)
import numpy.ma.testutils as mat
from numpy import array, arange, float32, float64, power
import numpy as np

import scipy.stats as stats
from common_tests import check_named_results

""" Numbers in docstrings beginning with 'W' refer to the section numbers
    and headings found in the STATISTICS QUIZ of Leland Wilkinson.  These are
    considered to be essential functionality.  True testing and
    evaluation of a statistics package requires use of the
    NIST Statistical test data.  See McCoullough(1999) Assessing The Reliability
    of Statistical Software for a test methodology and its
    implementation in testing SAS, SPSS, and S-Plus
"""

#  Datasets
#  These data sets are from the nasty.dat sets used by Wilkinson
#  For completeness, I should write the relevant tests and count them as failures
#  Somewhat acceptable, since this is still beta software.  It would count as a
#  good target for 1.0 status
X = array([1,2,3,4,5,6,7,8,9], float)
ZERO = array([0,0,0,0,0,0,0,0,0], float)
BIG = array([99999991,99999992,99999993,99999994,99999995,99999996,99999997,
             99999998,99999999], float)
LITTLE = array([0.99999991,0.99999992,0.99999993,0.99999994,0.99999995,0.99999996,
                0.99999997,0.99999998,0.99999999], float)
HUGE = array([1e+12,2e+12,3e+12,4e+12,5e+12,6e+12,7e+12,8e+12,9e+12], float)
TINY = array([1e-12,2e-12,3e-12,4e-12,5e-12,6e-12,7e-12,8e-12,9e-12], float)
ROUND = array([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5], float)


class TestTrimmedStats(TestCase):
    # TODO: write these tests to handle missing values properly
    dprec = np.finfo(np.float64).precision

    def test_tmean(self):
        y = stats.tmean(X, (2, 8), (True, True))
        assert_approx_equal(y, 5.0, significant=self.dprec)

        y1 = stats.tmean(X, limits=(2, 8), inclusive=(False, False))
        y2 = stats.tmean(X, limits=None)
        assert_approx_equal(y1, y2, significant=self.dprec)

    def test_tvar(self):
        y = stats.tvar(X, limits=(2, 8), inclusive=(True, True))
        assert_approx_equal(y, 4.6666666666666661, significant=self.dprec)

        y = stats.tvar(X, limits=None)
        assert_approx_equal(y, X.var(ddof=1), significant=self.dprec)

    def test_tstd(self):
        y = stats.tstd(X, (2, 8), (True, True))
        assert_approx_equal(y, 2.1602468994692865, significant=self.dprec)

        y = stats.tstd(X, limits=None)
        assert_approx_equal(y, X.std(ddof=1), significant=self.dprec)

    def test_tmin(self):
        assert_equal(stats.tmin(4), 4)

        x = np.arange(10)
        assert_equal(stats.tmin(x), 0)
        assert_equal(stats.tmin(x, lowerlimit=0), 0)
        assert_equal(stats.tmin(x, lowerlimit=0, inclusive=False), 1)

        x = x.reshape((5, 2))
        assert_equal(stats.tmin(x, lowerlimit=0, inclusive=False), [2, 1])
        assert_equal(stats.tmin(x, axis=1), [0, 2, 4, 6, 8])
        assert_equal(stats.tmin(x, axis=None), 0)

    def test_tmax(self):
        assert_equal(stats.tmax(4), 4)

        x = np.arange(10)
        assert_equal(stats.tmax(x), 9)
        assert_equal(stats.tmax(x, upperlimit=9), 9)
        assert_equal(stats.tmax(x, upperlimit=9, inclusive=False), 8)

        x = x.reshape((5, 2))
        assert_equal(stats.tmax(x, upperlimit=9, inclusive=False), [8, 7])
        assert_equal(stats.tmax(x, axis=1), [1, 3, 5, 7, 9])
        assert_equal(stats.tmax(x, axis=None), 9)

    def test_tsem(self):
        y = stats.tsem(X, limits=(3, 8), inclusive=(False, True))
        y_ref = np.array([4, 5, 6, 7, 8])
        assert_approx_equal(y, y_ref.std(ddof=1) / np.sqrt(y_ref.size),
                            significant=self.dprec)

        assert_approx_equal(stats.tsem(X, limits=[-1, 10]),
                            stats.tsem(X, limits=None),
                            significant=self.dprec)


class TestNanFunc(TestCase):
    def __init__(self, *args, **kw):
        TestCase.__init__(self, *args, **kw)
        self.X = X.copy()

        self.Xall = X.copy()
        self.Xall[:] = np.nan

        self.Xsome = X.copy()
        self.Xsomet = X.copy()
        self.Xsome[0] = np.nan
        self.Xsomet = self.Xsomet[1:]

    def test_nanmean_none(self):
        # Check nanmean when no values are nan.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', DeprecationWarning)
            m = stats.nanmean(X)
            assert_approx_equal(m, X[4])

    def test_nanmean_some(self):
        # Check nanmean when some values only are nan.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', DeprecationWarning)
            m = stats.nanmean(self.Xsome)
            assert_approx_equal(m, 5.5)

    def test_nanmean_all(self):
        # Check nanmean when all values are nan.
        with warnings.catch_warnings():
            warns = (DeprecationWarning, RuntimeWarning)
            warnings.simplefilter('ignore', warns)
            with np.errstate(invalid='ignore'):
                m = stats.nanmean(self.Xall)
        assert_(np.isnan(m))

    def test_nanstd_none(self):
        # Check nanstd when no values are nan.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', DeprecationWarning)
            s = stats.nanstd(self.X)
        assert_approx_equal(s, np.std(self.X, ddof=1))

    def test_nanstd_some(self):
        # Check nanstd when some values only are nan.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', DeprecationWarning)
            s = stats.nanstd(self.Xsome)
        assert_approx_equal(s, np.std(self.Xsomet, ddof=1))

    def test_nanstd_all(self):
        # Check nanstd when all values are nan.
        with warnings.catch_warnings():
            warns = (DeprecationWarning, RuntimeWarning)
            warnings.simplefilter('ignore', warns)
            with np.errstate(invalid='ignore'):
                s = stats.nanstd(self.Xall)
        assert_(np.isnan(s))

    def test_nanstd_bias_kw(self):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            s = stats.nanstd(self.X, bias=True)
        assert_approx_equal(s, np.std(self.X, ddof=0))

    def test_nanstd_negative_axis(self):
        x = np.array([1, 2, 3])
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            res = stats.nanstd(x, -1)
        assert_equal(res, 1)

    def test_nanmedian_none(self):
        # Check nanmedian when no values are nan.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', DeprecationWarning)
            m = stats.nanmedian(self.X)
        assert_approx_equal(m, np.median(self.X))

    def test_nanmedian_axis(self):
        # Check nanmedian with axis
        X = self.X.reshape(3,3)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', DeprecationWarning)
            m = stats.nanmedian(X, axis=0)
            assert_equal(m, np.median(X, axis=0))
            m = stats.nanmedian(X, axis=1)
            assert_equal(m, np.median(X, axis=1))

    def test_nanmedian_some(self):
        # Check nanmedian when some values only are nan.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', DeprecationWarning)
            m = stats.nanmedian(self.Xsome)
        assert_approx_equal(m, np.median(self.Xsomet))

    def test_nanmedian_all(self):
        # Check nanmedian when all values are nan.
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            m = stats.nanmedian(self.Xall)
            assert_(np.isnan(m))
            assert_equal(len(w), 2)  # Deprecation & RuntimeWarning
            assert_(issubclass(w[1].category, RuntimeWarning))

    def test_nanmedian_all_axis(self):
        # Check nanmedian when all values are nan.
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            m = stats.nanmedian(self.Xall.reshape(3,3), axis=1)
            assert_(np.isnan(m).all())
            assert_equal(len(w), 4)
            assert_(issubclass(w[-1].category, RuntimeWarning))

    def test_nanmedian_scalars(self):
        # Check nanmedian for scalar inputs. See ticket #1098.
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=RuntimeWarning)
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            assert_equal(stats.nanmedian(1), np.median(1))
            assert_equal(stats.nanmedian(True), np.median(True))
            assert_equal(stats.nanmedian(np.array(1)), np.median(np.array(1)))
            assert_equal(stats.nanmedian(np.nan), np.median(np.nan))


class TestCorrPearsonr(TestCase):
    """ W.II.D. Compute a correlation matrix on all the variables.

        All the correlations, except for ZERO and MISS, shoud be exactly 1.
        ZERO and MISS should have undefined or missing correlations with the
        other variables.  The same should go for SPEARMAN corelations, if
        your program has them.
    """
    def test_pXX(self):
        y = stats.pearsonr(X,X)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pXBIG(self):
        y = stats.pearsonr(X,BIG)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pXLITTLE(self):
        y = stats.pearsonr(X,LITTLE)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pXHUGE(self):
        y = stats.pearsonr(X,HUGE)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pXTINY(self):
        y = stats.pearsonr(X,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pXROUND(self):
        y = stats.pearsonr(X,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pBIGBIG(self):
        y = stats.pearsonr(BIG,BIG)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pBIGLITTLE(self):
        y = stats.pearsonr(BIG,LITTLE)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pBIGHUGE(self):
        y = stats.pearsonr(BIG,HUGE)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pBIGTINY(self):
        y = stats.pearsonr(BIG,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pBIGROUND(self):
        y = stats.pearsonr(BIG,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pLITTLELITTLE(self):
        y = stats.pearsonr(LITTLE,LITTLE)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pLITTLEHUGE(self):
        y = stats.pearsonr(LITTLE,HUGE)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pLITTLETINY(self):
        y = stats.pearsonr(LITTLE,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pLITTLEROUND(self):
        y = stats.pearsonr(LITTLE,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pHUGEHUGE(self):
        y = stats.pearsonr(HUGE,HUGE)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pHUGETINY(self):
        y = stats.pearsonr(HUGE,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pHUGEROUND(self):
        y = stats.pearsonr(HUGE,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pTINYTINY(self):
        y = stats.pearsonr(TINY,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pTINYROUND(self):
        y = stats.pearsonr(TINY,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_pROUNDROUND(self):
        y = stats.pearsonr(ROUND,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_r_exactly_pos1(self):
        a = arange(3.0)
        b = a
        r, prob = stats.pearsonr(a,b)
        assert_equal(r, 1.0)
        assert_equal(prob, 0.0)

    def test_r_exactly_neg1(self):
        a = arange(3.0)
        b = -a
        r, prob = stats.pearsonr(a,b)
        assert_equal(r, -1.0)
        assert_equal(prob, 0.0)

    def test_basic(self):
        # A basic test, with a correlation coefficient
        # that is not 1 or -1.
        a = array([-1, 0, 1])
        b = array([0, 0, 3])
        r, prob = stats.pearsonr(a, b)
        assert_approx_equal(r, np.sqrt(3)/2)
        assert_approx_equal(prob, 1.0/3)


class TestFisherExact(TestCase):
    """Some tests to show that fisher_exact() works correctly.

    Note that in SciPy 0.9.0 this was not working well for large numbers due to
    inaccuracy of the hypergeom distribution (see #1218). Fixed now.

    Also note that R and Scipy have different argument formats for their
    hypergeometric distribution functions.

    R:
    > phyper(18999, 99000, 110000, 39000, lower.tail = FALSE)
    [1] 1.701815e-09
    """
    def test_basic(self):
        fisher_exact = stats.fisher_exact

        res = fisher_exact([[14500, 20000], [30000, 40000]])[1]
        assert_approx_equal(res, 0.01106, significant=4)
        res = fisher_exact([[100, 2], [1000, 5]])[1]
        assert_approx_equal(res, 0.1301, significant=4)
        res = fisher_exact([[2, 7], [8, 2]])[1]
        assert_approx_equal(res, 0.0230141, significant=6)
        res = fisher_exact([[5, 1], [10, 10]])[1]
        assert_approx_equal(res, 0.1973244, significant=6)
        res = fisher_exact([[5, 15], [20, 20]])[1]
        assert_approx_equal(res, 0.0958044, significant=6)
        res = fisher_exact([[5, 16], [20, 25]])[1]
        assert_approx_equal(res, 0.1725862, significant=6)
        res = fisher_exact([[10, 5], [10, 1]])[1]
        assert_approx_equal(res, 0.1973244, significant=6)
        res = fisher_exact([[5, 0], [1, 4]])[1]
        assert_approx_equal(res, 0.04761904, significant=6)
        res = fisher_exact([[0, 1], [3, 2]])[1]
        assert_approx_equal(res, 1.0)
        res = fisher_exact([[0, 2], [6, 4]])[1]
        assert_approx_equal(res, 0.4545454545)
        res = fisher_exact([[2, 7], [8, 2]])
        assert_approx_equal(res[1], 0.0230141, significant=6)
        assert_approx_equal(res[0], 4.0 / 56)

    def test_precise(self):
        # results from R
        #
        # R defines oddsratio differently (see Notes section of fisher_exact
        # docstring), so those will not match.  We leave them in anyway, in
        # case they will be useful later on. We test only the p-value.
        tablist = [
            ([[100, 2], [1000, 5]], (2.505583993422285e-001, 1.300759363430016e-001)),
            ([[2, 7], [8, 2]], (8.586235135736206e-002, 2.301413756522114e-002)),
            ([[5, 1], [10, 10]], (4.725646047336584e+000, 1.973244147157190e-001)),
            ([[5, 15], [20, 20]], (3.394396617440852e-001, 9.580440012477637e-002)),
            ([[5, 16], [20, 25]], (3.960558326183334e-001, 1.725864953812994e-001)),
            ([[10, 5], [10, 1]], (2.116112781158483e-001, 1.973244147157190e-001)),
            ([[10, 5], [10, 0]], (0.000000000000000e+000, 6.126482213438734e-002)),
            ([[5, 0], [1, 4]], (np.inf, 4.761904761904762e-002)),
            ([[0, 5], [1, 4]], (0.000000000000000e+000, 1.000000000000000e+000)),
            ([[5, 1], [0, 4]], (np.inf, 4.761904761904758e-002)),
            ([[0, 1], [3, 2]], (0.000000000000000e+000, 1.000000000000000e+000))
            ]
        for table, res_r in tablist:
            res = stats.fisher_exact(np.asarray(table))
            np.testing.assert_almost_equal(res[1], res_r[1], decimal=11,
                                           verbose=True)

    @dec.slow
    def test_large_numbers(self):
        # Test with some large numbers. Regression test for #1401
        pvals = [5.56e-11, 2.666e-11, 1.363e-11]  # from R
        for pval, num in zip(pvals, [75, 76, 77]):
            res = stats.fisher_exact([[17704, 496], [1065, num]])[1]
            assert_approx_equal(res, pval, significant=4)

        res = stats.fisher_exact([[18000, 80000], [20000, 90000]])[1]
        assert_approx_equal(res, 0.2751, significant=4)

    def test_raises(self):
        # test we raise an error for wrong shape of input.
        assert_raises(ValueError, stats.fisher_exact,
                      np.arange(6).reshape(2, 3))

    def test_row_or_col_zero(self):
        tables = ([[0, 0], [5, 10]],
                  [[5, 10], [0, 0]],
                  [[0, 5], [0, 10]],
                  [[5, 0], [10, 0]])
        for table in tables:
            oddsratio, pval = stats.fisher_exact(table)
            assert_equal(pval, 1.0)
            assert_equal(oddsratio, np.nan)

    def test_less_greater(self):
        tables = (
            # Some tables to compare with R:
            [[2, 7], [8, 2]],
            [[200, 7], [8, 300]],
            [[28, 21], [6, 1957]],
            [[190, 800], [200, 900]],
            # Some tables with simple exact values
            # (includes regression test for ticket #1568):
            [[0, 2], [3, 0]],
            [[1, 1], [2, 1]],
            [[2, 0], [1, 2]],
            [[0, 1], [2, 3]],
            [[1, 0], [1, 4]],
            )
        pvals = (
            # from R:
            [0.018521725952066501, 0.9990149169715733],
            [1.0, 2.0056578803889148e-122],
            [1.0, 5.7284374608319831e-44],
            [0.7416227, 0.2959826],
            # Exact:
            [0.1, 1.0],
            [0.7, 0.9],
            [1.0, 0.3],
            [2./3, 1.0],
            [1.0, 1./3],
            )
        for table, pval in zip(tables, pvals):
            res = []
            res.append(stats.fisher_exact(table, alternative="less")[1])
            res.append(stats.fisher_exact(table, alternative="greater")[1])
            assert_allclose(res, pval, atol=0, rtol=1e-7)

    def test_gh3014(self):
        # check if issue #3014 has been fixed.
        # before, this would have risen a ValueError
        odds, pvalue = stats.fisher_exact([[1, 2], [9, 84419233]])


class TestCorrSpearmanr(TestCase):
    """ W.II.D. Compute a correlation matrix on all the variables.

        All the correlations, except for ZERO and MISS, shoud be exactly 1.
        ZERO and MISS should have undefined or missing correlations with the
        other variables.  The same should go for SPEARMAN corelations, if
        your program has them.
    """
    def test_scalar(self):
        y = stats.spearmanr(4., 2.)
        assert_(np.isnan(y).all())

    def test_uneven_lengths(self):
        assert_raises(ValueError, stats.spearmanr, [1, 2, 1], [8, 9])
        assert_raises(ValueError, stats.spearmanr, [1, 2, 1], 8)

    def test_sXX(self):
        y = stats.spearmanr(X,X)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sXBIG(self):
        y = stats.spearmanr(X,BIG)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sXLITTLE(self):
        y = stats.spearmanr(X,LITTLE)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sXHUGE(self):
        y = stats.spearmanr(X,HUGE)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sXTINY(self):
        y = stats.spearmanr(X,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sXROUND(self):
        y = stats.spearmanr(X,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sBIGBIG(self):
        y = stats.spearmanr(BIG,BIG)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sBIGLITTLE(self):
        y = stats.spearmanr(BIG,LITTLE)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sBIGHUGE(self):
        y = stats.spearmanr(BIG,HUGE)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sBIGTINY(self):
        y = stats.spearmanr(BIG,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sBIGROUND(self):
        y = stats.spearmanr(BIG,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sLITTLELITTLE(self):
        y = stats.spearmanr(LITTLE,LITTLE)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sLITTLEHUGE(self):
        y = stats.spearmanr(LITTLE,HUGE)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sLITTLETINY(self):
        y = stats.spearmanr(LITTLE,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sLITTLEROUND(self):
        y = stats.spearmanr(LITTLE,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sHUGEHUGE(self):
        y = stats.spearmanr(HUGE,HUGE)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sHUGETINY(self):
        y = stats.spearmanr(HUGE,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sHUGEROUND(self):
        y = stats.spearmanr(HUGE,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sTINYTINY(self):
        y = stats.spearmanr(TINY,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sTINYROUND(self):
        y = stats.spearmanr(TINY,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_sROUNDROUND(self):
        y = stats.spearmanr(ROUND,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)

    def test_spearmanr_result_attributes(self):
        res = stats.spearmanr(X, X)
        attributes = ('correlation', 'pvalue')
        check_named_results(res, attributes)


class TestCorrSpearmanrTies(TestCase):
    """Some tests of tie-handling by the spearmanr function."""

    def test_tie1(self):
        # Data
        x = [1.0, 2.0, 3.0, 4.0]
        y = [1.0, 2.0, 2.0, 3.0]
        # Ranks of the data, with tie-handling.
        xr = [1.0, 2.0, 3.0, 4.0]
        yr = [1.0, 2.5, 2.5, 4.0]
        # Result of spearmanr should be the same as applying
        # pearsonr to the ranks.
        sr = stats.spearmanr(x, y)
        pr = stats.pearsonr(xr, yr)
        assert_almost_equal(sr, pr)


#    W.II.E.  Tabulate X against X, using BIG as a case weight.  The values
#    should appear on the diagonal and the total should be 899999955.
#    If the table cannot hold these values, forget about working with
#    census data.  You can also tabulate HUGE against TINY.  There is no
#    reason a tabulation program should not be able to distinguish
#    different values regardless of their magnitude.

# I need to figure out how to do this one.


def test_kendalltau():
    # with some ties
    x1 = [12, 2, 1, 12, 2]
    x2 = [1, 4, 7, 1, 0]
    expected = (-0.47140452079103173, 0.24821309157521476)
    res = stats.kendalltau(x1, x2)
    assert_approx_equal(res[0], expected[0])
    assert_approx_equal(res[1], expected[1])

    # test for namedtuple attribute results
    attributes = ('correlation', 'pvalue')
    res = stats.kendalltau(x1, x2)
    check_named_results(res, attributes)

    # with only ties in one or both inputs
    assert_equal(stats.kendalltau([2,2,2], [2,2,2]), (np.nan, np.nan))
    assert_equal(stats.kendalltau([2,0,2], [2,2,2]), (np.nan, np.nan))
    assert_equal(stats.kendalltau([2,2,2], [2,0,2]), (np.nan, np.nan))

    # empty arrays provided as input
    assert_equal(stats.kendalltau([], []), (np.nan, np.nan))

    # check two different sort methods
    assert_approx_equal(stats.kendalltau(x1, x2, initial_lexsort=False)[1],
                        stats.kendalltau(x1, x2, initial_lexsort=True)[1])

    # and with larger arrays
    np.random.seed(7546)
    x = np.array([np.random.normal(loc=1, scale=1, size=500),
                np.random.normal(loc=1, scale=1, size=500)])
    corr = [[1.0, 0.3],
            [0.3, 1.0]]
    x = np.dot(np.linalg.cholesky(corr), x)
    expected = (0.19291382765531062, 1.1337108207276285e-10)
    res = stats.kendalltau(x[0], x[1])
    assert_approx_equal(res[0], expected[0])
    assert_approx_equal(res[1], expected[1])

    # and do we get a tau of 1 for identical inputs?
    assert_approx_equal(stats.kendalltau([1,1,2], [1,1,2])[0], 1.0)


class TestRegression(TestCase):
    def test_linregressBIGX(self):
        # W.II.F.  Regress BIG on X.
        # The constant should be 99999990 and the regression coefficient should be 1.
        y = stats.linregress(X,BIG)
        intercept = y[1]
        r = y[2]
        assert_almost_equal(intercept,99999990)
        assert_almost_equal(r,1.0)

    def test_regressXX(self):
        # W.IV.B.  Regress X on X.
        # The constant should be exactly 0 and the regression coefficient should be 1.
        # This is a perfectly valid regression.  The program should not complain.
        y = stats.linregress(X,X)
        intercept = y[1]
        r = y[2]
        assert_almost_equal(intercept,0.0)
        assert_almost_equal(r,1.0)
#     W.IV.C. Regress X on BIG and LITTLE (two predictors).  The program
#     should tell you that this model is "singular" because BIG and
#     LITTLE are linear combinations of each other.  Cryptic error
#     messages are unacceptable here.  Singularity is the most
#     fundamental regression error.
# Need to figure out how to handle multiple linear regression.  Not obvious

    def test_regressZEROX(self):
        # W.IV.D. Regress ZERO on X.
        # The program should inform you that ZERO has no variance or it should
        # go ahead and compute the regression and report a correlation and
        # total sum of squares of exactly 0.
        y = stats.linregress(X,ZERO)
        intercept = y[1]
        r = y[2]
        assert_almost_equal(intercept,0.0)
        assert_almost_equal(r,0.0)

    def test_regress_simple(self):
        # Regress a line with sinusoidal noise.
        x = np.linspace(0, 100, 100)
        y = 0.2 * np.linspace(0, 100, 100) + 10
        y += np.sin(np.linspace(0, 20, 100))

        res = stats.linregress(x, y)
        assert_almost_equal(res[4], 2.3957814497838803e-3)

    def test_regress_simple_onearg_rows(self):
        # Regress a line w sinusoidal noise, with a single input of shape (2, N).
        x = np.linspace(0, 100, 100)
        y = 0.2 * np.linspace(0, 100, 100) + 10
        y += np.sin(np.linspace(0, 20, 100))
        rows = np.vstack((x, y))

        res = stats.linregress(rows)
        assert_almost_equal(res[4], 2.3957814497838803e-3)

    def test_regress_simple_onearg_cols(self):
        x = np.linspace(0, 100, 100)
        y = 0.2 * np.linspace(0, 100, 100) + 10
        y += np.sin(np.linspace(0, 20, 100))
        cols = np.hstack((np.expand_dims(x, 1), np.expand_dims(y, 1)))

        res = stats.linregress(cols)
        assert_almost_equal(res[4], 2.3957814497838803e-3)

    def test_regress_shape_error(self):
        # Check that a single input argument to linregress with wrong shape
        # results in a ValueError.
        assert_raises(ValueError, stats.linregress, np.ones((3, 3)))

    def test_linregress(self):
        # compared with multivariate ols with pinv
        x = np.arange(11)
        y = np.arange(5,16)
        y[[(1),(-2)]] -= 1
        y[[(0),(-1)]] += 1

        res = (1.0, 5.0, 0.98229948625750, 7.45259691e-008, 0.063564172616372733)
        assert_array_almost_equal(stats.linregress(x,y),res,decimal=14)

    def test_regress_simple_negative_cor(self):
        # If the slope of the regression is negative the factor R tend to -1 not 1.
        # Sometimes rounding errors makes it < -1 leading to stderr being NaN
        a, n = 1e-71, 100000
        x = np.linspace(a, 2 * a, n)
        y = np.linspace(2 * a, a, n)
        stats.linregress(x, y)
        res = stats.linregress(x, y)
        assert_(res[2] >= -1)  # propagated numerical errors were not corrected
        assert_almost_equal(res[2], -1)  # perfect negative correlation case
        assert_(not np.isnan(res[4]))  # stderr should stay finite

    def test_linregress_result_attributes(self):
        # Regress a line with sinusoidal noise.
        x = np.linspace(0, 100, 100)
        y = 0.2 * np.linspace(0, 100, 100) + 10
        y += np.sin(np.linspace(0, 20, 100))

        res = stats.linregress(x, y)
        attributes = ('slope', 'intercept', 'rvalue', 'pvalue', 'stderr')
        check_named_results(res, attributes)


def test_theilslopes():
    # Basic slope test.
    slope, intercept, lower, upper = stats.theilslopes([0,1,1])
    assert_almost_equal(slope, 0.5)
    assert_almost_equal(intercept, 0.5)

    # Test of confidence intervals.
    x = [1, 2, 3, 4, 10, 12, 18]
    y = [9, 15, 19, 20, 45, 55, 78]
    slope, intercept, lower, upper = stats.theilslopes(y, x, 0.07)
    assert_almost_equal(slope, 4)
    assert_almost_equal(upper, 4.38, decimal=2)
    assert_almost_equal(lower, 3.71, decimal=2)


class TestHistogram(TestCase):
    # Tests that histogram works as it should, and keeps old behaviour
    #
    # what is untested:
    # - multidimensional arrays (since 'a' is ravel'd as the first line in the method)
    # - very large arrays
    # - Nans, Infs, empty and otherwise bad inputs

    # sample arrays to test the histogram with
    low_values = np.array([0.2, 0.3, 0.4, 0.5, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2],
                          dtype=float)  # 11 values
    high_range = np.array([2, 3, 4, 2, 21, 32, 78, 95, 65, 66, 66, 66, 66, 4],
                          dtype=float)  # 14 values
    low_range = np.array([2, 3, 3, 2, 3, 2.4, 2.1, 3.1, 2.9, 2.6, 2.7, 2.8, 2.2, 2.001],
                         dtype=float)  # 14 values
    few_values = np.array([2.0, 3.0, -1.0, 0.0], dtype=float)  # 4 values

    def test_simple(self):
        # Tests that each of the tests works as expected with default params
        #
        # basic tests, with expected results (no weighting)
        # results taken from the previous (slower) version of histogram
        basic_tests = ((self.low_values, (np.array([1., 1., 1., 2., 2.,
                                                     1., 1., 0., 1., 1.]),
                                          0.14444444444444446, 0.11111111111111112, 0)),
                       (self.high_range, (np.array([5., 0., 1., 1., 0.,
                                                     0., 5., 1., 0., 1.]),
                                          -3.1666666666666661, 10.333333333333332, 0)),
                       (self.low_range, (np.array([3., 1., 1., 1., 0., 1.,
                                                    1., 2., 3., 1.]),
                                         1.9388888888888889, 0.12222222222222223, 0)),
                       (self.few_values, (np.array([1., 0., 1., 0., 0., 0.,
                                                     0., 1., 0., 1.]),
                                          -1.2222222222222223, 0.44444444444444448, 0)),
                       )
        for inputs, expected_results in basic_tests:
            given_results = stats.histogram(inputs)
            assert_array_almost_equal(expected_results[0], given_results[0],
                                      decimal=2)
            for i in range(1, 4):
                assert_almost_equal(expected_results[i], given_results[i],
                                    decimal=2)

    def test_weighting(self):
        # Tests that weights give expected histograms

        # basic tests, with expected results, given a set of weights
        # weights used (first n are used for each test, where n is len of array) (14 values)
        weights = np.array([1., 3., 4.5, 0.1, -1.0, 0.0, 0.3, 7.0, 103.2, 2, 40, 0, 0, 1])
        # results taken from the numpy version of histogram
        basic_tests = ((self.low_values, (np.array([4.0, 0.0, 4.5, -0.9, 0.0,
                                                      0.3,110.2, 0.0, 0.0, 42.0]),
                                          0.2, 0.1, 0)),
                       (self.high_range, (np.array([9.6, 0., -1., 0., 0.,
                                                      0.,145.2, 0., 0.3, 7.]),
                                          2.0, 9.3, 0)),
                       (self.low_range, (np.array([2.4, 0., 0., 0., 0.,
                                                    2., 40., 0., 103.2, 13.5]),
                                         2.0, 0.11, 0)),
                       (self.few_values, (np.array([4.5, 0., 0.1, 0., 0., 0.,
                                                     0., 1., 0., 3.]),
                                          -1., 0.4, 0)),

                       )
        for inputs, expected_results in basic_tests:
            # use the first lot of weights for test
            # default limits given to reproduce output of numpy's test better
            given_results = stats.histogram(inputs, defaultlimits=(inputs.min(),
                                                                   inputs.max()),
                                            weights=weights[:len(inputs)])
            assert_array_almost_equal(expected_results[0], given_results[0],
                                      decimal=2)
            for i in range(1, 4):
                assert_almost_equal(expected_results[i], given_results[i],
                                    decimal=2)

    def test_reduced_bins(self):
        # Tests that reducing the number of bins produces expected results

        # basic tests, with expected results (no weighting),
        # except number of bins is halved to 5
        # results taken from the previous (slower) version of histogram
        basic_tests = ((self.low_values, (np.array([2., 3., 3., 1., 2.]),
                                          0.075000000000000011, 0.25, 0)),
                       (self.high_range, (np.array([5., 2., 0., 6., 1.]),
                                          -9.625, 23.25, 0)),
                       (self.low_range, (np.array([4., 2., 1., 3., 4.]),
                                         1.8625, 0.27500000000000002, 0)),
                       (self.few_values, (np.array([1., 1., 0., 1., 1.]),
                                          -1.5, 1.0, 0)),
                       )
        for inputs, expected_results in basic_tests:
            given_results = stats.histogram(inputs, numbins=5)
            assert_array_almost_equal(expected_results[0], given_results[0],
                                      decimal=2)
            for i in range(1, 4):
                assert_almost_equal(expected_results[i], given_results[i],
                                    decimal=2)

    def test_increased_bins(self):
        # Tests that increasing the number of bins produces expected results

        # basic tests, with expected results (no weighting),
        # except number of bins is double to 20
        # results taken from the previous (slower) version of histogram
        basic_tests = ((self.low_values, (np.array([1., 0., 1., 0., 1.,
                                                     0., 2., 0., 1., 0.,
                                                     1., 1., 0., 1., 0.,
                                                     0., 0., 1., 0., 1.]),
                                          0.1736842105263158, 0.052631578947368418, 0)),
                       (self.high_range, (np.array([5., 0., 0., 0., 1.,
                                                     0., 1., 0., 0., 0.,
                                                     0., 0., 0., 5., 0.,
                                                     0., 1., 0., 0., 1.]),
                                          -0.44736842105263142, 4.8947368421052628, 0)),
                       (self.low_range, (np.array([3., 0., 1., 1., 0., 0.,
                                                    0., 1., 0., 0., 1., 0.,
                                                    1., 0., 1., 0., 1., 3.,
                                                    0., 1.]),
                                         1.9710526315789474, 0.057894736842105263, 0)),
                       (self.few_values, (np.array([1., 0., 0., 0., 0., 1.,
                                                     0., 0., 0., 0., 0., 0.,
                                                     0., 0., 1., 0., 0., 0.,
                                                     0., 1.]),
                                          -1.1052631578947367, 0.21052631578947367, 0)),
                       )
        for inputs, expected_results in basic_tests:
            given_results = stats.histogram(inputs, numbins=20)
            assert_array_almost_equal(expected_results[0], given_results[0],
                                      decimal=2)
            for i in range(1, 4):
                assert_almost_equal(expected_results[i], given_results[i],
                                    decimal=2)

    def test_histogram_result_attributes(self):
        res = stats.histogram(self.low_range, numbins=20)
        attributes = ('count', 'lowerlimit', 'binsize', 'extrapoints')
        check_named_results(res, attributes)


def test_cumfreq():
    x = [1, 4, 2, 1, 3, 1]
    cumfreqs, lowlim, binsize, extrapoints = stats.cumfreq(x, numbins=4)
    assert_array_almost_equal(cumfreqs, np.array([3., 4., 5., 6.]))
    cumfreqs, lowlim, binsize, extrapoints = stats.cumfreq(x, numbins=4,
                                                      defaultreallimits=(1.5, 5))
    assert_(extrapoints == 3)

    # test for namedtuple attribute results
    attributes = ('cumcount', 'lowerlimit', 'binsize', 'extrapoints')
    res = stats.cumfreq(x, numbins=4, defaultreallimits=(1.5, 5))
    check_named_results(res, attributes)


def test_relfreq():
    a = np.array([1, 4, 2, 1, 3, 1])
    relfreqs, lowlim, binsize, extrapoints = stats.relfreq(a, numbins=4)
    assert_array_almost_equal(relfreqs,
                              array([0.5, 0.16666667, 0.16666667, 0.16666667]))

    # test for namedtuple attribute results
    attributes = ('frequency', 'lowerlimit', 'binsize', 'extrapoints')
    res = stats.relfreq(a, numbins=4)
    check_named_results(res, attributes)

    # check array_like input is accepted
    relfreqs2, lowlim, binsize, extrapoints = stats.relfreq([1, 4, 2, 1, 3, 1],
                                                            numbins=4)
    assert_array_almost_equal(relfreqs, relfreqs2)


class TestGMean(TestCase):

    def test_1D_list(self):
        a = (1,2,3,4)
        actual = stats.gmean(a)
        desired = power(1*2*3*4,1./4.)
        assert_almost_equal(actual, desired,decimal=14)

        desired1 = stats.gmean(a,axis=-1)
        assert_almost_equal(actual, desired1, decimal=14)

    def test_1D_array(self):
        a = array((1,2,3,4), float32)
        actual = stats.gmean(a)
        desired = power(1*2*3*4,1./4.)
        assert_almost_equal(actual, desired, decimal=7)

        desired1 = stats.gmean(a,axis=-1)
        assert_almost_equal(actual, desired1, decimal=7)

    def test_2D_array_default(self):
        a = array(((1,2,3,4),
                   (1,2,3,4),
                   (1,2,3,4)))
        actual = stats.gmean(a)
        desired = array((1,2,3,4))
        assert_array_almost_equal(actual, desired, decimal=14)

        desired1 = stats.gmean(a,axis=0)
        assert_array_almost_equal(actual, desired1, decimal=14)

    def test_2D_array_dim1(self):
        a = array(((1,2,3,4),
                   (1,2,3,4),
                   (1,2,3,4)))
        actual = stats.gmean(a, axis=1)
        v = power(1*2*3*4,1./4.)
        desired = array((v,v,v))
        assert_array_almost_equal(actual, desired, decimal=14)

    def test_large_values(self):
        a = array([1e100, 1e200, 1e300])
        actual = stats.gmean(a)
        assert_approx_equal(actual, 1e200, significant=14)


class TestHMean(TestCase):
    def test_1D_list(self):
        a = (1,2,3,4)
        actual = stats.hmean(a)
        desired = 4. / (1./1 + 1./2 + 1./3 + 1./4)
        assert_almost_equal(actual, desired, decimal=14)

        desired1 = stats.hmean(array(a),axis=-1)
        assert_almost_equal(actual, desired1, decimal=14)

    def test_1D_array(self):
        a = array((1,2,3,4), float64)
        actual = stats.hmean(a)
        desired = 4. / (1./1 + 1./2 + 1./3 + 1./4)
        assert_almost_equal(actual, desired, decimal=14)

        desired1 = stats.hmean(a,axis=-1)
        assert_almost_equal(actual, desired1, decimal=14)

    def test_2D_array_default(self):
        a = array(((1,2,3,4),
                   (1,2,3,4),
                   (1,2,3,4)))
        actual = stats.hmean(a)
        desired = array((1.,2.,3.,4.))
        assert_array_almost_equal(actual, desired, decimal=14)

        actual1 = stats.hmean(a,axis=0)
        assert_array_almost_equal(actual1, desired, decimal=14)

    def test_2D_array_dim1(self):
        a = array(((1,2,3,4),
                   (1,2,3,4),
                   (1,2,3,4)))

        v = 4. / (1./1 + 1./2 + 1./3 + 1./4)
        desired1 = array((v,v,v))
        actual1 = stats.hmean(a, axis=1)
        assert_array_almost_equal(actual1, desired1, decimal=14)


class TestScoreatpercentile(TestCase):
    def setUp(self):
        self.a1 = [3, 4, 5, 10, -3, -5, 6]
        self.a2 = [3, -6, -2, 8, 7, 4, 2, 1]
        self.a3 = [3., 4, 5, 10, -3, -5, -6, 7.0]

    def test_basic(self):
        x = arange(8) * 0.5
        assert_equal(stats.scoreatpercentile(x, 0), 0.)
        assert_equal(stats.scoreatpercentile(x, 100), 3.5)
        assert_equal(stats.scoreatpercentile(x, 50), 1.75)

    def test_fraction(self):
        scoreatperc = stats.scoreatpercentile

        # Test defaults
        assert_equal(scoreatperc(list(range(10)), 50), 4.5)
        assert_equal(scoreatperc(list(range(10)), 50, (2,7)), 4.5)
        assert_equal(scoreatperc(list(range(100)), 50, limit=(1, 8)), 4.5)
        assert_equal(scoreatperc(np.array([1, 10,100]), 50, (10,100)), 55)
        assert_equal(scoreatperc(np.array([1, 10,100]), 50, (1,10)), 5.5)

        # explicitly specify interpolation_method 'fraction' (the default)
        assert_equal(scoreatperc(list(range(10)), 50, interpolation_method='fraction'),
                     4.5)
        assert_equal(scoreatperc(list(range(10)), 50, limit=(2, 7),
                                 interpolation_method='fraction'),
                     4.5)
        assert_equal(scoreatperc(list(range(100)), 50, limit=(1, 8),
                                 interpolation_method='fraction'),
                     4.5)
        assert_equal(scoreatperc(np.array([1, 10,100]), 50, (10, 100),
                                 interpolation_method='fraction'),
                     55)
        assert_equal(scoreatperc(np.array([1, 10,100]), 50, (1,10),
                                 interpolation_method='fraction'),
                     5.5)

    def test_lower_higher(self):
        scoreatperc = stats.scoreatpercentile

        # interpolation_method 'lower'/'higher'
        assert_equal(scoreatperc(list(range(10)), 50,
                                 interpolation_method='lower'), 4)
        assert_equal(scoreatperc(list(range(10)), 50,
                                 interpolation_method='higher'), 5)
        assert_equal(scoreatperc(list(range(10)), 50, (2,7),
                                 interpolation_method='lower'), 4)
        assert_equal(scoreatperc(list(range(10)), 50, limit=(2,7),
                                 interpolation_method='higher'), 5)
        assert_equal(scoreatperc(list(range(100)), 50, (1,8),
                                 interpolation_method='lower'), 4)
        assert_equal(scoreatperc(list(range(100)), 50, (1,8),
                                 interpolation_method='higher'), 5)
        assert_equal(scoreatperc(np.array([1, 10, 100]), 50, (10, 100),
                                 interpolation_method='lower'), 10)
        assert_equal(scoreatperc(np.array([1, 10, 100]), 50, limit=(10, 100),
                                 interpolation_method='higher'), 100)
        assert_equal(scoreatperc(np.array([1, 10, 100]), 50, (1, 10),
                                 interpolation_method='lower'), 1)
        assert_equal(scoreatperc(np.array([1, 10, 100]), 50, limit=(1, 10),
                                 interpolation_method='higher'), 10)

    def test_sequence_per(self):
        x = arange(8) * 0.5
        expected = np.array([0, 3.5, 1.75])
        res = stats.scoreatpercentile(x, [0, 100, 50])
        assert_allclose(res, expected)
        assert_(isinstance(res, np.ndarray))
        # Test with ndarray.  Regression test for gh-2861
        assert_allclose(stats.scoreatpercentile(x, np.array([0, 100, 50])),
                        expected)
        # Also test combination of 2-D array, axis not None and array-like per
        res2 = stats.scoreatpercentile(np.arange(12).reshape((3,4)),
                                       np.array([0, 1, 100, 100]), axis=1)
        expected2 = array([[0, 4, 8],
                           [0.03, 4.03, 8.03],
                           [3, 7, 11],
                           [3, 7, 11]])
        assert_allclose(res2, expected2)

    def test_axis(self):
        scoreatperc = stats.scoreatpercentile
        x = arange(12).reshape(3, 4)

        assert_equal(scoreatperc(x, (25, 50, 100)), [2.75, 5.5, 11.0])

        r0 = [[2, 3, 4, 5], [4, 5, 6, 7], [8, 9, 10, 11]]
        assert_equal(scoreatperc(x, (25, 50, 100), axis=0), r0)

        r1 = [[0.75, 4.75, 8.75], [1.5, 5.5, 9.5], [3, 7, 11]]
        assert_equal(scoreatperc(x, (25, 50, 100), axis=1), r1)

        x = array([[1, 1, 1],
                   [1, 1, 1],
                   [4, 4, 3],
                   [1, 1, 1],
                   [1, 1, 1]])
        score = stats.scoreatpercentile(x, 50)
        assert_equal(score.shape, ())
        assert_equal(score, 1.0)
        score = stats.scoreatpercentile(x, 50, axis=0)
        assert_equal(score.shape, (3,))
        assert_equal(score, [1, 1, 1])

    def test_exception(self):
        assert_raises(ValueError, stats.scoreatpercentile, [1, 2], 56,
            interpolation_method='foobar')
        assert_raises(ValueError, stats.scoreatpercentile, [1], 101)
        assert_raises(ValueError, stats.scoreatpercentile, [1], -1)

    def test_empty(self):
        assert_equal(stats.scoreatpercentile([], 50), np.nan)
        assert_equal(stats.scoreatpercentile(np.array([[], []]), 50), np.nan)
        assert_equal(stats.scoreatpercentile([], [50, 99]), [np.nan, np.nan])


class TestItemfreq(object):
    a = [5, 7, 1, 2, 1, 5, 7] * 10
    b = [1, 2, 5, 7]

    def test_numeric_types(self):
        # Check itemfreq works for all dtypes (adapted from np.unique tests)
        def _check_itemfreq(dt):
            a = np.array(self.a, dt)
            v = stats.itemfreq(a)
            assert_array_equal(v[:, 0], [1, 2, 5, 7])
            assert_array_equal(v[:, 1], np.array([20, 10, 20, 20], dtype=dt))

        dtypes = [np.int32, np.int64, np.float32, np.float64,
                  np.complex64, np.complex128]
        for dt in dtypes:
            yield _check_itemfreq, dt

    def test_object_arrays(self):
        a, b = self.a, self.b
        dt = 'O'
        aa = np.empty(len(a), dt)
        aa[:] = a
        bb = np.empty(len(b), dt)
        bb[:] = b
        v = stats.itemfreq(aa)
        assert_array_equal(v[:, 0], bb)

    def test_structured_arrays(self):
        a, b = self.a, self.b
        dt = [('', 'i'), ('', 'i')]
        aa = np.array(list(zip(a, a)), dt)
        bb = np.array(list(zip(b, b)), dt)
        v = stats.itemfreq(aa)
        # Arrays don't compare equal because v[:,0] is object array
        assert_equal(tuple(v[2, 0]), tuple(bb[2]))


class TestMode(TestCase):
    def test_empty(self):
        vals, counts = stats.mode([])
        assert_equal(vals, np.array([]))
        assert_equal(counts, np.array([]))

    def test_scalar(self):
        vals, counts = stats.mode(4.)
        assert_equal(vals, np.array([4.]))
        assert_equal(counts, np.array([1]))

    def test_basic(self):
        data1 = [3, 5, 1, 10, 23, 3, 2, 6, 8, 6, 10, 6]
        vals = stats.mode(data1)
        assert_equal(vals[0][0], 6)
        assert_equal(vals[1][0], 3)

    def test_axes(self):
        data1 = [10, 10, 30, 40]
        data2 = [10, 10, 10, 10]
        data3 = [20, 10, 20, 20]
        data4 = [30, 30, 30, 30]
        data5 = [40, 30, 30, 30]
        arr = np.array([data1, data2, data3, data4, data5])

        vals = stats.mode(arr, axis=None)
        assert_equal(vals[0], np.array([30]))
        assert_equal(vals[1], np.array([8]))

        vals = stats.mode(arr, axis=0)
        assert_equal(vals[0], np.array([[10, 10, 30, 30]]))
        assert_equal(vals[1], np.array([[2, 3, 3, 2]]))

        vals = stats.mode(arr, axis=1)
        assert_equal(vals[0], np.array([[10], [10], [20], [30], [30]]))
        assert_equal(vals[1], np.array([[2], [4], [3], [4], [3]]))

    def test_strings(self):
        data1 = ['rain', 'showers', 'showers']
        vals = stats.mode(data1)
        assert_equal(vals[0][0], 'showers')
        assert_equal(vals[1][0], 2)

    @dec.knownfailureif(sys.version_info > (3,), 'numpy github issue 641')
    def test_mixed_objects(self):
        objects = [10, True, np.nan, 'hello', 10]
        arr = np.empty((5,), dtype=object)
        arr[:] = objects
        vals = stats.mode(arr)
        assert_equal(vals[0][0], 10)
        assert_equal(vals[1][0], 2)

    def test_objects(self):
        """Python objects must be sortable (le + eq) and have ne defined
        for np.unique to work. hash is for set.
        """
        class Point(object):
            def __init__(self, x):
                self.x = x

            def __eq__(self, other):
                return self.x == other.x

            def __ne__(self, other):
                return self.x != other.x

            def __lt__(self, other):
                return self.x < other.x

            def __hash__(self):
                return hash(self.x)

        points = [Point(x) for x in [1, 2, 3, 4, 3, 2, 2, 2]]
        arr = np.empty((8,), dtype=object)
        arr[:] = points
        assert len(set(points)) == 4
        assert_equal(np.unique(arr).shape, (4,))
        vals = stats.mode(arr)
        assert_equal(vals[0][0], Point(2))
        assert_equal(vals[1][0], 4)

    def test_mode_result_attributes(self):
        data1 = [3, 5, 1, 10, 23, 3, 2, 6, 8, 6, 10, 6]
        actual = stats.mode(data1)
        attributes = ('mode', 'count')
        check_named_results(actual, attributes)

class TestVariability(TestCase):

    testcase = [1,2,3,4]
    scalar_testcase = 4.

    def test_signaltonoise(self):
        # This is not in R, so used:
        #     mean(testcase, axis=0) / (sqrt(var(testcase) * 3/4))

        # y = stats.signaltonoise(self.shoes[0])
        # assert_approx_equal(y,4.5709967)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            y = stats.signaltonoise(self.testcase)
        assert_approx_equal(y, 2.236067977)

    def test_sem(self):
        # This is not in R, so used:
        #     sqrt(var(testcase)*3/4)/sqrt(3)

        # y = stats.sem(self.shoes[0])
        # assert_approx_equal(y,0.775177399)
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=RuntimeWarning)
            y = stats.sem(self.scalar_testcase)
        assert_(np.isnan(y))

        y = stats.sem(self.testcase)
        assert_approx_equal(y, 0.6454972244)
        n = len(self.testcase)
        assert_allclose(stats.sem(self.testcase, ddof=0) * np.sqrt(n/(n-2)),
                        stats.sem(self.testcase, ddof=2))

    def test_zmap(self):
        # not in R, so tested by using:
        #     (testcase[i] - mean(testcase, axis=0)) / sqrt(var(testcase) * 3/4)
        y = stats.zmap(self.testcase,self.testcase)
        desired = ([-1.3416407864999, -0.44721359549996, 0.44721359549996, 1.3416407864999])
        assert_array_almost_equal(desired,y,decimal=12)

    def test_zmap_axis(self):
        # Test use of 'axis' keyword in zmap.
        x = np.array([[0.0, 0.0, 1.0, 1.0],
                      [1.0, 1.0, 1.0, 2.0],
                      [2.0, 0.0, 2.0, 0.0]])

        t1 = 1.0/np.sqrt(2.0/3)
        t2 = np.sqrt(3.)/3
        t3 = np.sqrt(2.)

        z0 = stats.zmap(x, x, axis=0)
        z1 = stats.zmap(x, x, axis=1)

        z0_expected = [[-t1, -t3/2, -t3/2, 0.0],
                       [0.0, t3, -t3/2, t1],
                       [t1, -t3/2, t3, -t1]]
        z1_expected = [[-1.0, -1.0, 1.0, 1.0],
                       [-t2, -t2, -t2, np.sqrt(3.)],
                       [1.0, -1.0, 1.0, -1.0]]

        assert_array_almost_equal(z0, z0_expected)
        assert_array_almost_equal(z1, z1_expected)

    def test_zmap_ddof(self):
        # Test use of 'ddof' keyword in zmap.
        x = np.array([[0.0, 0.0, 1.0, 1.0],
                      [0.0, 1.0, 2.0, 3.0]])

        z = stats.zmap(x, x, axis=1, ddof=1)

        z0_expected = np.array([-0.5, -0.5, 0.5, 0.5])/(1.0/np.sqrt(3))
        z1_expected = np.array([-1.5, -0.5, 0.5, 1.5])/(np.sqrt(5./3))
        assert_array_almost_equal(z[0], z0_expected)
        assert_array_almost_equal(z[1], z1_expected)

    def test_zscore(self):
        # not in R, so tested by using:
        #    (testcase[i] - mean(testcase, axis=0)) / sqrt(var(testcase) * 3/4)
        y = stats.zscore(self.testcase)
        desired = ([-1.3416407864999, -0.44721359549996, 0.44721359549996, 1.3416407864999])
        assert_array_almost_equal(desired,y,decimal=12)

    def test_zscore_axis(self):
        # Test use of 'axis' keyword in zscore.
        x = np.array([[0.0, 0.0, 1.0, 1.0],
                      [1.0, 1.0, 1.0, 2.0],
                      [2.0, 0.0, 2.0, 0.0]])

        t1 = 1.0/np.sqrt(2.0/3)
        t2 = np.sqrt(3.)/3
        t3 = np.sqrt(2.)

        z0 = stats.zscore(x, axis=0)
        z1 = stats.zscore(x, axis=1)

        z0_expected = [[-t1, -t3/2, -t3/2, 0.0],
                       [0.0, t3, -t3/2, t1],
                       [t1, -t3/2, t3, -t1]]
        z1_expected = [[-1.0, -1.0, 1.0, 1.0],
                       [-t2, -t2, -t2, np.sqrt(3.)],
                       [1.0, -1.0, 1.0, -1.0]]

        assert_array_almost_equal(z0, z0_expected)
        assert_array_almost_equal(z1, z1_expected)

    def test_zscore_ddof(self):
        # Test use of 'ddof' keyword in zscore.
        x = np.array([[0.0, 0.0, 1.0, 1.0],
                      [0.0, 1.0, 2.0, 3.0]])

        z = stats.zscore(x, axis=1, ddof=1)

        z0_expected = np.array([-0.5, -0.5, 0.5, 0.5])/(1.0/np.sqrt(3))
        z1_expected = np.array([-1.5, -0.5, 0.5, 1.5])/(np.sqrt(5./3))
        assert_array_almost_equal(z[0], z0_expected)
        assert_array_almost_equal(z[1], z1_expected)


class TestMoments(TestCase):
    """
        Comparison numbers are found using R v.1.5.1
        note that length(testcase) = 4
        testmathworks comes from documentation for the
        Statistics Toolbox for Matlab and can be found at both
        http://www.mathworks.com/access/helpdesk/help/toolbox/stats/kurtosis.shtml
        http://www.mathworks.com/access/helpdesk/help/toolbox/stats/skewness.shtml
        Note that both test cases came from here.
    """
    testcase = [1,2,3,4]
    scalar_testcase = 4.
    np.random.seed(1234)
    testcase_moment_accuracy = np.random.rand(42)
    testmathworks = [1.165, 0.6268, 0.0751, 0.3516, -0.6965]

    def test_moment(self):
        # mean((testcase-mean(testcase))**power,axis=0),axis=0))**power))
        y = stats.moment(self.scalar_testcase)
        assert_approx_equal(y, 0.0)
        y = stats.moment(self.testcase, 1)
        assert_approx_equal(y, 0.0, 10)
        y = stats.moment(self.testcase, 2)
        assert_approx_equal(y, 1.25)
        y = stats.moment(self.testcase, 3)
        assert_approx_equal(y, 0.0)
        y = stats.moment(self.testcase, 4)
        assert_approx_equal(y, 2.5625)

    def test_variation(self):
        # variation = samplestd / mean
        y = stats.variation(self.scalar_testcase)
        assert_approx_equal(y, 0.0)
        y = stats.variation(self.testcase)
        assert_approx_equal(y, 0.44721359549996, 10)

    def test_skewness(self):
        # Scalar test case
        y = stats.skew(self.scalar_testcase)
        assert_approx_equal(y, 0.0)
        # sum((testmathworks-mean(testmathworks,axis=0))**3,axis=0) /
        #     ((sqrt(var(testmathworks)*4/5))**3)/5
        y = stats.skew(self.testmathworks)
        assert_approx_equal(y, -0.29322304336607, 10)
        y = stats.skew(self.testmathworks, bias=0)
        assert_approx_equal(y, -0.437111105023940, 10)
        y = stats.skew(self.testcase)
        assert_approx_equal(y, 0.0, 10)

    def test_skewness_scalar(self):
        # `skew` must return a scalar for 1-dim input
        assert_equal(stats.skew(arange(10)), 0.0)

    def test_kurtosis(self):
        # Scalar test case
        y = stats.kurtosis(self.scalar_testcase)
        assert_approx_equal(y, -3.0)
        #   sum((testcase-mean(testcase,axis=0))**4,axis=0)/((sqrt(var(testcase)*3/4))**4)/4
        #   sum((test2-mean(testmathworks,axis=0))**4,axis=0)/((sqrt(var(testmathworks)*4/5))**4)/5
        #   Set flags for axis = 0 and
        #   fisher=0 (Pearson's defn of kurtosis for compatiability with Matlab)
        y = stats.kurtosis(self.testmathworks, 0, fisher=0, bias=1)
        assert_approx_equal(y, 2.1658856802973, 10)

        # Note that MATLAB has confusing docs for the following case
        #  kurtosis(x,0) gives an unbiased estimate of Pearson's skewness
        #  kurtosis(x)  gives a biased estimate of Fisher's skewness (Pearson-3)
        #  The MATLAB docs imply that both should give Fisher's
        y = stats.kurtosis(self.testmathworks, fisher=0, bias=0)
        assert_approx_equal(y, 3.663542721189047, 10)
        y = stats.kurtosis(self.testcase, 0, 0)
        assert_approx_equal(y, 1.64)

    def test_kurtosis_array_scalar(self):
        assert_equal(type(stats.kurtosis([1,2,3])), float)

    def test_moment_accuracy(self):
        # 'moment' must have a small enough error compared to the slower
        #  but very accurate numpy.power() implementation.
        tc_no_mean = self.testcase_moment_accuracy - \
                     np.mean(self.testcase_moment_accuracy)
        assert_allclose(np.power(tc_no_mean, 42).mean(),
                            stats.moment(self.testcase_moment_accuracy, 42))


class TestThreshold(TestCase):
    def test_basic(self):
        a = [-1,2,3,4,5,-1,-2]
        assert_array_equal(stats.threshold(a),a)
        assert_array_equal(stats.threshold(a,3,None,0),
                           [0,0,3,4,5,0,0])
        assert_array_equal(stats.threshold(a,None,3,0),
                           [-1,2,3,0,0,-1,-2])
        assert_array_equal(stats.threshold(a,2,4,0),
                           [0,2,3,4,0,0,0])


class TestStudentTest(TestCase):
    X1 = np.array([-1, 0, 1])
    X2 = np.array([0, 1, 2])
    T1_0 = 0
    P1_0 = 1
    T1_1 = -1.732051
    P1_1 = 0.2254033
    T1_2 = -3.464102
    P1_2 = 0.0741799
    T2_0 = 1.732051
    P2_0 = 0.2254033

    def test_onesample(self):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=RuntimeWarning)
            t, p = stats.ttest_1samp(4., 3.)
        assert_(np.isnan(t))
        assert_(np.isnan(p))

        t, p = stats.ttest_1samp(self.X1, 0)

        assert_array_almost_equal(t, self.T1_0)
        assert_array_almost_equal(p, self.P1_0)

        res = stats.ttest_1samp(self.X1, 0)
        attributes = ('statistic', 'pvalue')
        check_named_results(res, attributes)

        t, p = stats.ttest_1samp(self.X2, 0)

        assert_array_almost_equal(t, self.T2_0)
        assert_array_almost_equal(p, self.P2_0)

        t, p = stats.ttest_1samp(self.X1, 1)

        assert_array_almost_equal(t, self.T1_1)
        assert_array_almost_equal(p, self.P1_1)

        t, p = stats.ttest_1samp(self.X1, 2)

        assert_array_almost_equal(t, self.T1_2)
        assert_array_almost_equal(p, self.P1_2)


def test_percentileofscore():
    pcos = stats.percentileofscore

    assert_equal(pcos([1,2,3,4,5,6,7,8,9,10],4), 40.0)

    for (kind, result) in [('mean', 35.0),
                           ('strict', 30.0),
                           ('weak', 40.0)]:
        yield assert_equal, pcos(np.arange(10) + 1,
                                                    4, kind=kind), \
                                                    result

    # multiple - 2
    for (kind, result) in [('rank', 45.0),
                           ('strict', 30.0),
                           ('weak', 50.0),
                           ('mean', 40.0)]:
        yield assert_equal, pcos([1,2,3,4,4,5,6,7,8,9],
                                                    4, kind=kind), \
                                                    result

    # multiple - 3
    assert_equal(pcos([1,2,3,4,4,4,5,6,7,8], 4), 50.0)
    for (kind, result) in [('rank', 50.0),
                           ('mean', 45.0),
                           ('strict', 30.0),
                           ('weak', 60.0)]:

        yield assert_equal, pcos([1,2,3,4,4,4,5,6,7,8],
                                                    4, kind=kind), \
                                                    result

    # missing
    for kind in ('rank', 'mean', 'strict', 'weak'):
        yield assert_equal, pcos([1,2,3,5,6,7,8,9,10,11],
                                                    4, kind=kind), \
                                                    30

    # larger numbers
    for (kind, result) in [('mean', 35.0),
                           ('strict', 30.0),
                           ('weak', 40.0)]:
        yield assert_equal, \
              pcos([10, 20, 30, 40, 50, 60, 70, 80, 90, 100], 40,
                   kind=kind), result

    for (kind, result) in [('mean', 45.0),
                           ('strict', 30.0),
                           ('weak', 60.0)]:
        yield assert_equal, \
              pcos([10, 20, 30, 40, 40, 40, 50, 60, 70, 80],
                   40, kind=kind), result

    for kind in ('rank', 'mean', 'strict', 'weak'):
        yield assert_equal, \
              pcos([10, 20, 30, 50, 60, 70, 80, 90, 100, 110],
                   40, kind=kind), 30.0

    # boundaries
    for (kind, result) in [('rank', 10.0),
                           ('mean', 5.0),
                           ('strict', 0.0),
                           ('weak', 10.0)]:
        yield assert_equal, \
              pcos([10, 20, 30, 50, 60, 70, 80, 90, 100, 110],
                   10, kind=kind), result

    for (kind, result) in [('rank', 100.0),
                           ('mean', 95.0),
                           ('strict', 90.0),
                           ('weak', 100.0)]:
        yield assert_equal, \
              pcos([10, 20, 30, 50, 60, 70, 80, 90, 100, 110],
                   110, kind=kind), result

    # out of bounds
    for (kind, score, result) in [('rank', 200, 100.0),
                                  ('mean', 200, 100.0),
                                  ('mean', 0, 0.0)]:
        yield assert_equal, \
              pcos([10, 20, 30, 50, 60, 70, 80, 90, 100, 110],
                   score, kind=kind), result

    assert_raises(ValueError, pcos, [1, 2, 3, 3, 4], 3, kind='unrecognized')


PowerDivCase = namedtuple('Case', ['f_obs', 'f_exp', 'ddof', 'axis',
                                   'chi2',     # Pearson's
                                   'log',      # G-test (log-likelihood)
                                   'mod_log',  # Modified log-likelihood
                                   'cr',       # Cressie-Read (lambda=2/3)
                                   ])

# The details of the first two elements in power_div_1d_cases are used
# in a test in TestPowerDivergence.  Check that code before making
# any changes here.
power_div_1d_cases = [
    # Use the default f_exp.
    PowerDivCase(f_obs=[4, 8, 12, 8], f_exp=None, ddof=0, axis=None,
                 chi2=4,
                 log=2*(4*np.log(4/8) + 12*np.log(12/8)),
                 mod_log=2*(8*np.log(8/4) + 8*np.log(8/12)),
                 cr=(4*((4/8)**(2/3) - 1) + 12*((12/8)**(2/3) - 1))/(5/9)),
    # Give a non-uniform f_exp.
    PowerDivCase(f_obs=[4, 8, 12, 8], f_exp=[2, 16, 12, 2], ddof=0, axis=None,
                 chi2=24,
                 log=2*(4*np.log(4/2) + 8*np.log(8/16) + 8*np.log(8/2)),
                 mod_log=2*(2*np.log(2/4) + 16*np.log(16/8) + 2*np.log(2/8)),
                 cr=(4*((4/2)**(2/3) - 1) + 8*((8/16)**(2/3) - 1) +
                     8*((8/2)**(2/3) - 1))/(5/9)),
    # f_exp is a scalar.
    PowerDivCase(f_obs=[4, 8, 12, 8], f_exp=8, ddof=0, axis=None,
                 chi2=4,
                 log=2*(4*np.log(4/8) + 12*np.log(12/8)),
                 mod_log=2*(8*np.log(8/4) + 8*np.log(8/12)),
                 cr=(4*((4/8)**(2/3) - 1) + 12*((12/8)**(2/3) - 1))/(5/9)),
    # f_exp equal to f_obs.
    PowerDivCase(f_obs=[3, 5, 7, 9], f_exp=[3, 5, 7, 9], ddof=0, axis=0,
                 chi2=0, log=0, mod_log=0, cr=0),
]


power_div_empty_cases = [
    # Shape is (0,)--a data set with length 0.  The computed
    # test statistic should be 0.
    PowerDivCase(f_obs=[],
                 f_exp=None, ddof=0, axis=0,
                 chi2=0, log=0, mod_log=0, cr=0),
    # Shape is (0, 3).  This is 3 data sets, but each data set has
    # length 0, so the computed test statistic should be [0, 0, 0].
    PowerDivCase(f_obs=np.array([[],[],[]]).T,
                 f_exp=None, ddof=0, axis=0,
                 chi2=[0, 0, 0],
                 log=[0, 0, 0],
                 mod_log=[0, 0, 0],
                 cr=[0, 0, 0]),
    # Shape is (3, 0).  This represents an empty collection of
    # data sets in which each data set has length 3.  The test
    # statistic should be an empty array.
    PowerDivCase(f_obs=np.array([[],[],[]]),
                 f_exp=None, ddof=0, axis=0,
                 chi2=[],
                 log=[],
                 mod_log=[],
                 cr=[]),
]


class TestPowerDivergence(object):

    def check_power_divergence(self, f_obs, f_exp, ddof, axis, lambda_,
                               expected_stat):
        f_obs = np.asarray(f_obs)
        if axis is None:
            num_obs = f_obs.size
        else:
            b = np.broadcast(f_obs, f_exp)
            num_obs = b.shape[axis]
        stat, p = stats.power_divergence(f_obs=f_obs, f_exp=f_exp, ddof=ddof,
                                         axis=axis, lambda_=lambda_)
        assert_allclose(stat, expected_stat)

        if lambda_ == 1 or lambda_ == "pearson":
            # Also test stats.chisquare.
            stat, p = stats.chisquare(f_obs=f_obs, f_exp=f_exp, ddof=ddof,
                                      axis=axis)
            assert_allclose(stat, expected_stat)

        ddof = np.asarray(ddof)
        expected_p = stats.chisqprob(expected_stat, num_obs - 1 - ddof)
        assert_allclose(p, expected_p)

    def test_basic(self):
        for case in power_div_1d_cases:
            yield (self.check_power_divergence,
                   case.f_obs, case.f_exp, case.ddof, case.axis,
                   None, case.chi2)
            yield (self.check_power_divergence,
                   case.f_obs, case.f_exp, case.ddof, case.axis,
                   "pearson", case.chi2)
            yield (self.check_power_divergence,
                   case.f_obs, case.f_exp, case.ddof, case.axis,
                   1, case.chi2)
            yield (self.check_power_divergence,
                   case.f_obs, case.f_exp, case.ddof, case.axis,
                   "log-likelihood", case.log)
            yield (self.check_power_divergence,
                   case.f_obs, case.f_exp, case.ddof, case.axis,
                   "mod-log-likelihood", case.mod_log)
            yield (self.check_power_divergence,
                   case.f_obs, case.f_exp, case.ddof, case.axis,
                   "cressie-read", case.cr)
            yield (self.check_power_divergence,
                   case.f_obs, case.f_exp, case.ddof, case.axis,
                   2/3, case.cr)

    def test_basic_masked(self):
        for case in power_div_1d_cases:
            mobs = np.ma.array(case.f_obs)
            yield (self.check_power_divergence,
                   mobs, case.f_exp, case.ddof, case.axis,
                   None, case.chi2)
            yield (self.check_power_divergence,
                   mobs, case.f_exp, case.ddof, case.axis,
                   "pearson", case.chi2)
            yield (self.check_power_divergence,
                   mobs, case.f_exp, case.ddof, case.axis,
                   1, case.chi2)
            yield (self.check_power_divergence,
                   mobs, case.f_exp, case.ddof, case.axis,
                   "log-likelihood", case.log)
            yield (self.check_power_divergence,
                   mobs, case.f_exp, case.ddof, case.axis,
                   "mod-log-likelihood", case.mod_log)
            yield (self.check_power_divergence,
                   mobs, case.f_exp, case.ddof, case.axis,
                   "cressie-read", case.cr)
            yield (self.check_power_divergence,
                   mobs, case.f_exp, case.ddof, case.axis,
                   2/3, case.cr)

    def test_axis(self):
        case0 = power_div_1d_cases[0]
        case1 = power_div_1d_cases[1]
        f_obs = np.vstack((case0.f_obs, case1.f_obs))
        f_exp = np.vstack((np.ones_like(case0.f_obs)*np.mean(case0.f_obs),
                           case1.f_exp))
        # Check the four computational code paths in power_divergence
        # using a 2D array with axis=1.
        yield (self.check_power_divergence,
               f_obs, f_exp, 0, 1,
               "pearson", [case0.chi2, case1.chi2])
        yield (self.check_power_divergence,
               f_obs, f_exp, 0, 1,
               "log-likelihood", [case0.log, case1.log])
        yield (self.check_power_divergence,
               f_obs, f_exp, 0, 1,
               "mod-log-likelihood", [case0.mod_log, case1.mod_log])
        yield (self.check_power_divergence,
               f_obs, f_exp, 0, 1,
               "cressie-read", [case0.cr, case1.cr])
        # Reshape case0.f_obs to shape (2,2), and use axis=None.
        # The result should be the same.
        yield (self.check_power_divergence,
               np.array(case0.f_obs).reshape(2, 2), None, 0, None,
               "pearson", case0.chi2)

    def test_ddof_broadcasting(self):
        # Test that ddof broadcasts correctly.
        # ddof does not affect the test statistic.  It is broadcast
        # with the computed test statistic for the computation of
        # the p value.

        case0 = power_div_1d_cases[0]
        case1 = power_div_1d_cases[1]
        # Create 4x2 arrays of observed and expected frequencies.
        f_obs = np.vstack((case0.f_obs, case1.f_obs)).T
        f_exp = np.vstack((np.ones_like(case0.f_obs)*np.mean(case0.f_obs),
                           case1.f_exp)).T

        expected_chi2 = [case0.chi2, case1.chi2]

        # ddof has shape (2, 1).  This is broadcast with the computed
        # statistic, so p will have shape (2,2).
        ddof = np.array([[0], [1]])

        stat, p = stats.power_divergence(f_obs, f_exp, ddof=ddof)
        assert_allclose(stat, expected_chi2)

        # Compute the p values separately, passing in scalars for ddof.
        stat0, p0 = stats.power_divergence(f_obs, f_exp, ddof=ddof[0,0])
        stat1, p1 = stats.power_divergence(f_obs, f_exp, ddof=ddof[1,0])

        assert_array_equal(p, np.vstack((p0, p1)))

    def test_empty_cases(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            for case in power_div_empty_cases:
                yield (self.check_power_divergence,
                       case.f_obs, case.f_exp, case.ddof, case.axis,
                       "pearson", case.chi2)
                yield (self.check_power_divergence,
                       case.f_obs, case.f_exp, case.ddof, case.axis,
                       "log-likelihood", case.log)
                yield (self.check_power_divergence,
                       case.f_obs, case.f_exp, case.ddof, case.axis,
                       "mod-log-likelihood", case.mod_log)
                yield (self.check_power_divergence,
                       case.f_obs, case.f_exp, case.ddof, case.axis,
                       "cressie-read", case.cr)

    def test_power_divergence_result_attributes(self):
        f_obs = power_div_1d_cases[0].f_obs
        f_exp = power_div_1d_cases[0].f_exp
        ddof = power_div_1d_cases[0].ddof
        axis = power_div_1d_cases[0].axis

        res = stats.power_divergence(f_obs=f_obs, f_exp=f_exp, ddof=ddof,
                                     axis=axis, lambda_="pearson")
        attributes = ('statistic', 'pvalue')
        check_named_results(res, attributes)


def test_chisquare_masked_arrays():
    # Test masked arrays.
    obs = np.array([[8, 8, 16, 32, -1], [-1, -1, 3, 4, 5]]).T
    mask = np.array([[0, 0, 0, 0, 1], [1, 1, 0, 0, 0]]).T
    mobs = np.ma.masked_array(obs, mask)
    expected_chisq = np.array([24.0, 0.5])
    expected_g = np.array([2*(2*8*np.log(0.5) + 32*np.log(2.0)),
                           2*(3*np.log(0.75) + 5*np.log(1.25))])

    chisq, p = stats.chisquare(mobs)
    mat.assert_array_equal(chisq, expected_chisq)
    mat.assert_array_almost_equal(p, stats.chisqprob(expected_chisq,
                                                     mobs.count(axis=0) - 1))

    g, p = stats.power_divergence(mobs, lambda_='log-likelihood')
    mat.assert_array_almost_equal(g, expected_g, decimal=15)
    mat.assert_array_almost_equal(p, stats.chisqprob(expected_g,
                                                     mobs.count(axis=0) - 1))

    chisq, p = stats.chisquare(mobs.T, axis=1)
    mat.assert_array_equal(chisq, expected_chisq)
    mat.assert_array_almost_equal(p,
                                  stats.chisqprob(expected_chisq,
                                                  mobs.T.count(axis=1) - 1))
    g, p = stats.power_divergence(mobs.T, axis=1, lambda_="log-likelihood")
    mat.assert_array_almost_equal(g, expected_g, decimal=15)
    mat.assert_array_almost_equal(p, stats.chisqprob(expected_g,
                                                     mobs.count(axis=0) - 1))

    obs1 = np.ma.array([3, 5, 6, 99, 10], mask=[0, 0, 0, 1, 0])
    exp1 = np.ma.array([2, 4, 8, 10, 99], mask=[0, 0, 0, 0, 1])
    chi2, p = stats.chisquare(obs1, f_exp=exp1)
    # Because of the mask at index 3 of obs1 and at index 4 of exp1,
    # only the first three elements are included in the calculation
    # of the statistic.
    mat.assert_array_equal(chi2, 1/2 + 1/4 + 4/8)

    # When axis=None, the two values should have type np.float64.
    chisq, p = stats.chisquare(np.ma.array([1,2,3]), axis=None)
    assert_(isinstance(chisq, np.float64))
    assert_(isinstance(p, np.float64))
    assert_equal(chisq, 1.0)
    assert_almost_equal(p, stats.chisqprob(1.0, 2))

    # Empty arrays:
    # A data set with length 0 returns a masked scalar.
    with np.errstate(invalid='ignore'):
        with warnings.catch_warnings(record=True):
            chisq, p = stats.chisquare(np.ma.array([]))
    assert_(isinstance(chisq, np.ma.MaskedArray))
    assert_equal(chisq.shape, ())
    assert_(chisq.mask)

    empty3 = np.ma.array([[],[],[]])

    # empty3 is a collection of 0 data sets (whose lengths would be 3, if
    # there were any), so the return value is an array with length 0.
    chisq, p = stats.chisquare(empty3)
    assert_(isinstance(chisq, np.ma.MaskedArray))
    mat.assert_array_equal(chisq, [])

    # empty3.T is an array containing 3 data sets, each with length 0,
    # so an array of size (3,) is returned, with all values masked.
    with np.errstate(invalid='ignore'):
        with warnings.catch_warnings(record=True):
            chisq, p = stats.chisquare(empty3.T)
    assert_(isinstance(chisq, np.ma.MaskedArray))
    assert_equal(chisq.shape, (3,))
    assert_(np.all(chisq.mask))


def test_power_divergence_against_cressie_read_data():
    # Test stats.power_divergence against tables 4 and 5 from
    # Cressie and Read, "Multimonial Goodness-of-Fit Tests",
    # J. R. Statist. Soc. B (1984), Vol 46, No. 3, pp. 440-464.
    # This tests the calculation for several values of lambda.

    # `table4` holds just the second and third columns from Table 4.
    table4 = np.array([
        # observed, expected,
        15, 15.171,
        11, 13.952,
        14, 12.831,
        17, 11.800,
        5, 10.852,
        11, 9.9796,
        10, 9.1777,
        4, 8.4402,
        8, 7.7620,
        10, 7.1383,
        7, 6.5647,
        9, 6.0371,
        11, 5.5520,
        3, 5.1059,
        6, 4.6956,
        1, 4.3183,
        1, 3.9713,
        4, 3.6522,
        ]).reshape(-1, 2)
    table5 = np.array([
        # lambda, statistic
        -10.0, 72.2e3,
        -5.0, 28.9e1,
        -3.0, 65.6,
        -2.0, 40.6,
        -1.5, 34.0,
        -1.0, 29.5,
        -0.5, 26.5,
        0.0, 24.6,
        0.5, 23.4,
        0.67, 23.1,
        1.0, 22.7,
        1.5, 22.6,
        2.0, 22.9,
        3.0, 24.8,
        5.0, 35.5,
        10.0, 21.4e1,
        ]).reshape(-1, 2)

    for lambda_, expected_stat in table5:
        stat, p = stats.power_divergence(table4[:,0], table4[:,1],
                                         lambda_=lambda_)
        assert_allclose(stat, expected_stat, rtol=5e-3)


def test_friedmanchisquare():
    # see ticket:113
    # verified with matlab and R
    # From Demsar "Statistical Comparisons of Classifiers over Multiple Data Sets"
    # 2006, Xf=9.28 (no tie handling, tie corrected Xf >=9.28)
    x1 = [array([0.763, 0.599, 0.954, 0.628, 0.882, 0.936, 0.661, 0.583,
                 0.775, 1.0, 0.94, 0.619, 0.972, 0.957]),
          array([0.768, 0.591, 0.971, 0.661, 0.888, 0.931, 0.668, 0.583,
                 0.838, 1.0, 0.962, 0.666, 0.981, 0.978]),
          array([0.771, 0.590, 0.968, 0.654, 0.886, 0.916, 0.609, 0.563,
                 0.866, 1.0, 0.965, 0.614, 0.9751, 0.946]),
          array([0.798, 0.569, 0.967, 0.657, 0.898, 0.931, 0.685, 0.625,
                 0.875, 1.0, 0.962, 0.669, 0.975, 0.970])]

    # From "Bioestadistica para las ciencias de la salud" Xf=18.95 p<0.001:
    x2 = [array([4,3,5,3,5,3,2,5,4,4,4,3]),
          array([2,2,1,2,3,1,2,3,2,1,1,3]),
          array([2,4,3,3,4,3,3,4,4,1,2,1]),
          array([3,5,4,3,4,4,3,3,3,4,4,4])]

    # From Jerrorl H. Zar, "Biostatistical Analysis"(example 12.6), Xf=10.68, 0.005 < p < 0.01:
    # Probability from this example is inexact using Chisquare aproximation of Friedman Chisquare.
    x3 = [array([7.0,9.9,8.5,5.1,10.3]),
          array([5.3,5.7,4.7,3.5,7.7]),
          array([4.9,7.6,5.5,2.8,8.4]),
          array([8.8,8.9,8.1,3.3,9.1])]

    assert_array_almost_equal(stats.friedmanchisquare(x1[0],x1[1],x1[2],x1[3]),
                              (10.2283464566929, 0.0167215803284414))
    assert_array_almost_equal(stats.friedmanchisquare(x2[0],x2[1],x2[2],x2[3]),
                              (18.9428571428571, 0.000280938375189499))
    assert_array_almost_equal(stats.friedmanchisquare(x3[0],x3[1],x3[2],x3[3]),
                              (10.68, 0.0135882729582176))
    np.testing.assert_raises(ValueError, stats.friedmanchisquare,x3[0],x3[1])

    # test for namedtuple attribute results
    attributes = ('statistic', 'pvalue')
    res = stats.friedmanchisquare(*x1)
    check_named_results(res, attributes)

    # test using mstats
    assert_array_almost_equal(stats.mstats.friedmanchisquare(x1[0],x1[1],x1[2],x1[3]),
                              (10.2283464566929, 0.0167215803284414))
    # the following fails
    # assert_array_almost_equal(stats.mstats.friedmanchisquare(x2[0],x2[1],x2[2],x2[3]),
    #                           (18.9428571428571, 0.000280938375189499))
    assert_array_almost_equal(stats.mstats.friedmanchisquare(x3[0],x3[1],x3[2],x3[3]),
                              (10.68, 0.0135882729582176))
    np.testing.assert_raises(ValueError,stats.mstats.friedmanchisquare,x3[0],x3[1])


def test_kstest():
    # from numpy.testing import assert_almost_equal

    # comparing with values from R
    x = np.linspace(-1,1,9)
    D,p = stats.kstest(x,'norm')
    assert_almost_equal(D, 0.15865525393145705, 12)
    assert_almost_equal(p, 0.95164069201518386, 1)

    x = np.linspace(-15,15,9)
    D,p = stats.kstest(x,'norm')
    assert_almost_equal(D, 0.44435602715924361, 15)
    assert_almost_equal(p, 0.038850140086788665, 8)

    # test for namedtuple attribute results
    attributes = ('statistic', 'pvalue')
    res = stats.kstest(x, 'norm')
    check_named_results(res, attributes)

    # the following tests rely on deterministicaly replicated rvs
    np.random.seed(987654321)
    x = stats.norm.rvs(loc=0.2, size=100)
    D,p = stats.kstest(x, 'norm', mode='asymp')
    assert_almost_equal(D, 0.12464329735846891, 15)
    assert_almost_equal(p, 0.089444888711820769, 15)
    assert_almost_equal(np.array(stats.kstest(x, 'norm', mode='asymp')),
                np.array((0.12464329735846891, 0.089444888711820769)), 15)
    assert_almost_equal(np.array(stats.kstest(x,'norm', alternative='less')),
                np.array((0.12464329735846891, 0.040989164077641749)), 15)
    # this 'greater' test fails with precision of decimal=14
    assert_almost_equal(np.array(stats.kstest(x,'norm', alternative='greater')),
                np.array((0.0072115233216310994, 0.98531158590396228)), 12)

    # missing: no test that uses *args


def test_ks_2samp():
    # exact small sample solution
    data1 = np.array([1.0,2.0])
    data2 = np.array([1.0,2.0,3.0])
    assert_almost_equal(np.array(stats.ks_2samp(data1+0.01,data2)),
                np.array((0.33333333333333337, 0.99062316386915694)))
    assert_almost_equal(np.array(stats.ks_2samp(data1-0.01,data2)),
                np.array((0.66666666666666674, 0.42490954988801982)))
    # these can also be verified graphically
    assert_almost_equal(
        np.array(stats.ks_2samp(np.linspace(1,100,100),
                              np.linspace(1,100,100)+2+0.1)),
        np.array((0.030000000000000027, 0.99999999996005062)))
    assert_almost_equal(
        np.array(stats.ks_2samp(np.linspace(1,100,100),
                              np.linspace(1,100,100)+2-0.1)),
        np.array((0.020000000000000018, 0.99999999999999933)))
    # these are just regression tests
    assert_almost_equal(
        np.array(stats.ks_2samp(np.linspace(1,100,100),
                              np.linspace(1,100,110)+20.1)),
        np.array((0.21090909090909091, 0.015880386730710221)))
    assert_almost_equal(
        np.array(stats.ks_2samp(np.linspace(1,100,100),
                              np.linspace(1,100,110)+20-0.1)),
        np.array((0.20818181818181825, 0.017981441789762638)))

    # test for namedtuple attribute results
    attributes = ('statistic', 'pvalue')
    res = stats.ks_2samp(data1 - 0.01, data2)
    check_named_results(res, attributes)


def test_ttest_rel():
    # regression test
    tr,pr = 0.81248591389165692, 0.41846234511362157
    tpr = ([tr,-tr],[pr,pr])

    rvs1 = np.linspace(1,100,100)
    rvs2 = np.linspace(1.01,99.989,100)
    rvs1_2D = np.array([np.linspace(1,100,100), np.linspace(1.01,99.989,100)])
    rvs2_2D = np.array([np.linspace(1.01,99.989,100), np.linspace(1,100,100)])

    t,p = stats.ttest_rel(rvs1, rvs2, axis=0)
    assert_array_almost_equal([t,p],(tr,pr))
    t,p = stats.ttest_rel(rvs1_2D.T, rvs2_2D.T, axis=0)
    assert_array_almost_equal([t,p],tpr)
    t,p = stats.ttest_rel(rvs1_2D, rvs2_2D, axis=1)
    assert_array_almost_equal([t,p],tpr)

    # test scalars
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        t, p = stats.ttest_rel(4., 3.)
    assert_(np.isnan(t))
    assert_(np.isnan(p))

    # test for namedtuple attribute results
    attributes = ('statistic', 'pvalue')
    res = stats.ttest_rel(rvs1, rvs2, axis=0)
    check_named_results(res, attributes)

    # test on 3 dimensions
    rvs1_3D = np.dstack([rvs1_2D,rvs1_2D,rvs1_2D])
    rvs2_3D = np.dstack([rvs2_2D,rvs2_2D,rvs2_2D])
    t,p = stats.ttest_rel(rvs1_3D, rvs2_3D, axis=1)
    assert_array_almost_equal(np.abs(t), tr)
    assert_array_almost_equal(np.abs(p), pr)
    assert_equal(t.shape, (2, 3))

    t,p = stats.ttest_rel(np.rollaxis(rvs1_3D,2), np.rollaxis(rvs2_3D,2), axis=2)
    assert_array_almost_equal(np.abs(t), tr)
    assert_array_almost_equal(np.abs(p), pr)
    assert_equal(t.shape, (3, 2))

    olderr = np.seterr(all='ignore')
    try:
        # test zero division problem
        t,p = stats.ttest_rel([0,0,0],[1,1,1])
        assert_equal((np.abs(t),p), (np.inf, 0))
        assert_equal(stats.ttest_rel([0,0,0], [0,0,0]), (np.nan, np.nan))

        # check that nan in input array result in nan output
        anan = np.array([[1, np.nan], [-1, 1]])
        assert_equal(stats.ttest_rel(anan, np.zeros((2, 2))),
                     ([0, np.nan], [1, np.nan]))
    finally:
        np.seterr(**olderr)

    # test incorrect input shape raise an error
    x = np.arange(24)
    assert_raises(ValueError, stats.ttest_rel, x.reshape((8, 3)),
                  x.reshape((2, 3, 4)))


def _desc_stats(x1, x2, axis=0):
    def _stats(x, axis=0):
        x = np.asarray(x)
        mu = np.mean(x, axis=axis)
        std = np.std(x, axis=axis, ddof=1)
        nobs = x.shape[axis]
        return mu, std, nobs
    return _stats(x1, axis) + _stats(x2, axis)


def test_ttest_ind():
    # regression test
    tr = 1.0912746897927283
    pr = 0.27647818616351882
    tpr = ([tr,-tr],[pr,pr])

    rvs2 = np.linspace(1,100,100)
    rvs1 = np.linspace(5,105,100)
    rvs1_2D = np.array([rvs1, rvs2])
    rvs2_2D = np.array([rvs2, rvs1])

    t,p = stats.ttest_ind(rvs1, rvs2, axis=0)
    assert_array_almost_equal([t,p],(tr,pr))
    # test from_stats API
    assert_array_almost_equal(stats.ttest_ind_from_stats(*_desc_stats(rvs1,
                                                                      rvs2)),
                              [t, p])
    t,p = stats.ttest_ind(rvs1_2D.T, rvs2_2D.T, axis=0)
    assert_array_almost_equal([t,p],tpr)
    args = _desc_stats(rvs1_2D.T, rvs2_2D.T)
    assert_array_almost_equal(stats.ttest_ind_from_stats(*args),
                              [t, p])
    t,p = stats.ttest_ind(rvs1_2D, rvs2_2D, axis=1)
    assert_array_almost_equal([t,p],tpr)
    args = _desc_stats(rvs1_2D, rvs2_2D, axis=1)
    assert_array_almost_equal(stats.ttest_ind_from_stats(*args),
                              [t, p])

    # test scalars
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        t, p = stats.ttest_ind(4., 3.)
    assert_(np.isnan(t))
    assert_(np.isnan(p))

    # test on 3 dimensions
    rvs1_3D = np.dstack([rvs1_2D,rvs1_2D,rvs1_2D])
    rvs2_3D = np.dstack([rvs2_2D,rvs2_2D,rvs2_2D])
    t,p = stats.ttest_ind(rvs1_3D, rvs2_3D, axis=1)
    assert_almost_equal(np.abs(t), np.abs(tr))
    assert_array_almost_equal(np.abs(p), pr)
    assert_equal(t.shape, (2, 3))

    t,p = stats.ttest_ind(np.rollaxis(rvs1_3D,2), np.rollaxis(rvs2_3D,2), axis=2)
    assert_array_almost_equal(np.abs(t), np.abs(tr))
    assert_array_almost_equal(np.abs(p), pr)
    assert_equal(t.shape, (3, 2))

    olderr = np.seterr(all='ignore')
    try:
        # test zero division problem
        t,p = stats.ttest_ind([0,0,0],[1,1,1])
        assert_equal((np.abs(t),p), (np.inf, 0))
        assert_equal(stats.ttest_ind([0,0,0], [0,0,0]), (np.nan, np.nan))

        # check that nan in input array result in nan output
        anan = np.array([[1,np.nan],[-1,1]])
        assert_equal(stats.ttest_ind(anan, np.zeros((2,2))),([0, np.nan], [1,np.nan]))
    finally:
        np.seterr(**olderr)


def test_ttest_ind_with_uneq_var():
    # check vs. R
    a = (1, 2, 3)
    b = (1.1, 2.9, 4.2)
    pr = 0.53619490753126731
    tr = -0.68649512735572582
    t, p = stats.ttest_ind(a, b, equal_var=False)
    assert_array_almost_equal([t,p], [tr, pr])
    # test from desc stats API
    assert_array_almost_equal(stats.ttest_ind_from_stats(*_desc_stats(a, b),
                                                         equal_var=False),
                              [t, p])

    a = (1, 2, 3, 4)
    pr = 0.84354139131608286
    tr = -0.2108663315950719
    t, p = stats.ttest_ind(a, b, equal_var=False)
    assert_array_almost_equal([t,p], [tr, pr])
    assert_array_almost_equal(stats.ttest_ind_from_stats(*_desc_stats(a, b),
                                                         equal_var=False),
                              [t, p])

    # regression test
    tr = 1.0912746897927283
    tr_uneq_n = 0.66745638708050492
    pr = 0.27647831993021388
    pr_uneq_n = 0.50873585065616544
    tpr = ([tr,-tr],[pr,pr])

    rvs3 = np.linspace(1,100, 25)
    rvs2 = np.linspace(1,100,100)
    rvs1 = np.linspace(5,105,100)
    rvs1_2D = np.array([rvs1, rvs2])
    rvs2_2D = np.array([rvs2, rvs1])

    t,p = stats.ttest_ind(rvs1, rvs2, axis=0, equal_var=False)
    assert_array_almost_equal([t,p],(tr,pr))
    assert_array_almost_equal(stats.ttest_ind_from_stats(*_desc_stats(rvs1,
                                                                      rvs2),
                                                         equal_var=False),
                              (t, p))

    t,p = stats.ttest_ind(rvs1, rvs3, axis=0, equal_var=False)
    assert_array_almost_equal([t,p], (tr_uneq_n, pr_uneq_n))
    assert_array_almost_equal(stats.ttest_ind_from_stats(*_desc_stats(rvs1,
                                                                      rvs3),
                                                         equal_var=False),
                              (t, p))

    t,p = stats.ttest_ind(rvs1_2D.T, rvs2_2D.T, axis=0, equal_var=False)
    assert_array_almost_equal([t,p],tpr)
    args = _desc_stats(rvs1_2D.T, rvs2_2D.T)
    assert_array_almost_equal(stats.ttest_ind_from_stats(*args,
                                                         equal_var=False),
                              (t, p))

    t,p = stats.ttest_ind(rvs1_2D, rvs2_2D, axis=1, equal_var=False)
    assert_array_almost_equal([t,p],tpr)
    args = _desc_stats(rvs1_2D, rvs2_2D, axis=1)
    assert_array_almost_equal(stats.ttest_ind_from_stats(*args,
                                                         equal_var=False),
                              (t, p))

    # test for namedtuple attribute results
    attributes = ('statistic', 'pvalue')
    res = stats.ttest_ind(rvs1, rvs2, axis=0, equal_var=False)
    check_named_results(res, attributes)

    # test on 3 dimensions
    rvs1_3D = np.dstack([rvs1_2D,rvs1_2D,rvs1_2D])
    rvs2_3D = np.dstack([rvs2_2D,rvs2_2D,rvs2_2D])
    t,p = stats.ttest_ind(rvs1_3D, rvs2_3D, axis=1, equal_var=False)
    assert_almost_equal(np.abs(t), np.abs(tr))
    assert_array_almost_equal(np.abs(p), pr)
    assert_equal(t.shape, (2, 3))
    args = _desc_stats(rvs1_3D, rvs2_3D, axis=1)
    t, p = stats.ttest_ind_from_stats(*args, equal_var=False)
    assert_almost_equal(np.abs(t), np.abs(tr))
    assert_array_almost_equal(np.abs(p), pr)
    assert_equal(t.shape, (2, 3))

    t,p = stats.ttest_ind(np.rollaxis(rvs1_3D,2), np.rollaxis(rvs2_3D,2),
                                   axis=2, equal_var=False)
    assert_array_almost_equal(np.abs(t), np.abs(tr))
    assert_array_almost_equal(np.abs(p), pr)
    assert_equal(t.shape, (3, 2))
    args = _desc_stats(np.rollaxis(rvs1_3D, 2),
                       np.rollaxis(rvs2_3D, 2), axis=2)
    t, p = stats.ttest_ind_from_stats(*args, equal_var=False)
    assert_array_almost_equal(np.abs(t), np.abs(tr))
    assert_array_almost_equal(np.abs(p), pr)
    assert_equal(t.shape, (3, 2))

    olderr = np.seterr(all='ignore')
    try:
        # test zero division problem
        t,p = stats.ttest_ind([0,0,0],[1,1,1], equal_var=False)
        assert_equal((np.abs(t),p), (np.inf, 0))
        assert_equal(stats.ttest_ind([0,0,0], [0,0,0], equal_var=False), (np.nan, np.nan))

        # check that nan in input array result in nan output
        anan = np.array([[1,np.nan],[-1,1]])
        assert_equal(stats.ttest_ind(anan, np.zeros((2,2)), equal_var=False),
                     ([0, np.nan], [1,np.nan]))
    finally:
        np.seterr(**olderr)


def test_ttest_1samp_new():
    n1, n2, n3 = (10,15,20)
    rvn1 = stats.norm.rvs(loc=5,scale=10,size=(n1,n2,n3))

    # check multidimensional array and correct axis handling
    # deterministic rvn1 and rvn2 would be better as in test_ttest_rel
    t1,p1 = stats.ttest_1samp(rvn1[:,:,:], np.ones((n2,n3)),axis=0)
    t2,p2 = stats.ttest_1samp(rvn1[:,:,:], 1,axis=0)
    t3,p3 = stats.ttest_1samp(rvn1[:,0,0], 1)
    assert_array_almost_equal(t1,t2, decimal=14)
    assert_almost_equal(t1[0,0],t3, decimal=14)
    assert_equal(t1.shape, (n2,n3))

    t1,p1 = stats.ttest_1samp(rvn1[:,:,:], np.ones((n1,n3)),axis=1)
    t2,p2 = stats.ttest_1samp(rvn1[:,:,:], 1,axis=1)
    t3,p3 = stats.ttest_1samp(rvn1[0,:,0], 1)
    assert_array_almost_equal(t1,t2, decimal=14)
    assert_almost_equal(t1[0,0],t3, decimal=14)
    assert_equal(t1.shape, (n1,n3))

    t1,p1 = stats.ttest_1samp(rvn1[:,:,:], np.ones((n1,n2)),axis=2)
    t2,p2 = stats.ttest_1samp(rvn1[:,:,:], 1,axis=2)
    t3,p3 = stats.ttest_1samp(rvn1[0,0,:], 1)
    assert_array_almost_equal(t1,t2, decimal=14)
    assert_almost_equal(t1[0,0],t3, decimal=14)
    assert_equal(t1.shape, (n1,n2))

    olderr = np.seterr(all='ignore')
    try:
        # test zero division problem
        t,p = stats.ttest_1samp([0,0,0], 1)
        assert_equal((np.abs(t),p), (np.inf, 0))
        assert_equal(stats.ttest_1samp([0,0,0], 0), (np.nan, np.nan))

        # check that nan in input array result in nan output
        anan = np.array([[1,np.nan],[-1,1]])
        assert_equal(stats.ttest_1samp(anan, 0),([0, np.nan], [1,np.nan]))
    finally:
        np.seterr(**olderr)


class TestDescribe(TestCase):
    def test_describe_scalar(self):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=RuntimeWarning)
            n, mm, m, v, sk, kurt = stats.describe(4.)
        assert_equal(n, 1)
        assert_equal(mm, (4.0, 4.0))
        assert_equal(m, 4.0)
        assert_(np.isnan(v))
        # not sure about precision with sk, skc
        assert_array_almost_equal(sk, 0.0, decimal=13)
        assert_array_almost_equal(kurt, -3.0, decimal=13)

    def test_describe_numbers(self):
        x = np.vstack((np.ones((3,4)), 2 * np.ones((2,4))))
        nc, mmc = (5, ([1., 1., 1., 1.], [2., 2., 2., 2.]))
        mc = np.array([1.4, 1.4, 1.4, 1.4])
        vc = np.array([0.3, 0.3, 0.3, 0.3])
        skc = [0.40824829046386357] * 4
        kurtc = [-1.833333333333333] * 4
        n, mm, m, v, sk, kurt = stats.describe(x)
        assert_equal(n, nc)
        assert_equal(mm, mmc)
        assert_equal(m, mc)
        assert_equal(v, vc)
        # not sure about precision with sk, skc
        assert_array_almost_equal(sk, skc, decimal=13)
        assert_array_almost_equal(kurt, kurtc, decimal=13)
        n, mm, m, v, sk, kurt = stats.describe(x.T, axis=1)
        assert_equal(n, nc)
        assert_equal(mm, mmc)
        assert_equal(m, mc)
        assert_equal(v, vc)
        # not sure about precision with sk, skc
        assert_array_almost_equal(sk, skc, decimal=13)
        assert_array_almost_equal(kurt, kurtc, decimal=13)

    def test_describe_result_attributes(self):
        actual = stats.describe(np.arange(5))
        attributes = ('nobs', 'minmax', 'mean', 'variance', 'skewness',
                      'kurtosis')
        check_named_results(actual, attributes)

def test_normalitytests():
    yield (assert_raises, ValueError, stats.skewtest, 4.)
    yield (assert_raises, ValueError, stats.kurtosistest, 4.)
    yield (assert_raises, ValueError, stats.normaltest, 4.)

    # numbers verified with R: dagoTest in package fBasics
    st_normal, st_skew, st_kurt = (3.92371918, 1.98078826, -0.01403734)
    pv_normal, pv_skew, pv_kurt = (0.14059673, 0.04761502, 0.98880019)
    x = np.array((-2,-1,0,1,2,3)*4)**2
    attributes = ('statistic', 'pvalue')

    yield assert_array_almost_equal, stats.normaltest(x), (st_normal, pv_normal)
    check_named_results(stats.normaltest(x), attributes)
    yield assert_array_almost_equal, stats.skewtest(x), (st_skew, pv_skew)
    check_named_results(stats.skewtest(x), attributes)
    yield assert_array_almost_equal, stats.kurtosistest(x), (st_kurt, pv_kurt)
    check_named_results(stats.kurtosistest(x), attributes)

    # Test axis=None (equal to axis=0 for 1-D input)
    yield (assert_array_almost_equal, stats.normaltest(x, axis=None),
           (st_normal, pv_normal))
    yield (assert_array_almost_equal, stats.skewtest(x, axis=None),
           (st_skew, pv_skew))
    yield (assert_array_almost_equal, stats.kurtosistest(x, axis=None),
           (st_kurt, pv_kurt))


class TestRankSums(TestCase):
    def test_ranksums_result_attributes(self):
        res = stats.ranksums(np.arange(5), np.arange(25))
        attributes = ('statistic', 'pvalue')
        check_named_results(res, attributes)


class TestJarqueBera(TestCase):
    def test_jarque_bera_stats(self):
        np.random.seed(987654321)
        x = np.random.normal(0, 1, 100000)
        y = np.random.chisquare(10000, 100000)
        z = np.random.rayleigh(1, 100000)

        assert_(stats.jarque_bera(x)[1] > stats.jarque_bera(y)[1])
        assert_(stats.jarque_bera(x)[1] > stats.jarque_bera(z)[1])
        assert_(stats.jarque_bera(y)[1] > stats.jarque_bera(z)[1])

    def test_jarque_bera_array_like(self):
        np.random.seed(987654321)
        x = np.random.normal(0, 1, 100000)

        JB1, p1 = stats.jarque_bera(list(x))
        JB2, p2 = stats.jarque_bera(tuple(x))
        JB3, p3 = stats.jarque_bera(x.reshape(2, 50000))

        assert_(JB1 == JB2 == JB3)
        assert_(p1 == p2 == p3)

    def test_jarque_bera_size(self):
        assert_raises(ValueError, stats.jarque_bera, [])


def test_skewtest_too_few_samples():
    # Regression test for ticket #1492.
    # skewtest requires at least 8 samples; 7 should raise a ValueError.
    x = np.arange(7.0)
    assert_raises(ValueError, stats.skewtest, x)


def test_kurtosistest_too_few_samples():
    # Regression test for ticket #1425.
    # kurtosistest requires at least 5 samples; 4 should raise a ValueError.
    x = np.arange(4.0)
    assert_raises(ValueError, stats.kurtosistest, x)


def test_mannwhitneyu():
    x = np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 2., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 2., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 2., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 2., 1., 1., 1., 1., 2., 1., 1., 2., 1., 1.,
        2., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 2., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        2., 1., 1., 1., 1., 1., 1., 1., 1., 1., 2., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 3., 1., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 1.])

    y = np.array([1., 1., 1., 1., 1., 1., 1., 2., 1., 2., 1., 1., 1.,
        1., 2., 1., 1., 1., 2., 1., 1., 1., 1., 1., 2., 1., 1., 3., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 2., 1., 2., 1., 1., 1., 1.,
        1., 1., 2., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 2., 1., 1., 1., 1., 1., 2., 2., 1., 1., 2., 1., 1., 2.,
        1., 2., 1., 1., 1., 1., 2., 2., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 2., 1., 1., 1., 1., 1., 2., 2., 2., 1.,
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
        1., 2., 1., 1., 2., 1., 1., 1., 1., 2., 1., 1., 1., 1., 1., 1.,
        1., 1., 1., 1., 1., 1., 2., 1., 1., 1., 2., 1., 1., 1., 1., 1.,
        1.])
    # p-value verified with matlab and R to 5 significant digits
    assert_array_almost_equal(stats.stats.mannwhitneyu(x,y),
                    (16980.5, 2.8214327656317373e-005), decimal=12)

    # test for namedtuple attribute results
    attributes = ('statistic', 'pvalue')
    res = stats.mannwhitneyu(x, y)
    check_named_results(res, attributes)


def test_pointbiserial():
    # same as mstats test except for the nan
    # Test data: http://support.sas.com/ctx/samples/index.jsp?sid=490&tab=output
    x = [1,0,1,1,1,1,0,1,0,0,0,1,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,
         0,0,0,0,1]
    y = [14.8,13.8,12.4,10.1,7.1,6.1,5.8,4.6,4.3,3.5,3.3,3.2,3.0,
         2.8,2.8,2.5,2.4,2.3,2.1,1.7,1.7,1.5,1.3,1.3,1.2,1.2,1.1,
         0.8,0.7,0.6,0.5,0.2,0.2,0.1]
    assert_almost_equal(stats.pointbiserialr(x, y)[0], 0.36149, 5)

    # test for namedtuple attribute results
    attributes = ('correlation', 'pvalue')
    res = stats.pointbiserialr(x, y)
    check_named_results(res, attributes)


def test_obrientransform():
    # A couple tests calculated by hand.
    x1 = np.array([0, 2, 4])
    t1 = stats.obrientransform(x1)
    expected = [7, -2, 7]
    assert_allclose(t1[0], expected)

    x2 = np.array([0, 3, 6, 9])
    t2 = stats.obrientransform(x2)
    expected = np.array([30, 0, 0, 30])
    assert_allclose(t2[0], expected)

    # Test two arguments.
    a, b = stats.obrientransform(x1, x2)
    assert_equal(a, t1[0])
    assert_equal(b, t2[0])

    # Test three arguments.
    a, b, c = stats.obrientransform(x1, x2, x1)
    assert_equal(a, t1[0])
    assert_equal(b, t2[0])
    assert_equal(c, t1[0])

    # This is a regression test to check np.var replacement.
    # The author of this test didn't separately verify the numbers.
    x1 = np.arange(5)
    result = np.array(
      [[5.41666667, 1.04166667, -0.41666667, 1.04166667, 5.41666667],
       [21.66666667, 4.16666667, -1.66666667, 4.16666667, 21.66666667]])
    assert_array_almost_equal(stats.obrientransform(x1, 2*x1), result, decimal=8)

    # Example from "O'Brien Test for Homogeneity of Variance"
    # by Herve Abdi.
    values = range(5, 11)
    reps = np.array([5, 11, 9, 3, 2, 2])
    data = np.repeat(values, reps)
    transformed_values = np.array([3.1828, 0.5591, 0.0344,
                                   1.6086, 5.2817, 11.0538])
    expected = np.repeat(transformed_values, reps)
    result = stats.obrientransform(data)
    assert_array_almost_equal(result[0], expected, decimal=4)


class HarMeanTestCase:
    def test_1dlist(self):
        #  Test a 1d list
        a = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        b = 34.1417152147
        self.do(a, b)

    def test_1darray(self):
        #  Test a 1d array
        a = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
        b = 34.1417152147
        self.do(a, b)

    def test_1dma(self):
        #  Test a 1d masked array
        a = np.ma.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
        b = 34.1417152147
        self.do(a, b)

    def test_1dmavalue(self):
        #  Test a 1d masked array with a masked value
        a = np.ma.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
                      mask=[0,0,0,0,0,0,0,0,0,1])
        b = 31.8137186141
        self.do(a, b)

    # Note the next tests use axis=None as default, not axis=0
    def test_2dlist(self):
        #  Test a 2d list
        a = [[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = 38.6696271841
        self.do(a, b)

    def test_2darray(self):
        #  Test a 2d array
        a = [[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = 38.6696271841
        self.do(np.array(a), b)

    def test_2dma(self):
        #  Test a 2d masked array
        a = [[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = 38.6696271841
        self.do(np.ma.array(a), b)

    def test_2daxis0(self):
        #  Test a 2d list with axis=0
        a = [[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = np.array([22.88135593, 39.13043478, 52.90076336, 65.45454545])
        self.do(a, b, axis=0)

    def test_2daxis1(self):
        #  Test a 2d list with axis=1
        a = [[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = np.array([19.2, 63.03939962, 103.80078637])
        self.do(a, b, axis=1)

    def test_2dmatrixdaxis0(self):
        #  Test a 2d list with axis=0
        a = [[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = np.matrix([[22.88135593, 39.13043478, 52.90076336, 65.45454545]])
        self.do(np.matrix(a), b, axis=0)

    def test_2dmatrixaxis1(self):
        #  Test a 2d list with axis=1
        a = [[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = np.matrix([[19.2, 63.03939962, 103.80078637]]).T
        self.do(np.matrix(a), b, axis=1)


class TestHarMean(HarMeanTestCase, TestCase):
    def do(self, a, b, axis=None, dtype=None):
        x = stats.hmean(a, axis=axis, dtype=dtype)
        assert_almost_equal(b, x)
        assert_equal(x.dtype, dtype)


class GeoMeanTestCase:
    def test_1dlist(self):
        #  Test a 1d list
        a = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        b = 45.2872868812
        self.do(a, b)

    def test_1darray(self):
        #  Test a 1d array
        a = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
        b = 45.2872868812
        self.do(a, b)

    def test_1dma(self):
        #  Test a 1d masked array
        a = np.ma.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
        b = 45.2872868812
        self.do(a, b)

    def test_1dmavalue(self):
        #  Test a 1d masked array with a masked value
        a = np.ma.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100], mask=[0,0,0,0,0,0,0,0,0,1])
        b = 41.4716627439
        self.do(a, b)

    # Note the next tests use axis=None as default, not axis=0
    def test_2dlist(self):
        #  Test a 2d list
        a = [[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = 52.8885199
        self.do(a, b)

    def test_2darray(self):
        #  Test a 2d array
        a = [[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = 52.8885199
        self.do(np.array(a), b)

    def test_2dma(self):
        #  Test a 2d masked array
        a = [[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = 52.8885199
        self.do(np.ma.array(a), b)

    def test_2daxis0(self):
        #  Test a 2d list with axis=0
        a = [[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = np.array([35.56893304, 49.32424149, 61.3579244, 72.68482371])
        self.do(a, b, axis=0)

    def test_2daxis1(self):
        #  Test a 2d list with axis=1
        a = [[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = np.array([22.13363839, 64.02171746, 104.40086817])
        self.do(a, b, axis=1)

    def test_2dmatrixdaxis0(self):
        #  Test a 2d list with axis=0
        a = [[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = np.matrix([[35.56893304, 49.32424149, 61.3579244, 72.68482371]])
        self.do(np.matrix(a), b, axis=0)

    def test_2dmatrixaxis1(self):
        #  Test a 2d list with axis=1
        a = [[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = np.matrix([[22.13363839, 64.02171746, 104.40086817]]).T
        self.do(np.matrix(a), b, axis=1)

    def test_1dlist0(self):
        #  Test a 1d list with zero element
        a = [10, 20, 30, 40, 50, 60, 70, 80, 90, 0]
        b = 0.0  # due to exp(-inf)=0
        olderr = np.seterr(all='ignore')
        try:
            self.do(a, b)
        finally:
            np.seterr(**olderr)

    def test_1darray0(self):
        #  Test a 1d array with zero element
        a = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 0])
        b = 0.0  # due to exp(-inf)=0
        olderr = np.seterr(all='ignore')
        try:
            self.do(a, b)
        finally:
            np.seterr(**olderr)

    def test_1dma0(self):
        #  Test a 1d masked array with zero element
        a = np.ma.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 0])
        b = 41.4716627439
        olderr = np.seterr(all='ignore')
        try:
            self.do(a, b)
        finally:
            np.seterr(**olderr)

    def test_1dmainf(self):
        #  Test a 1d masked array with negative element
        a = np.ma.array([10, 20, 30, 40, 50, 60, 70, 80, 90, -1])
        b = 41.4716627439
        olderr = np.seterr(all='ignore')
        try:
            self.do(a, b)
        finally:
            np.seterr(**olderr)


class TestGeoMean(GeoMeanTestCase, TestCase):
    def do(self, a, b, axis=None, dtype=None):
        # Note this doesn't test when axis is not specified
        x = stats.gmean(a, axis=axis, dtype=dtype)
        assert_almost_equal(b, x)
        assert_equal(x.dtype, dtype)


def test_binomtest():
    # precision tests compared to R for ticket:986
    pp = np.concatenate((np.linspace(0.1,0.2,5), np.linspace(0.45,0.65,5),
                          np.linspace(0.85,0.95,5)))
    n = 501
    x = 450
    results = [0.0, 0.0, 1.0159969301994141e-304,
    2.9752418572150531e-275, 7.7668382922535275e-250,
    2.3381250925167094e-099, 7.8284591587323951e-081,
    9.9155947819961383e-065, 2.8729390725176308e-050,
    1.7175066298388421e-037, 0.0021070691951093692,
    0.12044570587262322, 0.88154763174802508, 0.027120993063129286,
    2.6102587134694721e-006]

    for p, res in zip(pp,results):
        assert_approx_equal(stats.binom_test(x, n, p), res,
                            significant=12, err_msg='fail forp=%f' % p)

    assert_approx_equal(stats.binom_test(50,100,0.1), 5.8320387857343647e-024,
                            significant=12, err_msg='fail forp=%f' % p)


def test_binomtest2():
    # test added for issue #2384
    res2 = [
    [1.0, 1.0],
    [0.5,1.0,0.5],
    [0.25,1.00,1.00,0.25],
    [0.125,0.625,1.000,0.625,0.125],
    [0.0625,0.3750,1.0000,1.0000,0.3750,0.0625],
    [0.03125,0.21875,0.68750,1.00000,0.68750,0.21875,0.03125],
    [0.015625,0.125000,0.453125,1.000000,1.000000,0.453125,0.125000,0.015625],
    [0.0078125,0.0703125,0.2890625,0.7265625,1.0000000,0.7265625,0.2890625,
     0.0703125,0.0078125],
    [0.00390625,0.03906250,0.17968750,0.50781250,1.00000000,1.00000000,
     0.50781250,0.17968750,0.03906250,0.00390625],
    [0.001953125,0.021484375,0.109375000,0.343750000,0.753906250,1.000000000,
     0.753906250,0.343750000,0.109375000,0.021484375,0.001953125]
    ]

    for k in range(1, 11):
        res1 = [stats.binom_test(v, k, 0.5) for v in range(k + 1)]
        assert_almost_equal(res1, res2[k-1], decimal=10)


def test_binomtest3():
    # test added for issue #2384
    # test when x == n*p and neighbors
    res3 = [stats.binom_test(v, v*k, 1./k) for v in range(1, 11)
                                           for k in range(2, 11)]
    assert_equal(res3, np.ones(len(res3), int))

    #> bt=c()
    #> for(i in as.single(1:10)){for(k in as.single(2:10)){bt = c(bt, binom.test(i-1, k*i,(1/k))$p.value); print(c(i+1, k*i,(1/k)))}}
    binom_testm1 = np.array([
         0.5, 0.5555555555555556, 0.578125, 0.5904000000000003,
         0.5981224279835393, 0.603430543396034, 0.607304096221924,
         0.610255656871054, 0.612579511000001, 0.625, 0.670781893004115,
         0.68853759765625, 0.6980101120000006, 0.703906431368616,
         0.70793209416498, 0.7108561134173507, 0.713076544331419,
         0.714820192935702, 0.6875, 0.7268709038256367, 0.7418963909149174,
         0.74986110468096, 0.7548015520398076, 0.7581671424768577,
         0.760607984787832, 0.762459425024199, 0.7639120677676575, 0.7265625,
         0.761553963657302, 0.774800934828818, 0.7818005980538996,
         0.78613491480358, 0.789084353140195, 0.7912217659828884,
         0.79284214559524, 0.794112956558801, 0.75390625, 0.7856929451142176,
         0.7976688481430754, 0.8039848974727624, 0.807891868948366,
         0.8105487660137676, 0.812473307174702, 0.8139318233591120,
         0.815075399104785, 0.7744140625, 0.8037322594985427,
         0.814742863657656, 0.8205425178645808, 0.8241275984172285,
         0.8265645374416, 0.8283292196088257, 0.829666291102775,
         0.8307144686362666, 0.7905273437499996, 0.8178712053954738,
         0.828116983756619, 0.833508948940494, 0.8368403871552892,
         0.839104213210105, 0.840743186196171, 0.84198481438049,
         0.8429580531563676, 0.803619384765625, 0.829338573944648,
         0.8389591907548646, 0.84401876783902, 0.84714369697889,
         0.8492667010581667, 0.850803474598719, 0.851967542858308,
         0.8528799045949524, 0.8145294189453126, 0.838881732845347,
         0.847979024541911, 0.852760894015685, 0.8557134656773457,
         0.8577190131799202, 0.85917058278431, 0.860270010472127,
         0.861131648404582, 0.823802947998047, 0.846984756807511,
         0.855635653643743, 0.860180994825685, 0.86298688573253,
         0.864892525675245, 0.866271647085603, 0.867316125625004,
         0.8681346531755114
        ])

    # > bt=c()
    # > for(i in as.single(1:10)){for(k in as.single(2:10)){bt = c(bt, binom.test(i+1, k*i,(1/k))$p.value); print(c(i+1, k*i,(1/k)))}}

    binom_testp1 = np.array([
         0.5, 0.259259259259259, 0.26171875, 0.26272, 0.2632244513031551,
         0.2635138663069203, 0.2636951804161073, 0.2638162407564354,
         0.2639010709000002, 0.625, 0.4074074074074074, 0.42156982421875,
         0.4295746560000003, 0.43473045988554, 0.4383309503172684,
         0.4409884859402103, 0.4430309389962837, 0.444649849401104, 0.6875,
         0.4927602499618962, 0.5096031427383425, 0.5189636628480,
         0.5249280070771274, 0.5290623300865124, 0.5320974248125793,
         0.5344204730474308, 0.536255847400756, 0.7265625, 0.5496019313526808,
         0.5669248746708034, 0.576436455045805, 0.5824538812831795,
         0.5866053321547824, 0.589642781414643, 0.5919618019300193,
         0.593790427805202, 0.75390625, 0.590868349763505, 0.607983393277209,
         0.617303847446822, 0.623172512167948, 0.627208862156123,
         0.6301556891501057, 0.632401894928977, 0.6341708982290303,
         0.7744140625, 0.622562037497196, 0.639236102912278, 0.648263335014579,
         0.65392850011132, 0.657816519817211, 0.660650782947676,
         0.662808780346311, 0.6645068560246006, 0.7905273437499996,
         0.6478843304312477, 0.6640468318879372, 0.6727589686071775,
         0.6782129857784873, 0.681950188903695, 0.684671508668418,
         0.686741824999918, 0.688369886732168, 0.803619384765625,
         0.668716055304315, 0.684360013879534, 0.6927642396829181,
         0.6980155964704895, 0.701609591890657, 0.7042244320992127,
         0.7062125081341817, 0.707775152962577, 0.8145294189453126,
         0.686243374488305, 0.7013873696358975, 0.709501223328243,
         0.714563595144314, 0.718024953392931, 0.7205416252126137,
         0.722454130389843, 0.723956813292035, 0.823802947998047,
         0.701255953767043, 0.715928221686075, 0.723772209289768,
         0.7286603031173616, 0.7319999279787631, 0.7344267920995765,
         0.736270323773157, 0.737718376096348
        ])

    res4_p1 = [stats.binom_test(v+1, v*k, 1./k) for v in range(1, 11)
                                                for k in range(2, 11)]
    res4_m1 = [stats.binom_test(v-1, v*k, 1./k) for v in range(1, 11)
                                                for k in range(2, 11)]

    assert_almost_equal(res4_p1, binom_testp1, decimal=13)
    assert_almost_equal(res4_m1, binom_testm1, decimal=13)


class TestTrim(object):
    # test trim functions
    def test_trim1(self):
        a = np.arange(11)
        assert_equal(stats.trim1(a, 0.1), np.arange(10))
        assert_equal(stats.trim1(a, 0.2), np.arange(9))
        assert_equal(stats.trim1(a, 0.2, tail='left'), np.arange(2,11))
        assert_equal(stats.trim1(a, 3/11., tail='left'), np.arange(3,11))

    def test_trimboth(self):
        a = np.arange(11)
        assert_equal(stats.trimboth(a, 3/11.), np.arange(3,8))
        assert_equal(stats.trimboth(a, 0.2), np.array([2, 3, 4, 5, 6, 7, 8]))
        assert_equal(stats.trimboth(np.arange(24).reshape(6,4), 0.2),
                     np.arange(4,20).reshape(4,4))
        assert_equal(stats.trimboth(np.arange(24).reshape(4,6).T, 2/6.),
               np.array([[2, 8, 14, 20],[3, 9, 15, 21]]))
        assert_raises(ValueError, stats.trimboth,
               np.arange(24).reshape(4,6).T, 4/6.)

    def test_trim_mean(self):
        # don't use pre-sorted arrays
        a = np.array([4, 8, 2, 0, 9, 5, 10, 1, 7, 3, 6])
        idx = np.array([3, 5, 0, 1, 2, 4])
        a2 = np.arange(24).reshape(6, 4)[idx, :]
        a3 = np.arange(24).reshape(6, 4, order='F')[idx, :]
        assert_equal(stats.trim_mean(a3, 2/6.),
                     np.array([2.5, 8.5, 14.5, 20.5]))
        assert_equal(stats.trim_mean(a2, 2/6.),
                     np.array([10., 11., 12., 13.]))
        idx4 = np.array([1, 0, 3, 2])
        a4 = np.arange(24).reshape(4, 6)[idx4, :]
        assert_equal(stats.trim_mean(a4, 2/6.),
                     np.array([9., 10., 11., 12., 13., 14.]))
        # shuffled arange(24) as array_like
        a = [7, 11, 12, 21, 16, 6, 22, 1, 5, 0, 18, 10, 17, 9, 19, 15, 23,
             20, 2, 14, 4, 13, 8, 3]
        assert_equal(stats.trim_mean(a, 2/6.), 11.5)
        assert_equal(stats.trim_mean([5,4,3,1,2,0], 2/6.), 2.5)

        # check axis argument
        np.random.seed(1234)
        a = np.random.randint(20, size=(5, 6, 4, 7))
        for axis in [0, 1, 2, 3, -1]:
            res1 = stats.trim_mean(a, 2/6., axis=axis)
            res2 = stats.trim_mean(np.rollaxis(a, axis), 2/6.)
            assert_equal(res1, res2)

        res1 = stats.trim_mean(a, 2/6., axis=None)
        res2 = stats.trim_mean(a.ravel(), 2/6.)
        assert_equal(res1, res2)

        assert_raises(ValueError, stats.trim_mean, a, 0.6)


class TestSigamClip(object):
    def test_sigmaclip1(self):
        a = np.concatenate((np.linspace(9.5,10.5,31),np.linspace(0,20,5)))
        fact = 4  # default
        c, low, upp = stats.sigmaclip(a)
        assert_(c.min() > low)
        assert_(c.max() < upp)
        assert_equal(low, c.mean() - fact*c.std())
        assert_equal(upp, c.mean() + fact*c.std())
        assert_equal(c.size, a.size)

    def test_sigmaclip2(self):
        a = np.concatenate((np.linspace(9.5,10.5,31),np.linspace(0,20,5)))
        fact = 1.5
        c, low, upp = stats.sigmaclip(a, fact, fact)
        assert_(c.min() > low)
        assert_(c.max() < upp)
        assert_equal(low, c.mean() - fact*c.std())
        assert_equal(upp, c.mean() + fact*c.std())
        assert_equal(c.size, 4)
        assert_equal(a.size, 36)  # check original array unchanged

    def test_sigmaclip3(self):
        a = np.concatenate((np.linspace(9.5,10.5,11),np.linspace(-100,-50,3)))
        fact = 1.8
        c, low, upp = stats.sigmaclip(a, fact, fact)
        assert_(c.min() > low)
        assert_(c.max() < upp)
        assert_equal(low, c.mean() - fact*c.std())
        assert_equal(upp, c.mean() + fact*c.std())
        assert_equal(c, np.linspace(9.5,10.5,11))

    def test_sigmaclip_result_attributes(self):
        a = np.concatenate((np.linspace(9.5, 10.5, 11),
                            np.linspace(-100, -50, 3)))
        fact = 1.8
        res = stats.sigmaclip(a, fact, fact)
        attributes = ('clipped', 'lower', 'upper')
        check_named_results(res, attributes)

class TestFOneWay(TestCase):
    def test_trivial(self):
        # A trivial test of stats.f_oneway, with F=0.
        F, p = stats.f_oneway([0,2], [0,2])
        assert_equal(F, 0.0)

    def test_basic(self):
        # Despite being a floating point calculation, this data should
        # result in F being exactly 2.0.
        F, p = stats.f_oneway([0,2], [2,4])
        assert_equal(F, 2.0)

    def test_large_integer_array(self):
        a = np.array([655, 788], dtype=np.uint16)
        b = np.array([789, 772], dtype=np.uint16)
        F, p = stats.f_oneway(a, b)
        assert_almost_equal(F, 0.77450216931805538)

    def test_result_attributes(self):
        a = np.array([655, 788], dtype=np.uint16)
        b = np.array([789, 772], dtype=np.uint16)
        res = stats.f_oneway(a, b)
        attributes = ('statistic', 'pvalue')
        check_named_results(res, attributes)

class TestKruskal(TestCase):
    def test_simple(self):
        x = [1]
        y = [2]
        h, p = stats.kruskal(x, y)
        assert_equal(h, 1.0)
        assert_approx_equal(p, stats.chisqprob(h, 1))
        h, p = stats.kruskal(np.array(x), np.array(y))
        assert_equal(h, 1.0)
        assert_approx_equal(p, stats.chisqprob(h, 1))

    def test_basic(self):
        x = [1, 3, 5, 7, 9]
        y = [2, 4, 6, 8, 10]
        h, p = stats.kruskal(x, y)
        assert_approx_equal(h, 3./11, significant=10)
        assert_approx_equal(p, stats.chisqprob(3./11, 1))
        h, p = stats.kruskal(np.array(x), np.array(y))
        assert_approx_equal(h, 3./11, significant=10)
        assert_approx_equal(p, stats.chisqprob(3./11, 1))

    def test_simple_tie(self):
        x = [1]
        y = [1, 2]
        h_uncorr = 1.5**2 + 2*2.25**2 - 12
        corr = 0.75
        expected = h_uncorr / corr   # 0.5
        h, p = stats.kruskal(x, y)
        # Since the expression is simple and the exact answer is 0.5, it
        # should be safe to use assert_equal().
        assert_equal(h, expected)

    def test_another_tie(self):
        x = [1, 1, 1, 2]
        y = [2, 2, 2, 2]
        h_uncorr = (12. / 8. / 9.) * 4 * (3**2 + 6**2) - 3 * 9
        corr = 1 - float(3**3 - 3 + 5**3 - 5) / (8**3 - 8)
        expected = h_uncorr / corr
        h, p = stats.kruskal(x, y)
        assert_approx_equal(h, expected)

    def test_three_groups(self):
        # A test of stats.kruskal with three groups, with ties.
        x = [1, 1, 1]
        y = [2, 2, 2]
        z = [2, 2]
        h_uncorr = (12. / 8. / 9.) * (3*2**2 + 3*6**2 + 2*6**2) - 3 * 9  # 5.0
        corr = 1 - float(3**3 - 3 + 5**3 - 5) / (8**3 - 8)
        expected = h_uncorr / corr  # 7.0
        h, p = stats.kruskal(x, y, z)
        assert_approx_equal(h, expected)
        assert_approx_equal(p, stats.chisqprob(h, 2))

    def test_kruskal_result_attributes(self):
        x = [1, 3, 5, 7, 9]
        y = [2, 4, 6, 8, 10]
        res = stats.kruskal(x, y)
        attributes = ('statistic', 'pvalue')
        check_named_results(res, attributes)


class TestCombinePvalues(TestCase):

    def test_fisher(self):
        # Example taken from http://en.wikipedia.org/wiki/Fisher's_exact_test#Example
        xsq, p = stats.combine_pvalues([.01, .2, .3], method='fisher')
        assert_approx_equal(p, 0.02156, significant=4)

    def test_stouffer(self):
        Z, p = stats.combine_pvalues([.01, .2, .3], method='stouffer')
        assert_approx_equal(p, 0.01651, significant=4)

    def test_stouffer2(self):
        Z, p = stats.combine_pvalues([.5, .5, .5], method='stouffer')
        assert_approx_equal(p, 0.5, significant=4)

    def test_weighted_stouffer(self):
        Z, p = stats.combine_pvalues([.01, .2, .3], method='stouffer',
                                     weights=np.ones(3))
        assert_approx_equal(p, 0.01651, significant=4)

    def test_weighted_stouffer2(self):
        Z, p = stats.combine_pvalues([.01, .2, .3], method='stouffer',
                                     weights=np.array((1, 4, 9)))
        assert_approx_equal(p, 0.1464, significant=4)

if __name__ == "__main__":
    run_module_suite()
