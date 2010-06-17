""" Test functions for stats module

    WRITTEN BY LOUIS LUANGKESORN <lluang@yahoo.com> FOR THE STATS MODULE
    BASED ON WILKINSON'S STATISTICS QUIZ
    http://www.stanford.edu/~clint/bench/wilk.txt

"""

from numpy.testing import *
from numpy import array, arange, zeros, ravel, float32, float64, power
import numpy as np

import scipy.stats as stats


""" Numbers in docstrings begining with 'W' refer to the section numbers
    and headings found in the STATISTICS QUIZ of Leland Wilkinson.  These are
    considered to be essential functionality.  True testing and
    evaluation of a statistics package requires use of the
    NIST Statistical test data.  See McCoullough(1999) Assessing The Reliability
    of Statistical Software for a test methodology and its
    implementation in testing SAS, SPSS, and S-Plus
"""

##  Datasets
##  These data sets are from the nasty.dat sets used by Wilkinson
##  for MISS, need to be able to represent missing values
##  For completeness, I should write the relevant tests and count them as failures
##  Somewhat acceptable, since this is still beta software.  It would count as a
##  good target for 1.0 status
X = array([1,2,3,4,5,6,7,8,9],float)
ZERO= array([0,0,0,0,0,0,0,0,0], float)
#MISS=array([.,.,.,.,.,.,.,.,.], float)
BIG=array([99999991,99999992,99999993,99999994,99999995,99999996,99999997,99999998,99999999],float)
LITTLE=array([0.99999991,0.99999992,0.99999993,0.99999994,0.99999995,0.99999996,0.99999997,0.99999998,0.99999999],float)
HUGE=array([1e+12,2e+12,3e+12,4e+12,5e+12,6e+12,7e+12,8e+12,9e+12],float)
TINY=array([1e-12,2e-12,3e-12,4e-12,5e-12,6e-12,7e-12,8e-12,9e-12],float)
ROUND=array([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5],float)
X2 = X * X
X3 = X2 * X
X4 = X3 * X
X5 = X4 * X
X6 = X5 * X
X7 = X6 * X
X8 = X7 * X
X9 = X8 * X

class TestRound(TestCase):
    """ W.II. ROUND

        You should get the numbers 1 to 9.  Many language compilers,
        such as Turbo Pascal and Lattice C, fail this test (they round
        numbers inconsistently). Needless to say, statical packages
        written in these languages may fail the test as well.  You can
        also check the following expressions:
            Y = INT(2.6*7 -0.2)                   (Y should be 18)
            Y = 2-INT(EXP(LOG(SQR(2)*SQR(2))))    (Y should be 0)
            Y = INT(3-EXP(LOG(SQR(2)*SQR(2))))    (Y should be 1)
        INT is the integer function.  It converts decimal numbers to
        integers by throwing away numbers after the decimal point.  EXP
        is exponential, LOG is logarithm, and SQR is suqare root.  You may
        have to substitute similar names for these functions for different
        packages.  Since the square of a square root should return the same
        number, and the exponential of a log should return the same number,
        we should get back a 2 from this function of functions.  By taking
        the integer result and subtracting from 2, we are exposing the
        roundoff errors.  These simple functions are at the heart of
        statistical calculations.
    """

    def test_rounding0(self):
        """ W.II.A.0. Print ROUND with only one digit.

            You should get the numbers 1 to 9.  Many language compilers,
            such as Turbo Pascal and Lattice C, fail this test (they round
            numbers inconsistently). Needless to say, statical packages
            written in these languages may fail the test as well.
        """
        for i in range(0,9):
            y = round(ROUND[i])
            assert_equal(y,i+1)

    def test_rounding1(self):
        """ W.II.A.1. Y = INT(2.6*7 -0.2) (Y should be 18)"""
        y = int(2.6*7 -0.2)
        assert_equal(y, 18)

    def test_rounding2(self):
        """ W.II.A.2. Y = 2-INT(EXP(LOG(SQR(2)*SQR(2))))   (Y should be 0)"""
        y=2-int(np.exp(np.log(np.sqrt(2.)*np.sqrt(2.))))
        assert_equal(y,0)

    def test_rounding3(self):
        """ W.II.A.3. Y = INT(3-EXP(LOG(SQR(2)*SQR(2))))    (Y should be 1)"""
        y=(int(round((3-np.exp(np.log(np.sqrt(2.0)*np.sqrt(2.0)))))))
        assert_equal(y,1)

class TestBasicStats(TestCase):
    """ W.II.C. Compute basic statistic on all the variables.

        The means should be the fifth value of all the variables (case FIVE).
        The standard deviations should be "undefined" or missing for MISS,
        0 for ZERO, and 2.738612788 (times 10 to a power) for all the other variables.
        II. C. Basic Statistics
    """

    def test_meanX(self):
        y = np.mean(X)
        assert_almost_equal(y, 5.0)

    def test_stdX(self):
        y = np.std(X, ddof=1)
        assert_almost_equal(y, 2.738612788)

    def test_tmeanX(self):
        y = stats.tmean(X, (2, 8), (True, True))
        assert_almost_equal(y, 5.0)

    def test_tvarX(self):
        y = stats.tvar(X, (2, 8), (True, True))
        assert_almost_equal(y, 4.6666666666666661)

    def test_tstdX(self):
        y = stats.tstd(X, (2, 8), (True, True))
        assert_almost_equal(y, 2.1602468994692865)

    def test_meanZERO(self):
        y = np.mean(ZERO)
        assert_almost_equal(y, 0.0)

    def test_stdZERO(self):
        y = np.std(ZERO, ddof=1)
        assert_almost_equal(y, 0.0)

##    Really need to write these tests to handle missing values properly
##    def test_meanMISS(self):
##        y = np.mean(MISS)
##        assert_almost_equal(y, 0.0)
##
##    def test_stdMISS(self):
##        y = stats.stdev(MISS)
##        assert_almost_equal(y, 0.0)

    def test_meanBIG(self):
        y = np.mean(BIG)

        assert_almost_equal(y, 99999995.00)

    def test_stdBIG(self):
        y = np.std(BIG, ddof=1)
        assert_almost_equal(y, 2.738612788)

    def test_meanLITTLE(self):
        y = np.mean(LITTLE)
        assert_approx_equal(y, 0.999999950)

    def test_stdLITTLE(self):
        y = np.std(LITTLE, ddof=1)
        assert_approx_equal(y, 2.738612788e-8)

    def test_meanHUGE(self):
        y = np.mean(HUGE)
        assert_approx_equal(y, 5.00000e+12)

    def test_stdHUGE(self):
        y = np.std(HUGE, ddof=1)
        assert_approx_equal(y, 2.738612788e12)

    def test_meanTINY(self):
        y = np.mean(TINY)
        assert_almost_equal(y, 0.0)

    def test_stdTINY(self):
        y = np.std(TINY, ddof=1)
        assert_almost_equal(y, 0.0)

    def test_meanROUND(self):
        y = np.mean(ROUND)
        assert_approx_equal(y, 4.500000000)

    def test_stdROUND(self):
        y = np.std(ROUND, ddof=1)
        assert_approx_equal(y, 2.738612788)

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
        """Check nanmean when no values are nan."""
        m = stats.nanmean(X)
        assert_approx_equal(m, X[4])

    def test_nanmean_some(self):
        """Check nanmean when some values only are nan."""
        m = stats.nanmean(self.Xsome)
        assert_approx_equal(m, 5.5)

    def test_nanmean_all(self):
        """Check nanmean when all values are nan."""
        m = stats.nanmean(self.Xall)
        assert np.isnan(m)

    def test_nanstd_none(self):
        """Check nanstd when no values are nan."""
        s = stats.nanstd(self.X)
        assert_approx_equal(s, np.std(self.X, ddof=1))

    def test_nanstd_some(self):
        """Check nanstd when some values only are nan."""
        s = stats.nanstd(self.Xsome)
        assert_approx_equal(s, np.std(self.Xsomet, ddof=1))

    def test_nanstd_all(self):
        """Check nanstd when all values are nan."""
        s = stats.nanstd(self.Xall)
        assert np.isnan(s)

    def test_nanstd_negative_axis(self):
        x = np.array([1, 2, 3])
        assert_equal(stats.nanstd(x, -1), 1)

    def test_nanmedian_none(self):
        """Check nanmedian when no values are nan."""
        m = stats.nanmedian(self.X)
        assert_approx_equal(m, np.median(self.X))

    def test_nanmedian_some(self):
        """Check nanmedian when some values only are nan."""
        m = stats.nanmedian(self.Xsome)
        assert_approx_equal(m, np.median(self.Xsomet))

    def test_nanmedian_all(self):
        """Check nanmedian when all values are nan."""
        m = stats.nanmedian(self.Xall)
        assert np.isnan(m)

class TestCorr(TestCase):
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

##    W.II.E.  Tabulate X against X, using BIG as a case weight.  The values
##    should appear on the diagonal and the total should be 899999955.
##    If the table cannot hold these values, forget about working with
##    census data.  You can also tabulate HUGE against TINY.  There is no
##    reason a tabulation program should not be able to digtinguish
##    different values regardless of their magnitude.

### I need to figure out how to do this one.


class TestRegression(TestCase):
    def test_linregressBIGX(self):
        """ W.II.F.  Regress BIG on X.

            The constant should be 99999990 and the regression coefficient should be 1.
        """
        y = stats.linregress(X,BIG)
        intercept = y[1]
        r=y[2]
        assert_almost_equal(intercept,99999990)
        assert_almost_equal(r,1.0)

##     W.IV.A. Take the NASTY dataset above.  Use the variable X as a
##     basis for computing polynomials.  Namely, compute X1=X, X2=X*X,
##     X3=X*X*X, and so on up to 9 products.  Use the algebraic
##     transformation language within the statistical package itself.  You
##     will end up with 9 variables.  Now regress X1 on X2-X9 (a perfect
##     fit).  If the package balks (singular or roundoff error messages),
##     try X1 on X2-X8, and so on.  Most packages cannot handle more than
##     a few polynomials.
##     Scipy's stats.py does not seem to handle multiple linear regression
##     The datasets X1 . . X9 are at the top of the file.


    def test_regressXX(self):
        """ W.IV.B.  Regress X on X.

            The constant should be exactly 0 and the regression coefficient should be 1.
            This is a perfectly valid regression.  The program should not complain.
        """
        y = stats.linregress(X,X)
        intercept = y[1]
        r=y[2]
        assert_almost_equal(intercept,0.0)
        assert_almost_equal(r,1.0)
##     W.IV.C. Regress X on BIG and LITTLE (two predictors).  The program
##     should tell you that this model is "singular" because BIG and
##     LITTLE are linear combinations of each other.  Cryptic error
##     messages are unacceptable here.  Singularity is the most
##     fundamental regression error.
### Need to figure out how to handle multiple linear regression.  Not obvious

    def test_regressZEROX(self):
        """ W.IV.D. Regress ZERO on X.

            The program should inform you that ZERO has no variance or it should
            go ahead and compute the regression and report a correlation and
            total sum of squares of exactly 0.
        """
        y = stats.linregress(X,ZERO)
        intercept = y[1]
        r=y[2]
        assert_almost_equal(intercept,0.0)
        assert_almost_equal(r,0.0)

    def test_regress_simple(self):
        """Regress a line with sinusoidal noise."""
        x = np.linspace(0, 100, 100)
        y = 0.2 * np.linspace(0, 100, 100) + 10
        y += np.sin(np.linspace(0, 20, 100))

        res = stats.linregress(x, y)
        assert_almost_equal(res[4], 2.3957814497838803e-3) #4.3609875083149268e-3)

    def test_linregress(self):
        '''compared with multivariate ols with pinv'''
        x = np.arange(11)
        y = np.arange(5,16)
        y[[(1),(-2)]] -= 1
        y[[(0),(-1)]] += 1

        res = (1.0, 5.0, 0.98229948625750, 7.45259691e-008, 0.063564172616372733)
        assert_array_almost_equal(stats.linregress(x,y),res,decimal=14)

class TestHistogram(TestCase):
    """ Tests that histogram works as it should, and keeps old behaviour
    """
    # what is untested:
    # - multidimensional arrays (since 'a' is ravel'd as the first line in the method)
    # - very large arrays
    # - Nans, Infs, empty and otherwise bad inputs

    # sample arrays to test the histogram with
    low_values = np.array([0.2, 0.3, 0.4, 0.5, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2],
                          dtype=float) # 11 values
    high_range = np.array([2, 3, 4, 2, 21, 32, 78, 95, 65, 66, 66, 66, 66, 4],
                          dtype=float) # 14 values
    low_range = np.array([2, 3, 3, 2, 3, 2.4, 2.1, 3.1, 2.9, 2.6, 2.7, 2.8, 2.2, 2.001],
                         dtype=float) # 14 values
    few_values = np.array([2.0, 3.0, -1.0, 0.0], dtype=float) # 4 values

    def test_simple(self):
        """ Tests that each of the tests works as expected with default params
        """
        # basic tests, with expected results (no weighting)
        # results taken from the previous (slower) version of histogram
        basic_tests = ((self.low_values, (np.array([ 1.,  1.,  1.,  2.,  2.,
                                                     1.,  1.,  0.,  1.,  1.]),
                                          0.14444444444444446, 0.11111111111111112, 0)),
                       (self.high_range, (np.array([ 5.,  0.,  1.,  1.,  0.,
                                                     0.,  5.,  1.,  0.,  1.]),
                                          -3.1666666666666661, 10.333333333333332, 0)),
                       (self.low_range, (np.array([ 3.,  1.,  1.,  1.,  0.,  1.,
                                                    1.,  2.,  3.,  1.]),
                                         1.9388888888888889, 0.12222222222222223, 0)),
                       (self.few_values, (np.array([ 1.,  0.,  1.,  0.,  0.,  0.,
                                                     0.,  1.,  0.,  1.]),
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
        """ Tests that weights give expected histograms
        """
        # basic tests, with expected results, given a set of weights
        # weights used (first n are used for each test, where n is len of array) (14 values)
        weights = np.array([1., 3., 4.5, 0.1, -1.0, 0.0, 0.3, 7.0, 103.2, 2, 40, 0, 0, 1])
        # results taken from the numpy version of histogram
        basic_tests = ((self.low_values, (np.array([  4.0,  0.0,  4.5,  -0.9,  0.0,
                                                      0.3,110.2,  0.0,  0.0,  42.0]),
                                          0.2, 0.1, 0)),
                       (self.high_range, (np.array([  9.6,  0. ,  -1. ,  0. ,  0. ,
                                                      0. ,145.2,   0. ,  0.3,  7. ]),
                                          2.0, 9.3, 0)),
                       (self.low_range, (np.array([ 2.4,  0. ,  0. ,   0. ,  0. ,
                                                    2. , 40. ,  0. , 103.2, 13.5]),
                                         2.0, 0.11, 0)),
                       (self.few_values, (np.array([ 4.5,  0. ,  0.1,  0. ,  0. ,  0. ,
                                                     0. ,  1. ,  0. ,  3. ]),
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
        """ Tests that reducing the number of bins produces expected results
        """
        # basic tests, with expected results (no weighting),
        # except number of bins is halved to 5
        # results taken from the previous (slower) version of histogram
        basic_tests = ((self.low_values, (np.array([ 2.,  3.,  3.,  1.,  2.]),
                                          0.075000000000000011, 0.25, 0)),
                       (self.high_range, (np.array([ 5.,  2.,  0.,  6.,  1.]),
                                          -9.625, 23.25, 0)),
                       (self.low_range, (np.array([ 4.,  2.,  1.,  3.,  4.]),
                                         1.8625, 0.27500000000000002, 0)),
                       (self.few_values, (np.array([ 1.,  1.,  0.,  1.,  1.]),
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
        """ Tests that increasing the number of bins produces expected results
        """
        # basic tests, with expected results (no weighting),
        # except number of bins is double to 20
        # results taken from the previous (slower) version of histogram
        basic_tests = ((self.low_values, (np.array([ 1.,  0.,  1.,  0.,  1.,
                                                     0.,  2.,  0.,  1.,  0.,
                                                     1.,  1.,  0.,  1.,  0.,
                                                     0.,  0.,  1.,  0.,  1.]),
                                          0.1736842105263158, 0.052631578947368418, 0)),
                       (self.high_range, (np.array([ 5.,  0.,  0.,  0.,  1.,
                                                     0.,  1.,  0.,  0.,  0.,
                                                     0.,  0.,  0.,  5.,  0.,
                                                     0.,  1.,  0.,  0.,  1.]),
                                          -0.44736842105263142, 4.8947368421052628, 0)),
                       (self.low_range, (np.array([ 3.,  0.,  1.,  1.,  0.,  0.,
                                                    0.,  1.,  0.,  0.,  1.,  0.,
                                                    1., 0.,  1.,  0.,  1.,  3.,
                                                    0.,  1.]),
                                         1.9710526315789474, 0.057894736842105263, 0)),
                       (self.few_values, (np.array([ 1.,  0.,  0.,  0.,  0.,  1.,
                                                     0.,  0.,  0.,  0.,  0.,  0.,
                                                     0.,  0.,  1.,  0.,  0.,  0.,
                                                     0.,  1.]),
                                          -1.1052631578947367, 0.21052631578947367, 0)),
                       )
        for inputs, expected_results in basic_tests:
            given_results = stats.histogram(inputs, numbins=20)
            assert_array_almost_equal(expected_results[0], given_results[0],
                                      decimal=2)
            for i in range(1, 4):
                assert_almost_equal(expected_results[i], given_results[i],
                                    decimal=2)


# Utility

def compare_results(res,desired):
    for i in range(len(desired)):
        assert_array_equal(res[i],desired[i])


##################################################
### Test for sum

class TestGMean(TestCase):

    def test_1D_list(self):
        a = (1,2,3,4)
        actual= stats.gmean(a)
        desired = power(1*2*3*4,1./4.)
        assert_almost_equal(actual, desired,decimal=14)

        desired1 = stats.gmean(a,axis=-1)
        assert_almost_equal(actual, desired1, decimal=14)

    def test_1D_array(self):
        a = array((1,2,3,4), float32)
        actual= stats.gmean(a)
        desired = power(1*2*3*4,1./4.)
        assert_almost_equal(actual, desired, decimal=7)

        desired1 = stats.gmean(a,axis=-1)
        assert_almost_equal(actual, desired1, decimal=7)

    def test_2D_array_default(self):
        a = array(((1,2,3,4),
                   (1,2,3,4),
                   (1,2,3,4)))
        actual= stats.gmean(a)
        desired = array((1,2,3,4))
        assert_array_almost_equal(actual, desired, decimal=14)

        desired1 = stats.gmean(a,axis=0)
        assert_array_almost_equal(actual, desired1, decimal=14)

    def test_2D_array_dim1(self):
        a = array(((1,2,3,4),
                   (1,2,3,4),
                   (1,2,3,4)))
        actual= stats.gmean(a, axis=1)
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
        actual= stats.hmean(a)
        desired =  4. / (1./1 + 1./2 + 1./3 + 1./4)
        assert_almost_equal(actual, desired, decimal=14)

        desired1 = stats.hmean(array(a),axis=-1)
        assert_almost_equal(actual, desired1, decimal=14)
    def test_1D_array(self):
        a = array((1,2,3,4), float64)
        actual= stats.hmean(a)
        desired =  4. / (1./1 + 1./2 + 1./3 + 1./4)
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


class TestMean(TestCase):
    def test_basic(self):
        a = [3,4,5,10,-3,-5,6]
        af = [3.,4,5,10,-3,-5,-6]
        Na = len(a)
        Naf = len(af)
        mn1 = 0.0
        for el in a:
            mn1 += el / float(Na)
        assert_almost_equal(np.mean(a),mn1,11)
        mn2 = 0.0
        for el in af:
            mn2 += el / float(Naf)
        assert_almost_equal(np.mean(af),mn2,11)

    def test_2d(self):
        a = [[1.0, 2.0, 3.0],
             [2.0, 4.0, 6.0],
             [8.0, 12.0, 7.0]]
        A = array(a)
        N1, N2 = (3, 3)
        mn1 = zeros(N2, dtype=float)
        for k in range(N1):
            mn1 += A[k,:] / N1
        assert_almost_equal(np.mean(a, axis=0), mn1, decimal=13)
        mn2 = zeros(N1, dtype=float)
        for k in range(N2):
            mn2 += A[:,k]
        mn2 /= N2
        assert_almost_equal(np.mean(a, axis=1), mn2, decimal=13)

    def test_ravel(self):
        a = rand(5,3,5)
        A = 0
        for val in ravel(a):
            A += val
        assert_almost_equal(np.mean(a,axis=None),A/(5*3.0*5))

class TestPercentile(TestCase):
    def setUp(self):
        self.a1 = [3,4,5,10,-3,-5,6]
        self.a2 = [3,-6,-2,8,7,4,2,1]
        self.a3 = [3.,4,5,10,-3,-5,-6,7.0]

    def test_median(self):
        assert_equal(np.median(self.a1), 4)
        assert_equal(np.median(self.a2), 2.5)
        assert_equal(np.median(self.a3), 3.5)

    def test_percentile(self):
        x = arange(8) * 0.5
        assert_equal(stats.scoreatpercentile(x, 0), 0.)
        assert_equal(stats.scoreatpercentile(x, 100), 3.5)
        assert_equal(stats.scoreatpercentile(x, 50), 1.75)

    def test_2D(self):
        x = array([[1, 1, 1],
                   [1, 1, 1],
                   [4, 4, 3],
                   [1, 1, 1],
                   [1, 1, 1]])
        assert_array_equal(stats.scoreatpercentile(x,50),
                           [1,1,1])


class TestStd(TestCase):
    def test_basic(self):
        a = [3,4,5,10,-3,-5,6]
        b = [3,4,5,10,-3,-5,-6]
        assert_almost_equal(np.std(a, ddof=1),5.2098807225172772,11)
        assert_almost_equal(np.std(b, ddof=1),5.9281411203561225,11)

    def test_2d(self):
        a = [[1.0, 2.0, 3.0],
             [2.0, 4.0, 6.0],
             [8.0, 12.0, 7.0]]
        b1 = array((3.7859388972001824, 5.2915026221291814,
                    2.0816659994661335))
        b2 = array((1.0,2.0,2.64575131106))
        assert_array_almost_equal(np.std(a,ddof=1,axis=0),b1,11)
        assert_array_almost_equal(np.std(a,ddof=1,axis=1),b2,11)


class TestCMedian(TestCase):
    def test_basic(self):
        data = [1,2,3,1,5,3,6,4,3,2,4,3,5,2.0]
        assert_almost_equal(stats.cmedian(data,5),3.2916666666666665)
        assert_almost_equal(stats.cmedian(data,3),3.083333333333333)
        assert_almost_equal(stats.cmedian(data),3.0020020020020022)

class TestMedian(TestCase):
    def test_basic(self):
        data1 = [1,3,5,2,3,1,19,-10,2,4.0]
        data2 = [3,5,1,10,23,-10,3,-2,6,8,15]
        assert_almost_equal(np.median(data1),2.5)
        assert_almost_equal(np.median(data2),5)

    def test_basic2(self):
        a1 = [3,4,5,10,-3,-5,6]
        a2 = [3,-6,-2,8,7,4,2,1]
        a3 = [3.,4,5,10,-3,-5,-6,7.0]
        assert_equal(np.median(a1),4)
        assert_equal(np.median(a2),2.5)
        assert_equal(np.median(a3),3.5)

    def test_axis(self):
        """Regression test for #760."""
        a1 = np.array([[3,4,5], [10,-3,-5]])
        assert_equal(np.median(a1), 3.5)
        assert_equal(np.median(a1, axis=0), np.array([6.5, 0.5, 0.]))
        assert_equal(np.median(a1, axis=-1), np.array([4., -3]))

class TestMode(TestCase):
    def test_basic(self):
        data1 = [3,5,1,10,23,3,2,6,8,6,10,6]
        vals = stats.mode(data1)
        assert_almost_equal(vals[0][0],6)
        assert_almost_equal(vals[1][0],3)


class TestVariability(TestCase):
    """  Comparison numbers are found using R v.1.5.1
         note that length(testcase) = 4
    """
    testcase = [1,2,3,4]
    def test_std(self):
        y = np.std(self.testcase, ddof=1)
        assert_approx_equal(y,1.290994449)

    def test_var(self):
        """
        var(testcase) = 1.666666667 """
        #y = stats.var(self.shoes[0])
        #assert_approx_equal(y,6.009)
        y = np.var(self.testcase)
        assert_approx_equal(y,1.25)
        y = np.var(self.testcase, ddof=1)
        assert_approx_equal(y,1.666666667)

    def test_samplevar(self):
        """
        R does not have 'samplevar' so the following was used
        var(testcase)*(4-1)/4  where 4 = length(testcase)
        """
        #y = stats.samplevar(self.shoes[0])
        #assert_approx_equal(y,5.4081)
        y = stats.samplevar(self.testcase)
        assert_approx_equal(y,1.25)

    def test_samplestd(self):
        #y = stats.samplestd(self.shoes[0])
        #assert_approx_equal(y,2.325532197)
        y = stats.samplestd(self.testcase)
        assert_approx_equal(y,1.118033989)

    def test_signaltonoise(self):
        """
        this is not in R, so used
        mean(testcase,axis=0)/(sqrt(var(testcase)*3/4)) """
        #y = stats.signaltonoise(self.shoes[0])
        #assert_approx_equal(y,4.5709967)
        y = stats.signaltonoise(self.testcase)
        assert_approx_equal(y,2.236067977)

    def test_stderr(self):
        """
        this is not in R, so used
        sqrt(var(testcase))/sqrt(4)
        """
##        y = stats.stderr(self.shoes[0])
##        assert_approx_equal(y,0.775177399)
        y = stats.stderr(self.testcase)
        assert_approx_equal(y,0.6454972244)
    def test_sem(self):
        """
        this is not in R, so used
        sqrt(var(testcase)*3/4)/sqrt(3)
        """
        #y = stats.sem(self.shoes[0])
        #assert_approx_equal(y,0.775177399)
        y = stats.sem(self.testcase)
        assert_approx_equal(y,0.6454972244)

    def test_z(self):
        """
        not in R, so used
        (10-mean(testcase,axis=0))/sqrt(var(testcase)*3/4)
        """
        y = stats.z(self.testcase,np.mean(self.testcase, axis=0))
        assert_almost_equal(y,0.0)

    def test_zs(self):
        """
        not in R, so tested by using
        (testcase[i]-mean(testcase,axis=0))/sqrt(var(testcase)*3/4)
        """
        y = stats.zs(self.testcase)
        desired = ([-1.3416407864999, -0.44721359549996 , 0.44721359549996 , 1.3416407864999])
        assert_array_almost_equal(desired,y,decimal=12)

    def test_zmap(self):
        """
        not in R, so tested by using
        (testcase[i]-mean(testcase,axis=0))/sqrt(var(testcase)*3/4)
        copied from test_zs
        """
        y = stats.zmap(self.testcase,self.testcase)
        desired = ([-1.3416407864999, -0.44721359549996 , 0.44721359549996 , 1.3416407864999])
        assert_array_almost_equal(desired,y,decimal=12)

    def test_zscore(self):
        """
        not in R, so tested by using
        (testcase[i]-mean(testcase,axis=0))/sqrt(var(testcase)*3/4)
        copied from test_zs as regression test for new function
        """
        y = stats.zscore(self.testcase)
        desired = ([-1.3416407864999, -0.44721359549996 , 0.44721359549996 , 1.3416407864999])
        assert_array_almost_equal(desired,y,decimal=12)

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
    testmathworks = [1.165 , 0.6268, 0.0751, 0.3516, -0.6965]
    def test_moment(self):
        """
        mean((testcase-mean(testcase))**power,axis=0),axis=0))**power))"""
        y = stats.moment(self.testcase,1)
        assert_approx_equal(y,0.0,10)
        y = stats.moment(self.testcase,2)
        assert_approx_equal(y,1.25)
        y = stats.moment(self.testcase,3)
        assert_approx_equal(y,0.0)
        y = stats.moment(self.testcase,4)
        assert_approx_equal(y,2.5625)
    def test_variation(self):
        """
        variation = samplestd/mean """
##        y = stats.variation(self.shoes[0])
##        assert_approx_equal(y,21.8770668)
        y = stats.variation(self.testcase)
        assert_approx_equal(y,0.44721359549996, 10)

    def test_skewness(self):
        """
        sum((testmathworks-mean(testmathworks,axis=0))**3,axis=0)/
            ((sqrt(var(testmathworks)*4/5))**3)/5
        """
        y = stats.skew(self.testmathworks)
        assert_approx_equal(y,-0.29322304336607,10)
        y = stats.skew(self.testmathworks,bias=0)
        assert_approx_equal(y,-0.437111105023940,10)
        y = stats.skew(self.testcase)
        assert_approx_equal(y,0.0,10)

    def test_skewness_scalar(self):
        """
        `skew` must return a scalar for 1-dim input
        """
        assert_equal(stats.skew(arange(10)), 0.0)

    def test_kurtosis(self):
        """
            sum((testcase-mean(testcase,axis=0))**4,axis=0)/((sqrt(var(testcase)*3/4))**4)/4
            sum((test2-mean(testmathworks,axis=0))**4,axis=0)/((sqrt(var(testmathworks)*4/5))**4)/5
            Set flags for axis = 0 and
            fisher=0 (Pearson's defn of kurtosis for compatiability with Matlab)
        """
        y = stats.kurtosis(self.testmathworks,0,fisher=0,bias=1)
        assert_approx_equal(y, 2.1658856802973,10)

        # Note that MATLAB has confusing docs for the following case
        #  kurtosis(x,0) gives an unbiased estimate of Pearson's skewness
        #  kurtosis(x)  gives a biased estimate of Fisher's skewness (Pearson-3)
        #  The MATLAB docs imply that both should give Fisher's
        y = stats.kurtosis(self.testmathworks,fisher=0,bias=0)
        assert_approx_equal(y, 3.663542721189047,10)
        y = stats.kurtosis(self.testcase,0,0)
        assert_approx_equal(y,1.64)

    def test_kurtosis_array_scalar(self):
        assert_equal(type(stats.kurtosis([1,2,3])), float)

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

# Hypothesis test tests
class TestStudentTest(TestCase):
    X1 = np.array([-1, 0, 1])
    X2 = np.array([0, 1, 2])
    T1_0 = 0
    P1_0 = 1
    T1_1 = -1.732051
    P1_1 = 0.2254033
    T1_2 = -3.464102
    P1_2 =  0.0741799
    T2_0 = 1.732051
    P2_0 = 0.2254033
    def test_onesample(self):
        t, p = stats.ttest_1samp(self.X1, 0)

        assert_array_almost_equal(t, self.T1_0)
        assert_array_almost_equal(p, self.P1_0)

        t, p = stats.ttest_1samp(self.X2, 0)

        assert_array_almost_equal(t, self.T2_0)
        assert_array_almost_equal(p, self.P2_0)

        t, p = stats.ttest_1samp(self.X1, 1)

        assert_array_almost_equal(t, self.T1_1)
        assert_array_almost_equal(p, self.P1_1)

        t, p = stats.ttest_1samp(self.X1, 2)

        assert_array_almost_equal(t, self.T1_2)
        assert_array_almost_equal(p, self.P1_2)

def test_scoreatpercentile():
    assert_equal(stats.scoreatpercentile(range(10), 50), 4.5)
    assert_equal(stats.scoreatpercentile(range(10), 50, (2,7)), 4.5)
    assert_equal(stats.scoreatpercentile(range(100), 50, (1,8)), 4.5)

    assert_equal(stats.scoreatpercentile(np.array([1, 10 ,100]),
                                         50, (10,100)),
                 55)
    assert_equal(stats.scoreatpercentile(np.array([1, 10 ,100]),
                                         50, (1,10)),
                 5.5)

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

    #larger numbers
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

    #boundaries
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

    #out of bounds
    for (kind, score, result) in [('rank', 200, 100.0),
                                  ('mean', 200, 100.0),
                                  ('mean', 0, 0.0)]:
        yield assert_equal, \
              pcos([10, 20, 30, 50, 60, 70, 80, 90, 100, 110],
                   score, kind=kind), result


def test_friedmanchisquare():
    # see ticket:113
    # verified with matlab and R
    #From Demsar "Statistical Comparisons of Classifiers over Multiple Data Sets"
    #2006, Xf=9.28 (no tie handling, tie corrected Xf >=9.28)
    x1 = [array([0.763, 0.599, 0.954, 0.628, 0.882, 0.936, 0.661, 0.583,
                 0.775, 1.0, 0.94, 0.619, 0.972, 0.957]),
          array([0.768, 0.591, 0.971, 0.661, 0.888, 0.931, 0.668, 0.583,
                 0.838, 1.0, 0.962, 0.666, 0.981, 0.978]),
          array([0.771, 0.590, 0.968, 0.654, 0.886, 0.916, 0.609, 0.563,
                 0.866, 1.0, 0.965, 0.614, 0.9751, 0.946]),
          array([0.798, 0.569, 0.967, 0.657, 0.898, 0.931, 0.685, 0.625,
                 0.875, 1.0, 0.962, 0.669, 0.975, 0.970])]

    #From "Bioestadistica para las ciencias de la salud" Xf=18.95 p<0.001:
    x2 = [array([4,3,5,3,5,3,2,5,4,4,4,3]),
          array([2,2,1,2,3,1,2,3,2,1,1,3]),
          array([2,4,3,3,4,3,3,4,4,1,2,1]),
          array([3,5,4,3,4,4,3,3,3,4,4,4])]

    #From Jerrorl H. Zar, "Biostatistical Analysis"(example 12.6), Xf=10.68, 0.005 < p < 0.01:
    #Probability from this example is inexact using Chisquare aproximation of Friedman Chisquare.
    x3 = [array([7.0,9.9,8.5,5.1,10.3]),
          array([5.3,5.7,4.7,3.5,7.7]),
          array([4.9,7.6,5.5,2.8,8.4]),
          array([8.8,8.9,8.1,3.3,9.1])]


    assert_array_almost_equal(stats.friedmanchisquare(x1[0],x1[1],x1[2],x1[3]),(10.2283464566929, 0.0167215803284414))
    assert_array_almost_equal(stats.friedmanchisquare(x2[0],x2[1],x2[2],x2[3]),(18.9428571428571, 0.000280938375189499))
    assert_array_almost_equal(stats.friedmanchisquare(x3[0],x3[1],x3[2],x3[3]),(10.68, 0.0135882729582176))
    np.testing.assert_raises(ValueError, stats.friedmanchisquare,x3[0],x3[1])

    # test using mstats
    assert_array_almost_equal(stats.mstats.friedmanchisquare(x1[0],x1[1],x1[2],x1[3]),(10.2283464566929, 0.0167215803284414))
    # the following fails
    #assert_array_almost_equal(stats.mstats.friedmanchisquare(x2[0],x2[1],x2[2],x2[3]),(18.9428571428571, 0.000280938375189499))
    assert_array_almost_equal(stats.mstats.friedmanchisquare(x3[0],x3[1],x3[2],x3[3]),(10.68, 0.0135882729582176))
    np.testing.assert_raises(ValueError,stats.mstats.friedmanchisquare,x3[0],x3[1])

def test_kstest():
    #from numpy.testing import assert_almost_equal

    # comparing with values from R
    x = np.linspace(-1,1,9)
    D,p = stats.kstest(x,'norm')
    assert_almost_equal( D, 0.15865525393145705, 12)
    assert_almost_equal( p, 0.95164069201518386, 1)

    x = np.linspace(-15,15,9)
    D,p = stats.kstest(x,'norm')
    assert_almost_equal( D, 0.44435602715924361, 15)
    assert_almost_equal( p, 0.038850140086788665, 8)

    # the following tests rely on deterministicaly replicated rvs
    np.random.seed(987654321)
    x = stats.norm.rvs(loc=0.2, size=100)
    D,p = stats.kstest(x, 'norm', mode='asymp')
    assert_almost_equal( D, 0.12464329735846891, 15)
    assert_almost_equal( p, 0.089444888711820769, 15)
    assert_almost_equal( np.array(stats.kstest(x, 'norm', mode='asymp')),
                np.array((0.12464329735846891, 0.089444888711820769)), 15)
    assert_almost_equal( np.array(stats.kstest(x,'norm', alternative = 'less')),
                np.array((0.12464329735846891, 0.040989164077641749)), 15)
    # this 'greater' test fails with precision of decimal=14
    assert_almost_equal( np.array(stats.kstest(x,'norm', alternative = 'greater')),
                np.array((0.0072115233216310994, 0.98531158590396228)), 12)

    #missing: no test that uses *args


def test_ks_2samp():
    #exact small sample solution
    data1 = np.array([1.0,2.0])
    data2 = np.array([1.0,2.0,3.0])
    assert_almost_equal(np.array(stats.ks_2samp(data1+0.01,data2)),
                np.array((0.33333333333333337, 0.99062316386915694)))
    assert_almost_equal(np.array(stats.ks_2samp(data1-0.01,data2)),
                np.array((0.66666666666666674, 0.42490954988801982)))
    #these can also be verified graphically
    assert_almost_equal(
        np.array(stats.ks_2samp(np.linspace(1,100,100),
                              np.linspace(1,100,100)+2+0.1)),
        np.array((0.030000000000000027, 0.99999999996005062)))
    assert_almost_equal(
        np.array(stats.ks_2samp(np.linspace(1,100,100),
                              np.linspace(1,100,100)+2-0.1)),
        np.array((0.020000000000000018, 0.99999999999999933)))
    #these are just regression tests
    assert_almost_equal(
        np.array(stats.ks_2samp(np.linspace(1,100,100),
                              np.linspace(1,100,110)+20.1)),
        np.array((0.21090909090909091, 0.015880386730710221)))
    assert_almost_equal(
        np.array(stats.ks_2samp(np.linspace(1,100,100),
                              np.linspace(1,100,110)+20-0.1)),
        np.array((0.20818181818181825, 0.017981441789762638)))

def test_ttest_rel():
    #regression test
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

    #test on 3 dimensions
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

    #test zero division problem
    t,p = stats.ttest_rel([0,0,0],[1,1,1])
    assert_equal((np.abs(t),p), (np.inf, 0))
    assert_almost_equal(stats.ttest_rel([0,0,0], [0,0,0]), (1.0, 0.42264973081037421))

    #check that nan in input array result in nan output
    anan = np.array([[1,np.nan],[-1,1]])
    assert_equal(stats.ttest_ind(anan, np.zeros((2,2))),([0, np.nan], [1,np.nan]))


def test_ttest_ind():
    #regression test
    tr = 1.0912746897927283
    pr = 0.27647818616351882
    tpr = ([tr,-tr],[pr,pr])

    rvs2 = np.linspace(1,100,100)
    rvs1 = np.linspace(5,105,100)
    rvs1_2D = np.array([rvs1, rvs2])
    rvs2_2D = np.array([rvs2, rvs1])

    t,p = stats.ttest_ind(rvs1, rvs2, axis=0)
    assert_array_almost_equal([t,p],(tr,pr))
    t,p = stats.ttest_ind(rvs1_2D.T, rvs2_2D.T, axis=0)
    assert_array_almost_equal([t,p],tpr)
    t,p = stats.ttest_ind(rvs1_2D, rvs2_2D, axis=1)
    assert_array_almost_equal([t,p],tpr)

    #test on 3 dimensions
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

    #test zero division problem
    t,p = stats.ttest_ind([0,0,0],[1,1,1])
    assert_equal((np.abs(t),p), (np.inf, 0))
    assert_almost_equal(stats.ttest_ind([0,0,0], [0,0,0]), (1.0, 0.37390096630005898))

    #check that nan in input array result in nan output
    anan = np.array([[1,np.nan],[-1,1]])
    assert_equal(stats.ttest_ind(anan, np.zeros((2,2))),([0, np.nan], [1,np.nan]))




def test_ttest_1samp_new():
    n1, n2, n3 = (10,15,20)
    rvn1 = stats.norm.rvs(loc=5,scale=10,size=(n1,n2,n3))
    rvn2 = stats.norm.rvs(loc=5,scale=10,size=(n1,n2,n3))

    #check multidimensional array and correct axis handling
    #deterministic rvn1 and rvn2 would be better as in test_ttest_rel
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

    #test zero division problem
    t,p = stats.ttest_1samp([0,0,0], 1)
    assert_equal((np.abs(t),p), (np.inf, 0))
    assert_almost_equal(stats.ttest_1samp([0,0,0], 0), (1.0, 0.42264973081037421))

    #check that nan in input array result in nan output
    anan = np.array([[1,np.nan],[-1,1]])
    assert_equal(stats.ttest_1samp(anan, 0),([0, np.nan], [1,np.nan]))

def test_describe():
    x = np.vstack((np.ones((3,4)),2*np.ones((2,4))))
    nc, mmc = (5, ([ 1.,  1.,  1.,  1.], [ 2.,  2.,  2.,  2.]))
    mc = np.array([ 1.4,  1.4,  1.4,  1.4])
    vc = np.array([ 0.3,  0.3,  0.3,  0.3])
    skc = [0.40824829046386357]*4
    kurtc = [-1.833333333333333]*4
    n, mm, m, v, sk, kurt = stats.describe(x)
    assert_equal(n, nc)
    assert_equal(mm, mmc)
    assert_equal(m, mc)
    assert_equal(v, vc)
    assert_array_almost_equal(sk, skc, decimal=13) #not sure about precision
    assert_array_almost_equal(kurt, kurtc, decimal=13)
    n, mm, m, v, sk, kurt = stats.describe(x.T, axis=1)
    assert_equal(n, nc)
    assert_equal(mm, mmc)
    assert_equal(m, mc)
    assert_equal(v, vc)
    assert_array_almost_equal(sk, skc, decimal=13) #not sure about precision
    assert_array_almost_equal(kurt, kurtc, decimal=13)

def test_normalitytests():
    # numbers verified with R: dagoTest in package fBasics
    st_normal, st_skew, st_kurt = (3.92371918, 1.98078826, -0.01403734)
    pv_normal, pv_skew, pv_kurt = (0.14059673, 0.04761502,  0.98880019)
    x = np.array((-2,-1,0,1,2,3)*4)**2
    yield assert_array_almost_equal, stats.normaltest(x), (st_normal, pv_normal)
    yield assert_array_almost_equal, stats.skewtest(x), (st_skew, pv_skew)
    yield assert_array_almost_equal, stats.kurtosistest(x), (st_kurt, pv_kurt)


def mannwhitneyu():
    x = np.array([ 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
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

    y = np.array([ 1., 1., 1., 1., 1., 1., 1., 2., 1., 2., 1., 1., 1.,
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
    #p-value verified with matlab and R to 5 significant digits
    assert_array_almost_equal(stats.stats.mannwhitneyu(x,y),
                    (16980.5, 2.8214327656317373e-005), decimal=12)



def test_pointbiserial():
    # copied from mstats tests removing nans
    x = [1,0,1,1,1,1,0,1,0,0,0,1,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,
         0,0,0,0,1]
    y = [14.8,13.8,12.4,10.1,7.1,6.1,5.8,4.6,4.3,3.5,3.3,3.2,3.0,
         2.8,2.8,2.5,2.4,2.3,2.1,1.7,1.7,1.5,1.3,1.3,1.2,1.2,1.1,
         0.8,0.7,0.6,0.5,0.2,0.2,0.1]
    assert_almost_equal(stats.pointbiserialr(x, y)[0], 0.36149, 5)


def test_obrientransform():
    #this is a regression test to check np.var replacement
    #I didn't separately verigy the numbers
    x1 = np.arange(5)
    result = np.array(
      [[  5.41666667,   1.04166667,  -0.41666667,   1.04166667,  5.41666667],
       [ 21.66666667,   4.16666667,  -1.66666667,   4.16666667, 21.66666667]])
    assert_array_almost_equal(stats.obrientransform(x1, 2*x1), result, decimal=8)


class HarMeanTestCase:
    def test_1dlist(self):
        ''' Test a 1d list'''
        a=[10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        b = 34.1417152147
        self.do(a, b)
    def test_1darray(self):
        ''' Test a 1d array'''
        a=np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
        b = 34.1417152147
        self.do(a, b)
    def test_1dma(self):
        ''' Test a 1d masked array'''
        a=np.ma.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
        b = 34.1417152147
        self.do(a, b)
    def test_1dmavalue(self):
        ''' Test a 1d masked array with a masked value'''
        a=np.ma.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
                      mask=[0,0,0,0,0,0,0,0,0,1])
        b = 31.8137186141
        self.do(a, b)

    # Note the next tests use axis=None as default, not axis=0
    def test_2dlist(self):
        ''' Test a 2d list'''
        a=[[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = 38.6696271841
        self.do(a, b)
    def test_2darray(self):
        ''' Test a 2d array'''
        a=[[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = 38.6696271841
        self.do(np.array(a), b)
    def test_2dma(self):
        ''' Test a 2d masked array'''
        a=[[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = 38.6696271841
        self.do(np.ma.array(a), b)
    def test_2daxis0(self):
        ''' Test a 2d list with axis=0'''
        a=[[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = np.array([ 22.88135593,  39.13043478,  52.90076336,  65.45454545])
        self.do(a, b, axis=0)
    def test_2daxis1(self):
        ''' Test a 2d list with axis=1'''
        a=[[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = np.array([  19.2       ,   63.03939962,  103.80078637])
        self.do(a, b, axis=1)
    def test_2dmatrixdaxis0(self):
        ''' Test a 2d list with axis=0'''
        a=[[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = np.matrix([[ 22.88135593,  39.13043478,  52.90076336,  65.45454545]])
        self.do(np.matrix(a), b, axis=0)
    def test_2dmatrixaxis1(self):
        ''' Test a 2d list with axis=1'''
        a=[[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = np.matrix([[  19.2       ,   63.03939962,  103.80078637]]).T
        self.do(np.matrix(a), b, axis=1)
##    def test_dtype(self):
##        ''' Test a 1d list with a new dtype'''
##        a=[10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
##        b = 34.1417152147
##        self.do(a, b, dtype=np.float128)  # does not work on Win32

class TestHarMean(HarMeanTestCase, TestCase):
    def do(self, a, b, axis=None, dtype=None):
        x = stats.hmean(a, axis=axis, dtype=dtype)
        assert_almost_equal(b, x)
	assert_equal(x.dtype, dtype)

class GeoMeanTestCase:
    def test_1dlist(self):
        ''' Test a 1d list'''
        a=[10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        b = 45.2872868812
        self.do(a, b)
    def test_1darray(self):
        ''' Test a 1d array'''
        a=np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
        b = 45.2872868812
        self.do(a, b)
    def test_1dma(self):
        ''' Test a 1d masked array'''
        a=np.ma.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
        b = 45.2872868812
        self.do(a, b)
    def test_1dmavalue(self):
        ''' Test a 1d masked array with a masked value'''
        a=np.ma.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100], mask=[0,0,0,0,0,0,0,0,0,1])
        b = 41.4716627439
        self.do(a, b)

    # Note the next tests use axis=None as default, not axis=0
    def test_2dlist(self):
        ''' Test a 2d list'''
        a=[[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = 52.8885199
        self.do(a, b)
    def test_2darray(self):
        ''' Test a 2d array'''
        a=[[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = 52.8885199
        self.do(np.array(a), b)
    def test_2dma(self):
        ''' Test a 2d masked array'''
        a=[[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = 52.8885199
        self.do(np.ma.array(a), b)
    def test_2daxis0(self):
        ''' Test a 2d list with axis=0'''
        a=[[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = np.array([35.56893304,  49.32424149,  61.3579244 ,  72.68482371])
        self.do(a, b, axis=0)
    def test_2daxis1(self):
        ''' Test a 2d list with axis=1'''
        a=[[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = np.array([  22.13363839,   64.02171746,  104.40086817])
        self.do(a, b, axis=1)
    def test_2dmatrixdaxis0(self):
        ''' Test a 2d list with axis=0'''
        a=[[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = np.matrix([[35.56893304,  49.32424149,  61.3579244 ,  72.68482371]])
        self.do(np.matrix(a), b, axis=0)
    def test_2dmatrixaxis1(self):
        ''' Test a 2d list with axis=1'''
        a=[[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
        b = np.matrix([[  22.13363839,   64.02171746,  104.40086817]]).T
        self.do(np.matrix(a), b, axis=1)
##    def test_dtype(self):
##        ''' Test a 1d list with a new dtype'''
##        a=[10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
##        b = 45.2872868812
##        self.do(a, b, dtype=np.float128)  # does not exist on win32
    def test_1dlist0(self):
        ''' Test a 1d list with zero element'''
        a=[10, 20, 30, 40, 50, 60, 70, 80, 90, 0]
        b = 0.0 # due to exp(-inf)=0
        self.do(a, b)
    def test_1darray0(self):
        ''' Test a 1d array with zero element'''
        a=np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 0])
        b = 0.0 # due to exp(-inf)=0
        self.do(a, b)
    def test_1dma0(self):
        ''' Test a 1d masked array with zero element'''
        a=np.ma.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 0])
        b = 41.4716627439
        self.do(a, b)
    def test_1dmainf(self):
        ''' Test a 1d masked array with negative element'''
        a=np.ma.array([10, 20, 30, 40, 50, 60, 70, 80, 90, -1])
        b = 41.4716627439
        self.do(a, b)

class TestGeoMean(GeoMeanTestCase, TestCase):
    def do(self, a, b, axis=None, dtype=None):
        #Note this doesn't test when axis is not specified
        x = stats.gmean(a, axis=axis, dtype=dtype)
        assert_almost_equal(b, x)
	assert_equal(x.dtype, dtype)


def test_binomtest():
    # precision tests compared to R for ticket:986
    pp = np.concatenate(( np.linspace(0.1,0.2,5), np.linspace(0.45,0.65,5),
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
                            significant=12, err_msg='fail forp=%f'%p)

    assert_approx_equal(stats.binom_test(50,100,0.1), 5.8320387857343647e-024,
                            significant=12, err_msg='fail forp=%f'%p)

class Test_Trim(object):
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
               np.array([[ 2,  8, 14, 20],[ 3,  9, 15, 21]]))
        assert_raises(ValueError, stats.trimboth,
               np.arange(24).reshape(4,6).T, 4/6.)

    def test_trim_mean(self):
        a = np.arange(11)
        assert_equal(stats.trim_mean(np.arange(24).reshape(4,6).T, 2/6.),
                        np.array([  2.5,   8.5,  14.5,  20.5]))
        assert_equal(stats.trim_mean(np.arange(24).reshape(4,6), 2/6.),
                        np.array([  9.,  10.,  11.,  12.,  13.,  14.]))
        assert_equal(stats.trim_mean(np.arange(24), 2/6.), 11.5)
        assert_equal(stats.trim_mean([5,4,3,1,2,0], 2/6.), 2.5)


class TestSigamClip(object):
    def test_sigmaclip1(self):
        a = np.concatenate((np.linspace(9.5,10.5,31),np.linspace(0,20,5)))
        fact = 4  #default
        c, low, upp = stats.sigmaclip(a)
        assert_(c.min()>low)
        assert_(c.max()<upp)
        assert_equal(low, c.mean() - fact*c.std())
        assert_equal(upp, c.mean() + fact*c.std())
        assert_equal(c.size, a.size)

    def test_sigmaclip2(self):
        a = np.concatenate((np.linspace(9.5,10.5,31),np.linspace(0,20,5)))
        fact = 1.5
        c, low, upp = stats.sigmaclip(a, fact, fact)
        assert_(c.min()>low)
        assert_(c.max()<upp)
        assert_equal(low, c.mean() - fact*c.std())
        assert_equal(upp, c.mean() + fact*c.std())
        assert_equal(c.size, 4)
        assert_equal(a.size, 36) #check original array unchanged

    def test_sigmaclip3(self):
        a = np.concatenate((np.linspace(9.5,10.5,11),np.linspace(-100,-50,3)))
        fact = 1.8
        c, low, upp = stats.sigmaclip(a, fact, fact)
        assert_(c.min()>low)
        assert_(c.max()<upp)
        assert_equal(low, c.mean() - fact*c.std())
        assert_equal(upp, c.mean() + fact*c.std())
        assert_equal(c, np.linspace(9.5,10.5,11))


if __name__ == "__main__":
    run_module_suite()
