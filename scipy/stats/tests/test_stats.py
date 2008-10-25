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
        y = stats.mean(X)
        assert_almost_equal(y, 5.0)

    def test_stdX(self):
        y = stats.std(X)
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
        y = stats.mean(ZERO)
        assert_almost_equal(y, 0.0)

    def test_stdZERO(self):
        y = stats.std(ZERO)
        assert_almost_equal(y, 0.0)

##    Really need to write these tests to handle missing values properly
##    def test_meanMISS(self):
##        y = stats.mean(MISS)
##        assert_almost_equal(y, 0.0)
##
##    def test_stdMISS(self):
##        y = stats.stdev(MISS)
##        assert_almost_equal(y, 0.0)

    def test_meanBIG(self):
        y = stats.mean(BIG)

        assert_almost_equal(y, 99999995.00)

    def test_stdBIG(self):
        y = stats.std(BIG)
        assert_almost_equal(y, 2.738612788)

    def test_meanLITTLE(self):
        y = stats.mean(LITTLE)
        assert_approx_equal(y, 0.999999950)

    def test_stdLITTLE(self):
        y = stats.std(LITTLE)
        assert_approx_equal(y, 2.738612788e-8)

    def test_meanHUGE(self):
        y = stats.mean(HUGE)
        assert_approx_equal(y, 5.00000e+12)

    def test_stdHUGE(self):
        y = stats.std(HUGE)
        assert_approx_equal(y, 2.738612788e12)

    def test_meanTINY(self):
        y = stats.mean(TINY)
        assert_almost_equal(y, 0.0)

    def test_stdTINY(self):
        y = stats.std(TINY)
        assert_almost_equal(y, 0.0)

    def test_meanROUND(self):
        y = stats.mean(ROUND)
        assert_approx_equal(y, 4.500000000)

    def test_stdROUND(self):
        y = stats.std(ROUND)
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
        assert_approx_equal(s, stats.std(self.X))

    def test_nanstd_some(self):
        """Check nanstd when some values only are nan."""
        s = stats.nanstd(self.Xsome)
        assert_approx_equal(s, stats.std(self.Xsomet))

    def test_nanstd_all(self):
        """Check nanstd when all values are nan."""
        s = stats.nanstd(self.Xall)
        assert np.isnan(s)

    def test_nanmedian_none(self):
        """Check nanmedian when no values are nan."""
        m = stats.nanmedian(self.X)
        assert_approx_equal(m, stats.median(self.X))

    def test_nanmedian_some(self):
        """Check nanmedian when some values only are nan."""
        m = stats.nanmedian(self.Xsome)
        assert_approx_equal(m, stats.median(self.Xsomet))

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
        assert_almost_equal(res[4], 4.3609875083149268e-3)

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
        assert_almost_equal(stats.mean(a),mn1,11)
        mn2 = 0.0
        for el in af:
            mn2 += el / float(Naf)
        assert_almost_equal(stats.mean(af),mn2,11)

    def test_2d(self):
        a = [[1.0, 2.0, 3.0],
             [2.0, 4.0, 6.0],
             [8.0, 12.0, 7.0]]
        A = array(a)
        N1, N2 = (3, 3)
        mn1 = zeros(N2, dtype=float)
        for k in range(N1):
            mn1 += A[k,:] / N1
        assert_almost_equal(stats.mean(a, axis=0), mn1, decimal=13)
        assert_almost_equal(stats.mean(a), mn1, decimal=13)
        mn2 = zeros(N1, dtype=float)
        for k in range(N2):
            mn2 += A[:,k]
        mn2 /= N2
        assert_almost_equal(stats.mean(a, axis=1), mn2, decimal=13)

    def test_ravel(self):
        a = rand(5,3,5)
        A = 0
        for val in ravel(a):
            A += val
        assert_almost_equal(stats.mean(a,axis=None),A/(5*3.0*5))

class TestPercentile(TestCase):
    def setUp(self):
        self.a1 = [3,4,5,10,-3,-5,6]
        self.a2 = [3,-6,-2,8,7,4,2,1]
        self.a3 = [3.,4,5,10,-3,-5,-6,7.0]

    def test_median(self):
        assert_equal(stats.median(self.a1), 4)
        assert_equal(stats.median(self.a2), 2.5)
        assert_equal(stats.median(self.a3), 3.5)

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
        assert_almost_equal(stats.std(a),5.2098807225172772,11)
        assert_almost_equal(stats.std(b),5.9281411203561225,11)

    def test_2d(self):
        a = [[1.0, 2.0, 3.0],
             [2.0, 4.0, 6.0],
             [8.0, 12.0, 7.0]]
        b1 = array((3.7859388972001824, 5.2915026221291814,
                    2.0816659994661335))
        b2 = array((1.0,2.0,2.64575131106))
        assert_array_almost_equal(stats.std(a),b1,11)
        assert_array_almost_equal(stats.std(a,axis=0),b1,11)
        assert_array_almost_equal(stats.std(a,axis=1),b2,11)


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
        assert_almost_equal(stats.median(data1),2.5)
        assert_almost_equal(stats.median(data2),5)

    def test_basic2(self):
        a1 = [3,4,5,10,-3,-5,6]
        a2 = [3,-6,-2,8,7,4,2,1]
        a3 = [3.,4,5,10,-3,-5,-6,7.0]
        assert_equal(stats.median(a1),4)
        assert_equal(stats.median(a2),2.5)
        assert_equal(stats.median(a3),3.5)

    def test_axis(self):
        """Regression test for #760."""
        a1 = np.array([[3,4,5], [10,-3,-5]])
        assert_equal(stats.median(a1), np.array([6.5, 0.5, 0.]))
        assert_equal(stats.median(a1, axis=-1), np.array([4., -3]))

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
        y = stats.std(self.testcase)
        assert_approx_equal(y,1.290994449)

    def test_var(self):
        """
        var(testcase) = 1.666666667 """
        #y = stats.var(self.shoes[0])
        #assert_approx_equal(y,6.009)
        y = stats.var(self.testcase)
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
        y = stats.z(self.testcase,stats.mean(self.testcase))
        assert_almost_equal(y,0.0)

    def test_zs(self):
        """
        not in R, so tested by using
        (testcase[i]-mean(testcase,axis=0))/sqrt(var(testcase)*3/4)
        """
        y = stats.zs(self.testcase)
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
            sum((testmathworks-mean(testmathworks,axis=0))**3,axis=0)/((sqrt(var(testmathworks)*4/5))**3)/5
        """
        y = stats.skew(self.testmathworks)
        assert_approx_equal(y,-0.29322304336607,10)
        y = stats.skew(self.testmathworks,bias=0)
        assert_approx_equal(y,-0.437111105023940,10)
        y = stats.skew(self.testcase)
        assert_approx_equal(y,0.0,10)
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

if __name__ == "__main__":
    run_module_suite()
