""" Test functions for stats module

    WRITTEN BY LOUIS LUANGKESORN <lluang@yahoo.com> FOR THE STATS MODULE
    BASED ON WILKINSON'S STATISTICS QUIZ
    http://www.stanford.edu/~clint/bench/wilk.txt

"""

import unittest
import sys
from scipy_test.testing import *
set_package_path()
#XXX: fix scipy.stats to stats
import scipy
import scipy_base.fastumath as math
del sys.path[0]

import Numeric
N = Numeric


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
##  For completeness, I should write the relavant tests and count them as failures
##  Somewhat acceptable, since this is still beta software.  It would count as a
##  good target for 1.0 status
X = N.array([1,2,3,4,5,6,7,8,9],N.Float)
ZERO= N.array([0,0,0,0,0,0,0,0,0], N.Float)
#MISS=N.array([.,.,.,.,.,.,.,.,.], N.Float)
BIG=N.array([99999991,99999992,99999993,99999994,99999995,99999996,99999997,99999998,99999999],N.Float)
LITTLE=N.array([0.99999991,0.99999992,0.99999993,0.99999994,0.99999995,0.99999996,0.99999997,0.99999998,0.99999999],N.Float)
HUGE=N.array([1e+12,2e+12,3e+12,4e+12,5e+12,6e+12,7e+12,8e+12,9e+12],N.Float)
TINY=N.array([1e-12,2e-12,3e-12,4e-12,5e-12,6e-12,7e-12,8e-12,9e-12],N.Float)
ROUND=N.array([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5],N.Float)
X2 = X * X
X3 = X2 * X
X4 = X3 * X
X5 = X4 * X
X6 = X5 * X
X7 = X6 * X
X8 = X7 * X
X9 = X8 * X

class test_round(unittest.TestCase):
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

    def check_rounding0(self):
        """ W.II.A.0. Print ROUND with only one digit.
        
            You should get the numbers 1 to 9.  Many language compilers,
            such as Turbo Pascal and Lattice C, fail this test (they round
            numbers inconsistently). Needless to say, statical packages
            written in these languages may fail the test as well.
        """
 
        for i in range(0,9):
            y = scipy.round(ROUND[i])
            assert_equal(y,i+1)
        
    def check_rounding1(self):
        """ W.II.A.1. Y = INT(2.6*7 -0.2) (Y should be 18)"""
        y = int(2.6*7 -0.2)
        assert_equal(y, 18)
        
    def check_rounding2(self):
        """ W.II.A.2. Y = 2-INT(EXP(LOG(SQR(2)*SQR(2))))   (Y should be 0)"""
        y=2-int(scipy.exp(scipy.log(scipy.sqrt(2.)*scipy.sqrt(2.))))
        assert_equal(y,0)
        
    def check_rounding3(self):
        """ W.II.A.3. Y = INT(3-EXP(LOG(SQR(2)*SQR(2))))    (Y should be 1)"""
        y=(int(round((3-scipy.exp(scipy.log(scipy.sqrt(2.0)*scipy.sqrt(2.0)))))))
        assert_equal(y,1)

class test_basicstats(unittest.TestCase):
    """ W.II.C. Compute basic statistic on all the variables.

        The means should be the fifth value of all the variables (case FIVE).
        The standard deviations should be "undefined" or missing for MISS,
        0 for ZERO, and 2.738612788 (times 10 to a power) for all the other variables.
        II. C. Basic Statistics
    """

    def check_meanX(self):
        y = scipy.stats.mean(X)
        assert_almost_equal(y,5.0)
    def check_stdX(self):
        y = scipy.stats.std(X)
        assert_almost_equal(y,2.738612788)

    def check_meanZERO(self):
        y = scipy.stats.mean(ZERO)
        assert_almost_equal(y,0.0)

    def check_stdZERO(self):
        y = scipy.stats.std(ZERO)
        assert_almost_equal(y,0.0)

##    Really need to write these tests to handle missing values properly
##    def check_meanMISS(self):
##        y = scipy.stats.mean(MISS)
##        assert_almost_equal(y,0.0)
##
##    def check_stdMISS(self):
##        y = scipy.stats.stdev(MISS)
##        assert_almost_equal(y,0.0)

    def check_meanBIG(self):
        y = scipy.stats.mean(BIG)
        assert_almost_equal(y,99999995.00)

    def check_stdBIG(self):
        y = scipy.stats.std(BIG)
        assert_almost_equal(y,2.738612788)

    def check_meanLITTLE(self):
        y = scipy.stats.mean(LITTLE)
        assert_approx_equal(y,0.999999950)

    def check_stdLITTLE(self):
        y = scipy.stats.std(LITTLE)
        assert_approx_equal(y,2.738612788e-8)

    def check_meanHUGE(self):
        y = scipy.stats.mean(HUGE)
        assert_approx_equal(y,5.00000e+12)

    def check_stdHUGE(self):
        y = scipy.stats.std(HUGE)
        assert_approx_equal(y,2.738612788e12)

    def check_meanTINY(self):
        y = scipy.stats.mean(TINY)
        assert_almost_equal(y,0.0)

    def check_stdTINY(self):
        y = scipy.stats.std(TINY)
        assert_almost_equal(y,0.0)

    def check_meanROUND(self):
        y = scipy.stats.mean(ROUND)
        assert_approx_equal(y,4.500000000)

    def check_stdROUND(self):
        y = scipy.stats.std(ROUND)
        assert_approx_equal(y,2.738612788)
            
class test_corr(unittest.TestCase):
    """ W.II.D. Compute a correlation matrix on all the variables.

        All the correlations, except for ZERO and MISS, shoud be exactly 1.
        ZERO and MISS should have undefined or missing correlations with the
        other variables.  The same should go for SPEARMAN corelations, if
        your program has them.
    """
    def check_pXX(self):    
        y = scipy.stats.pearsonr(X,X)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pXBIG(self):    
        y = scipy.stats.pearsonr(X,BIG)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pXLITTLE(self):    
        y = scipy.stats.pearsonr(X,LITTLE)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pXHUGE(self):    
        y = scipy.stats.pearsonr(X,HUGE)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pXTINY(self):    
        y = scipy.stats.pearsonr(X,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pXROUND(self):    
        y = scipy.stats.pearsonr(X,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pBIGBIG(self):    
        y = scipy.stats.pearsonr(BIG,BIG)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pBIGLITTLE(self):    
        y = scipy.stats.pearsonr(BIG,LITTLE)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pBIGHUGE(self):    
        y = scipy.stats.pearsonr(BIG,HUGE)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pBIGTINY(self):    
        y = scipy.stats.pearsonr(BIG,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pBIGROUND(self):    
        y = scipy.stats.pearsonr(BIG,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pLITTLELITTLE(self):    
        y = scipy.stats.pearsonr(LITTLE,LITTLE)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pLITTLEHUGE(self):    
        y = scipy.stats.pearsonr(LITTLE,HUGE)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pLITTLETINY(self):    
        y = scipy.stats.pearsonr(LITTLE,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pLITTLEROUND(self):    
        y = scipy.stats.pearsonr(LITTLE,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pHUGEHUGE(self):    
        y = scipy.stats.pearsonr(HUGE,HUGE)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pHUGETINY(self):    
        y = scipy.stats.pearsonr(HUGE,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pHUGEROUND(self):    
        y = scipy.stats.pearsonr(HUGE,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pTINYTINY(self):    
        y = scipy.stats.pearsonr(TINY,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pTINYROUND(self):    
        y = scipy.stats.pearsonr(TINY,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_pROUNDROUND(self):    
        y = scipy.stats.pearsonr(ROUND,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sXX(self):    
        y = scipy.stats.spearmanr(X,X)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sXBIG(self):    
        y = scipy.stats.spearmanr(X,BIG)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sXLITTLE(self):    
        y = scipy.stats.spearmanr(X,LITTLE)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sXHUGE(self):    
        y = scipy.stats.spearmanr(X,HUGE)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sXTINY(self):    
        y = scipy.stats.spearmanr(X,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sXROUND(self):    
        y = scipy.stats.spearmanr(X,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sBIGBIG(self):    
        y = scipy.stats.spearmanr(BIG,BIG)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sBIGLITTLE(self):    
        y = scipy.stats.spearmanr(BIG,LITTLE)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sBIGHUGE(self):    
        y = scipy.stats.spearmanr(BIG,HUGE)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sBIGTINY(self):    
        y = scipy.stats.spearmanr(BIG,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sBIGROUND(self):    
        y = scipy.stats.spearmanr(BIG,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sLITTLELITTLE(self):    
        y = scipy.stats.spearmanr(LITTLE,LITTLE)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sLITTLEHUGE(self):    
        y = scipy.stats.spearmanr(LITTLE,HUGE)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sLITTLETINY(self):    
        y = scipy.stats.spearmanr(LITTLE,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sLITTLEROUND(self):    
        y = scipy.stats.spearmanr(LITTLE,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sHUGEHUGE(self):    
        y = scipy.stats.spearmanr(HUGE,HUGE)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sHUGETINY(self):    
        y = scipy.stats.spearmanr(HUGE,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sHUGEROUND(self):    
        y = scipy.stats.spearmanr(HUGE,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sTINYTINY(self):    
        y = scipy.stats.spearmanr(TINY,TINY)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sTINYROUND(self):    
        y = scipy.stats.spearmanr(TINY,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)
    def check_sROUNDROUND(self):    
        y = scipy.stats.spearmanr(ROUND,ROUND)
        r = y[0]
        assert_approx_equal(r,1.0)

##    W.II.E.  Tabulate X against X, using BIG as a case weight.  The values 
##    should appear on the diagonal and the total should be 899999955.
##    If the table cannot hold these values, forget about working with 
##    census data.  You can also tabulate HUGE against TINY.  There is no
##    reason a tabulation program should not be able to digtinguish 
##    different values regardless of their magnitude.

### I need to figure out how to do this one.


class test_regression(unittest.TestCase):
    def check_linregressBIGX(self):
        """ W.II.F.  Regress BIG on X.

            The constant should be 99999990 and the regression coefficient should be 1.
        """
        y = scipy.stats.linregress(X,BIG)
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


    def check_regressXX(self):
        """ W.IV.B.  Regress X on X.

            The constant should be exactly 0 and the regression coefficient should be 1.
            This is a perfectly valid regression.  The program should not complain.
        """
        y = scipy.stats.linregress(X,X)
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

    def check_regressZEROX(self):
        """ W.IV.D. Regress ZERO on X.

            The program should inform you that ZERO has no variance or it should
            go ahead and compute the regression and report a correlation and
            total sum of squares of exactly 0.
        """
        y = scipy.stats.linregress(X,ZERO)
        intercept = y[1]
        r=y[2]
        assert_almost_equal(intercept,0.0)
        assert_almost_equal(r,0.0)

anovadata = ([[1,'A1','B1',2],
              [2,'A1','B1',1],
              [3,'A1','B1',3],
              [4,'A1','B1',2],
              [5,'A1','B2',3],
              [6,'A1','B2',4],
              [7,'A1','B2',5],
              [8,'A2','B1',4],
              [9,'A2','B1',6],
              [10,'A2','B1',5],
              [11,'A1','B2',2],
              [12,'A1','B2',4],
              [13,'A1','B2',5]])

##class test_anova(unittest.TestCase):
##    """ W.V.A.  Simple contrasts.
##
##        The following data contain an unbalanced design with a significant
##        interaction.  A least squares analysis shows the main efefct for A
##        is not significant, but this test is not particularlrly meaningful
##        because of the interaction.  Test, therefore, the simple constrast
##        between A1 and A2 within B1.  Then test A1 vs. A2 within B2.  Both
##        tests should use the same residual error term (separate t-tests are 
##        unacceptable).  Several widely used mainframe programs fail this
##        test.  Unless the program can constrast any terms in a model (not
##        just main effects), it cannot handle this frequently encountered
##        type of problem.
##                    B1        B2
##               ---------------------
##               |    2     |   3    |
##        A1     |    1     |   4    |
##               |    3     |   5    |
##               |    2     |        |
##               ---------------------
##               |    4     |   2    |
##        A2     |    6     |   4    |
##               |    5     |   5    |
##               ---------------------
##        The numbers in Wilkinson's statistic quiz are from SYSTAT output,
##        not certified to be the correct values.
##    """
##    y=scipy.stats.anova(anovadata)
##    for i in [0,len(y)]:
##        if (y[i][0]=='A'):
##            assert_approx_equal(y[i][1],5.689,3)
##            """ Sum Squares - A"""
##            assert_equal(y[i][2],1)
##            """ Degrees Freedom - A """
##            assert_approx_equal(y[i][3],5.689,3)
##            """ Mean Squares - A"""
##            assert_approx_equal(y[i][4],4.800,3)
##            """ F ratio - A"""
##        if (y[i][0]=='B'):
##            assert_approx_equal(y[i][1],0.356,2)
##            """ Sum Squares - B """
##            assert_equal(y[i][2],1)
##            """ Degrees Freedom - B """
##            assert_approx_equal(y[i][3],0.356,2)
##            """ Mean Squares - B """
##            assert_approx_equal(y[i][4],0.300,2)
##            """ F ratio - B """
##        if (y[i][0]=='AB'):
##            assert_approx_equal(y[i][1],8.889,3)
##            """ Sum Squares - AB """
##            assert_equal(y[i][2],1)
##            """ Degrees Freedom - AB """
##            assert_approx_equal(y[i][3],8.889,3)
##            """ Mean Squares - AB """
##            assert_approx_equal(y[i][4],7.500,3)
##            """ F ratio - AB """

# Utility

def compare_results(res,desired):
    for i in range(len(desired)):
        assert_array_equal(res[i],desired[i])


from Numeric import *
import Numeric
from scipy_base.fastumath import *
import unittest
from scipy_test.testing import assert_array_equal, assert_equal
from scipy_test.testing import assert_almost_equal
from scipy_test.testing import assert_array_almost_equal, rand
from scipy_test.testing import ScipyTestCase as TestCase
import scipy
from scipy import stats


##################################################
### Test for sum

class test_gmean(TestCase):

    def check_1D_list(self):
        a = (1,2,3,4)
        actual= stats.gmean(a)
        desired = power(1*2*3*4,1./4.)
        assert_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.gmean(a,axis=-1)
        assert_almost_equal(desired1,actual,decimal=14)
    def check_1D_array(self):
        a = array((1,2,3,4),Float32)
        actual= stats.gmean(a)
        desired = power(1*2*3*4,1./4.)
        assert_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.gmean(a,axis=-1)
        assert_almost_equal(desired1,actual,decimal=14)

    def check_2D_array_default(self):
        a = array(((1,2,3,4),
                   (1,2,3,4),
                   (1,2,3,4)))
        actual= stats.gmean(a)
        v = power(1*2*3*4,1./4.)
        desired = array((v,v,v))
        assert_array_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.gmean(a,axis=-1)
        assert_array_almost_equal(desired1,actual,decimal=14)

    def check_2D_array_dim0(self):
        a = array(((1,2,3,4),
                   (1,2,3,4),
                   (1,2,3,4)))
        actual= stats.gmean(a,axis=0)
        desired = array((1,2,3,4))
        assert_array_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.gmean(a,axis=0)
        assert_array_almost_equal(desired1,actual,decimal=14)

class test_hmean(TestCase):
    def check_1D_list(self):
        a = (1,2,3,4)
        actual= stats.hmean(a)
        desired =  4. / (1./1 + 1./2 + 1./3 + 1./4)
        assert_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.hmean(array(a),axis=-1)
        assert_almost_equal(desired1,actual,decimal=14)
    def check_1D_array(self):
        a = array((1,2,3,4),Float64)
        actual= stats.hmean(a)
        desired =  4. / (1./1 + 1./2 + 1./3 + 1./4)
        assert_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.hmean(a,axis=-1)
        assert_almost_equal(desired1,actual,decimal=14)

    def check_2D_array_default(self):
        a = array(((1,2,3,4),
                   (1,2,3,4),
                   (1,2,3,4)))
        actual= stats.hmean(a)
        v = 4. / (1./1 + 1./2 + 1./3 + 1./4)
        desired = array((v,v,v))      
        assert_array_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.hmean(a,axis=-1)
        assert_array_almost_equal(desired1,actual,decimal=14)

    def check_2D_array_dim0(self):
        a = array(((1,2,3,4),
                   (1,2,3,4),
                   (1,2,3,4)))
        actual= stats.hmean(a,axis=0)
        desired = array((1.,2.,3.,4.))        
        assert_array_almost_equal(desired,actual,decimal=14)

        desired1 = scipy.stats.stats.hmean(a,axis=0)
        assert_array_almost_equal(desired1,actual,decimal=14)

class test_mean(TestCase):
    def check_basic(self):
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

    def check_2d(self):
        a = [[1.0, 2.0, 3.0],
             [2.0, 4.0, 6.0],
             [8.0, 12.0, 7.0]]
        A = array(a,'d')
        N1,N2 = (3,3)
        mn1 = zeros(N2,'d')
        for k in range(N1):
            mn1 += A[k,:] / N1
        Numeric.allclose(stats.mean(a),mn1,rtol=1e-13,atol=1e-13)
        mn2 = zeros(N1,'d')
        for k in range(N2):
            mn2 += A[:,k] / N2
        Numeric.allclose(stats.mean(a,axis=0),mn2,rtol=1e-13,atol=1e-13)

    def check_ravel(self):
        a = rand(5,3,5)
        A = 0
        for val in ravel(a):
            A += val
        assert_almost_equal(stats.mean(a,axis=None),A/(5*3.0*5))

class test_median(TestCase):
    def check_basic(self):
        a1 = [3,4,5,10,-3,-5,6]
        a2 = [3,-6,-2,8,7,4,2,1]
        a3 = [3.,4,5,10,-3,-5,-6,7.0]
        assert_equal(stats.median(a1),4)
        assert_equal(stats.median(a2),2.5)
        assert_equal(stats.median(a3),3.5)        

class test_std(TestCase):
    def check_basic(self):
        a = [3,4,5,10,-3,-5,6]
        af = [3.,4,5,10,-3,-5,-6]
        Na = len(a)
        Naf = len(af)
        assert_almost_equal(stats.std(a),5.2098807225172772,11)
        assert_almost_equal(stats.std(af),5.9281411203561225,11)

    def check_2d(self):
        a = [[1.0, 2.0, 3.0],
             [2.0, 4.0, 6.0],
             [8.0, 12.0, 7.0]]
        b1 = array((3.7859388972001824, 5.2915026221291814,
                    2.0816659994661335))
        b2 = array((1.0,2.0,2.64575131106))
        assert_array_almost_equal(stats.std(a),b2,11)
        assert_array_almost_equal(stats.std(a,axis=0),b1,11)


class test_cmedian(TestCase):
    def check_basic(self):
        data = [1,2,3,1,5,3,6,4,3,2,4,3,5,2.0]
        assert_almost_equal(stats.cmedian(data,5),2.7)
        assert_almost_equal(stats.cmedian(data,3),3.1296296296296298)
        assert_almost_equal(stats.cmedian(data),2.999665)

class test_median(TestCase):
    def check_basic(self):
        data1 = [1,3,5,2,3,1,19,-10,2,4.0]
        data2 = [3,5,1,10,23,-10,3,-2,6,8,15]
        assert_almost_equal(stats.median(data1),2.5)
        assert_almost_equal(stats.median(data2),5)

class test_mode(TestCase):
    def check_basic(self):
        data1 = [3,5,1,10,23,3,2,6,8,6,10,6]
        vals = stats.mode(data1)
        assert_almost_equal(vals[0][0],6)
        assert_almost_equal(vals[1][0],3)


class test_variability(TestCase):
    """  Comparison numbers are found using R v.1.5.1
         note that length(testcase) = 4
    """
    testcase = [1,2,3,4]
    def check_std(self):
        y = scipy.stats.std(self.testcase)
        assert_approx_equal(y,1.290994449)

    def check_var(self):
        """
        var(testcase) = 1.666666667 """
        #y = scipy.stats.var(self.shoes[0])
        #assert_approx_equal(y,6.009)
        y = scipy.stats.var(self.testcase)
        assert_approx_equal(y,1.666666667)

    def check_samplevar(self):
        """
        R does not have 'samplevar' so the following was used
        var(testcase)*(4-1)/4  where 4 = length(testcase)
        """
        #y = scipy.stats.samplevar(self.shoes[0])
        #assert_approx_equal(y,5.4081)
        y = scipy.stats.samplevar(self.testcase)
        assert_approx_equal(y,1.25)
        
    def check_samplestd(self):
        #y = scipy.stats.samplestd(self.shoes[0])
        #assert_approx_equal(y,2.325532197)
        y = scipy.stats.samplestd(self.testcase)
        assert_approx_equal(y,1.118033989)

    def check_signaltonoise(self):
        """
        this is not in R, so used
        mean(testcase)/(sqrt(var(testcase)*3/4)) """
        #y = scipy.stats.signaltonoise(self.shoes[0])
        #assert_approx_equal(y,4.5709967)
        y = scipy.stats.signaltonoise(self.testcase)
        assert_approx_equal(y,2.236067977)

    def check_stderr(self):
        """
        this is not in R, so used
        sqrt(var(testcase))/sqrt(4)
        """
##        y = scipy.stats.stderr(self.shoes[0])
##        assert_approx_equal(y,0.775177399)
        y = scipy.stats.stderr(self.testcase)
        assert_approx_equal(y,0.6454972244)
    def check_sem(self):
        """
        this is not in R, so used
        sqrt(var(testcase)*3/4)/sqrt(3)
        """
        #y = scipy.stats.sem(self.shoes[0])
        #assert_approx_equal(y,0.775177399)
        y = scipy.stats.sem(self.testcase)
        assert_approx_equal(y,0.6454972244)
        
    def check_z(self):
        """ 
        not in R, so used
        (10-mean(testcase))/sqrt(var(testcase)*3/4)
        """
        y = scipy.stats.z(self.testcase,scipy.stats.mean(self.testcase))
        assert_almost_equal(y,0.0)

    def check_zs(self):
        """
        not in R, so tested by using
        (testcase[i]-mean(testcase))/sqrt(var(testcase)*3/4)
        """
        y = scipy.stats.zs(self.testcase)
        desired = ([-1.3416407864999, -0.44721359549996 , 0.44721359549996 , 1.3416407864999])
        assert_array_almost_equal(desired,y,decimal=12)
        


class test_moments(TestCase):
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
    def check_moment(self):
        """
        mean((testcase-mean(testcase))**power))"""
        y = scipy.stats.moment(self.testcase,1)
        assert_approx_equal(y,0.0,10)
        y = scipy.stats.moment(self.testcase,2)
        assert_approx_equal(y,1.25)
        y = scipy.stats.moment(self.testcase,3)
        assert_approx_equal(y,0.0)
        y = scipy.stats.moment(self.testcase,4)
        assert_approx_equal(y,2.5625)
    def check_variation(self):
        """
        variation = samplestd/mean """
##        y = scipy.stats.variation(self.shoes[0])
##        assert_approx_equal(y,21.8770668)
        y = scipy.stats.variation(self.testcase)
        assert_approx_equal(y,0.44721359549996, 10)

    def check_skewness(self):
        """
            sum((testmathworks-mean(testmathworks))**3)/((sqrt(var(testmathworks)*4/5))**3)/5
        """
        y = scipy.stats.skew(self.testmathworks)
        assert_approx_equal(y,-0.29322304336607,10)
        y = scipy.stats.skew(self.testmathworks,bias=0)
        assert_approx_equal(y,-0.437111105023940,10)
        y = scipy.stats.skew(self.testcase)
        assert_approx_equal(y,0.0,10)
    def check_kurtosis(self):
        """
            sum((testcase-mean(testcase))**4)/((sqrt(var(testcase)*3/4))**4)/4
            sum((test2-mean(testmathworks))**4)/((sqrt(var(testmathworks)*4/5))**4)/5
            Set flags for axis = 0 and
            fisher=0 (Pearson's defn of kurtosis for compatiability with Matlab)
        """
        y = scipy.stats.kurtosis(self.testmathworks,0,fisher=0,bias=1)
        assert_approx_equal(y, 2.1658856802973,10)

        # Note that MATLAB has confusing docs for the following case
        #  kurtosis(x,0) gives an unbiased estimate of Pearson's skewness
        #  kurtosis(x)  gives a biased estimate of Fisher's skewness (Pearson-3)
        #  The MATLAB docs imply that both should give Fisher's
        y = scipy.stats.kurtosis(self.testmathworks,fisher=0,bias=0)
        assert_approx_equal(y, 3.663542721189047,10)
        y = scipy.stats.kurtosis(self.testcase,0,0)
        assert_approx_equal(y,1.64)

if __name__ == "__main__":
    ScipyTest('stats.stats').run()


