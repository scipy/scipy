""" Test functions for stats module

    THIS IS BEING REWRITTEN BY LOUIS LUANGKESORN <lluang@yahoo.com> FOR THE STATS MODULE
    BASED ON WILKINSON'S STATISTICS QUIZ
"""

import unittest
import scipy
from scipy.scipy_test import assert_array_equal, assert_equal
from scipy.scipy_test import assert_almost_equal, assert_array_almost_equal
#from scipy import *
from scipy import stats
import math
import Numeric
N = Numeric

##  Numbers in comments refer to the section numbers and headings
##  found in the STATISTICS QUIZ of Leland Wilkinson.  These are
##  considered to be essential functionality.  True testing and
##  evaluation of a statistics package requires use of the
##  NIST Statistical test data.  See McCoullough(1999) Assessing The Reliability
##  of Statistical Software for a test methodology and its
##  implementation in testing SAS, SPSS, and S-Plus

##################################################
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
##X = [1,2,3,4,5,6,7,8,9]
##ZERO= [0,0,0,0,0,0,0,0,0]
###MISS=[.,.,.,.,.,.,.,.,.]
##BIG=[99999991,99999992,99999993,99999994,99999995,99999996,99999997,99999998,99999999]
##LITTLE=[0.99999991,0.99999992,0.99999993,0.99999994,0.99999995,0.99999996,0.99999997,0.99999998,0.99999999]
##HUGE=[1e+12,2e+12,3e+12,4e+12,5e+12,6e+12,7e+12,8e+12,9e+12]
##TINY=[1e-12,2e-12,3e-12,4e-12,5e-12,6e-12,7e-12,8e-12,9e-12]
##ROUND=[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5]

##  II.A  Round
##         A. Print ROUND with only one digit.  You should get the 
##    numbers 1 to 9.  Many language compilers, such as Turbo Pascal and
##    Lattice C, fail this test (they round numbers inconsistently).
##    Needless to say, statical packages written in these languages
##    may fail the test as well.  You can also check the following
##    expressions:
##     Y = INT(2.6*7 -0.2)                   (Y should be 18)
##     Y = 2-INT(EXP(LOG(SQR(2)*SQR(2))))    (Y should be 0)
##     Y = INT(3-EXP(LOG(SQR(2)*SQR(2))))    (Y should be 1)
##
##    INT is the integer function.  It converts decimal numbers to
##    integers by throwing away numbers after the decimal point.  EXP
##    is exponential, LOG is logarithm, and SQR is suqare root.  You may
##    have to substitute similar names for these functions for different
##    packages.  Since the square of a square root should return the same
##    number, and the exponential of a log should return the same number,
##    we should get back a 2 from this function of functions.  By taking
##    the integer result and subtracting from 2, we are exposing the 
##    roundoff errors.  These simple functions are at the heart of  
##    statistical calculations.
class test_round(unittest.TestCase):
    def check_rounding0(self):
        for i in range(0,9):
            y = scipy.round(ROUND[i])
            assert_almost_equal(y,i+1)
        
    def check_rounding1(self):
        y = math.floor(2.6*7 -0.2)
        assert_almost_equal(y, 18.0)
        
    def check_rounding2(self):
        y=2-math.floor(math.exp(math.log(math.sqrt(2)*math.sqrt(2))))
        assert_almost_equal(y,0.0)
        
    def check_rounding3(self):
        y=(math.floor(3-math.exp(math.log(math.sqrt(2.0)*math.sqrt(2.0)))))
        assert_almost_equal(y,1.00000000)
##    II.C. Compute basic statistic on all the variables.  The means should 
##    be the fifth value of all the variables (case FIVE).  The standard
##    deviations should be "undefined" or missing for MISS, 0 for ZERO,
##    and 2.738612788 (times 10 to a power) for all the other variables.
##    II. C. Basic Statistics
##
##    variable        mean            standard deviation
##           X        5.000000000     2.738612788
##        ZERO        0.000000000     0.000000000
##        MISS        .               .
##         BIG        99999995.00     2.738612788
##      LITTLE        0.999999950     2.73861E-08
##        HUGE        5.00000E+12     0.22739E+13
##        TINY        0.000000000     0.000000000
##       ROUND        4.500000000     2.738612788

class test_basicstats(unittest.TestCase):
    def check_meanX(self):
        y = scipy.stats.mean(X)
        assert_almost_equal(y,5.0)
    def check_stdX(self):
        y = scipy.stats.stdev(X)
        assert_almost_equal(y,2.738612788)

    def check_meanZERO(self):
        y = scipy.stats.mean(ZERO)
        assert_almost_equal(y,0.0)

    def check_stdZERO(self):
        y = scipy.stats.stdev(ZERO)
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
        y = scipy.stats.stdev(BIG)
        assert_almost_equal(y,2.738612788)

    def check_meanLITTLE(self):
        y = scipy.stats.mean(LITTLE)
        assert_almost_equal(y,0.999999950)

    def check_stdLITTLE(self):
        y = scipy.stats.stdev(LITTLE)
        assert_almost_equal(y,2.738612788e-08)

    def check_meanHUGE(self):
        y = scipy.stats.mean(HUGE)
        assert_almost_equal(y,5.00000e+12)

    def check_stdHUGE(self):
        y = scipy.stats.stdev(HUGE)
        assert_almost_equal(y,2.738612788e+12)

    def check_meanTINY(self):
        y = scipy.stats.mean(TINY)
        assert_almost_equal(y,0.0)

    def check_stdTINY(self):
        y = scipy.stats.stdev(TINY)
        assert_almost_equal(y,0.0)

    def check_meanROUND(self):
        y = scipy.stats.mean(ROUND)
        assert_almost_equal(y,4.500000000)

    def check_stdROUND(self):
        y = scipy.stats.stdev(ROUND)
        assert_almost_equal(y,2.738612788)
            
##    II.D. Compute a correlation matrix on all the variables.  All the
##    correlations, except for ZERO and MISS, shoud be exactly 1.  ZERO
##    and MISS should have undefined or missing correlations with the
##    other variables.  The same should go for SPEARMAN corelations, if
##    your program has them.

testcorr = N.array([[ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
       [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
       [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
       [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
       [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
       [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
       [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
       [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.],
       [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.]])    

class test_corr(unittest.TestCase):

##    check_ZERO should end up with UNDEFINED    
##    def check_ZERO(self):
##        Z = N.array((X,ZERO))
##        y = scipy.stats.acorrelation(Z)
##        assert_array_almost_equal(y,testcorr):
    
##    check_MISS should result with undefined
##    def check_MISS(self):
##        Z = N.array((X,MISS))
##        y = scipy.stats.acorrelation(Z)
##        assert_array_almost_equal(y,testcorr):

    def check_BIG(self):
        Z = N.array((X,BIG))
        y = scipy.stats.acorrelation(Z)
        assert_array_almost_equal(y,testcorr)
  
    def check_LITTLE(self):
        Z = N.array((X,LITTLE))
        y = scipy.stats.acorrelation(Z)
        assert_array_almost_equal(y,testcorr)

    def check_HUGE(self):
        Z = N.array((X,HUGE))
        y = scipy.stats.acorrelation(Z)
        assert_array_almost_equal(y,testcorr)

    def check_TINY(self):
        Z = N.array((X,TINY))
        y = scipy.stats.acorrelation(Z)
        assert_array_almost_equal(y,testcorr)

    def check_ROUND(self):
        Z = N.array((X,ROUND))
        y = scipy.stats.acorrelation(Z)
        assert_array_almost_equal(y,testcorr)


##    II.E.  Tabulate X against X, using BIG as a case weight.  The values 
##    should appear on the diagonal and the total should be 899999955.
##    If the table cannot hold these values, forget about working with 
##    census data.  You can also tabulate HUGE against TINY.  There is no
##    reason a tabulation program should not be able to digtinguish 
##    different values regardless of their magnitude.

### I need to figure out how to do this one.
X1 = X
X2 = N.matrixmultiply(X1,X)
X3 = N.matrixmultiply(X2,X)
X4 = N.matrixmultiply(X3,X)
X5 = N.matrixmultiply(X4,X)
X6 = N.matrixmultiply(X5,X)
X7 = N.matrixmultiply(X6,X)
X8 = N.matrixmultiply(X7,X)
X9 = N.matrixmultiply(X8,X)

##    II.F.  Regress BIG on X.  The constant should be 99999990 and the
##    regression coefficient should be 1.

class test_regression(unittest.TestCase):
    def check_linregressBIGX(self):
        y = scipy.stats.linregress(BIG,X)
        intercept = y[1]
        r=y[2]
        assert_almost_equal(intercept,99999990)
        assert_almost_equal(r,1.0)
        
##     IV.A. Take the NASTY dataset above.  Use the variable X as a
##     basis for computing polynomials.  Namely, compute X1=X, X2=X*X,
##     X3=X*X*X, and so on up to 9 products.  Use the algebraic
##     transformation language within the statistical package itself.  You
##     will end up with 9 variables.  Now regress X1 on X2-X9 (a perfect
##     fit).  If the package balks (singular or roundoff error messages),
##     try X1 on X2-X8, and so on.  Most packages cannot handle more than
##     a few polynomials.

##     B.  Regress X on X.  The constant should be exactly 0 and the
##     regression coefficient should be 1.  This is a perfectly valid
##     regression.  The program should not complain.
    def check_regressXX(self):
        y = scipy.stats.linregress(X,X)
        intercept = y[1]
        r=y[2]
        assert_almost_equal(intercept,0.0)
        assert_almost_equal(r,1.0)
##     C. Regress X on BIG and LITTLE (two predictors).  The program
##     should tell you that this model is "singular" because BIG and
##     LITTLE are linear combinations of each other.  Cryptic error
##     messages are unacceptable here.  Singularity is the most
##     fundamental regression error.
### Need to figure out how to handle multiple linear regression.  Not obvious

##     D. Regress ZERO on X.  The program should inform you that ZERO has
##     no variance or it should go ahead and compute the regression
##     and report a correlation and total sum of squares of exactly 0.
    def check_regressZEROX(self):
        y = scipy.stats.linregress(ZERO,X)
        intercept = y[1]
        r=y[2]
        assert_almost_equal(intercept,0.0)
        assert_almost_equal(r,0.0)

##anovadata = ([[1,'A1','B1',2],
##              [2,'A1','B1',1],
##              [3,'A1','B1',3],
##              [4,'A1','B1',2],
##              [5,'A1','B2',3],
##              [6,'A1','B2',4],
##              [7,'A1','B2',5],
##              [8,'A2','B1',4],
##              [9,'A2','B1',6],
##              [10,'A2','B1',5],
##              [11,'A1','B2',2],
##              [12,'A1','B2',4],
##              [13,'A1','B2',5]])
##
##class test_anova(unittest.TestCase):
##    y=scipy.stats.anova(anovadata)
    
    
        


# Utility

def compare_results(res,desired):
    for i in range(len(desired)):
        assert_array_equal(res[i],desired[i])

##################################################


def test_suite():
    suites = []
    suites.append( unittest.makeSuite(test_round, 'check_') )
    suites.append( unittest.makeSuite(test_basicstats, 'check_') )
    suites.append( unittest.makeSuite(test_corr, 'check_') )
    suites.append( unittest.makeSuite(test_regression, 'check_') )
##    suites.append( unittest.makeSuite(test_anova, 'check_') )
    
    
    total_suite = unittest.TestSuite(suites)
    return total_suite

def test():
    all_tests = test_suite()
    runner = unittest.TextTestRunner()
    runner.run(all_tests)
    return runner


if __name__ == "__main__":
    test()
