# Author:  Travis Oliphant, 2002
#
# Further enhancements and tests added by numerous SciPy developers.
#

import warnings

from numpy.testing import TestCase, run_module_suite, assert_array_equal, \
    assert_almost_equal, assert_array_less, assert_array_almost_equal, \
    assert_raises, assert_

import scipy.stats as stats

import numpy as np
from numpy.random import RandomState


g1 = [1.006, 0.996, 0.998, 1.000, 0.992, 0.993, 1.002, 0.999, 0.994, 1.000]
g2 = [0.998, 1.006, 1.000, 1.002, 0.997, 0.998, 0.996, 1.000, 1.006, 0.988]
g3 = [0.991, 0.987, 0.997, 0.999, 0.995, 0.994, 1.000, 0.999, 0.996, 0.996]
g4 = [1.005, 1.002, 0.994, 1.000, 0.995, 0.994, 0.998, 0.996, 1.002, 0.996]
g5 = [0.998, 0.998, 0.982, 0.990, 1.002, 0.984, 0.996, 0.993, 0.980, 0.996]
g6 = [1.009, 1.013, 1.009, 0.997, 0.988, 1.002, 0.995, 0.998, 0.981, 0.996]
g7 = [0.990, 1.004, 0.996, 1.001, 0.998, 1.000, 1.018, 1.010, 0.996, 1.002]
g8 = [0.998, 1.000, 1.006, 1.000, 1.002, 0.996, 0.998, 0.996, 1.002, 1.006]
g9 = [1.002, 0.998, 0.996, 0.995, 0.996, 1.004, 1.004, 0.998, 0.999, 0.991]
g10= [0.991, 0.995, 0.984, 0.994, 0.997, 0.997, 0.991, 0.998, 1.004, 0.997]


class TestShapiro(TestCase):
    def test_basic(self):
        x1 = [0.11,7.87,4.61,10.14,7.95,3.14,0.46,
              4.43,0.21,4.75,0.71,1.52,3.24,
              0.93,0.42,4.97,9.53,4.55,0.47,6.66]
        w,pw = stats.shapiro(x1)
        assert_almost_equal(w,0.90047299861907959,6)
        assert_almost_equal(pw,0.042089745402336121,6)
        x2 = [1.36,1.14,2.92,2.55,1.46,1.06,5.27,-1.11,
              3.48,1.10,0.88,-0.51,1.46,0.52,6.20,1.69,
              0.08,3.67,2.81,3.49]
        w,pw = stats.shapiro(x2)
        assert_almost_equal(w,0.9590270,6)
        assert_almost_equal(pw,0.52460,3)

    def test_bad_arg(self):
        # Length of x is less than 3.
        x = [1]
        assert_raises(ValueError, stats.shapiro, x)


class TestAnderson(TestCase):
    def test_normal(self):
        rs = RandomState(1234567890)
        x1 = rs.standard_exponential(size=50)
        x2 = rs.standard_normal(size=50)
        A,crit,sig = stats.anderson(x1)
        assert_array_less(crit[:-1], A)
        A,crit,sig = stats.anderson(x2)
        assert_array_less(A, crit[-2:])

    def test_expon(self):
        rs = RandomState(1234567890)
        x1 = rs.standard_exponential(size=50)
        x2 = rs.standard_normal(size=50)
        A,crit,sig = stats.anderson(x1,'expon')
        assert_array_less(A, crit[-2:])
        olderr = np.seterr(all='ignore')
        try:
            A,crit,sig = stats.anderson(x2,'expon')
        finally:
            np.seterr(**olderr)
        assert_(A > crit[-1])

    def test_bad_arg(self):
        assert_raises(ValueError, stats.anderson, [1], dist='plate_of_shrimp')


class TestAnsari(TestCase):

    def test_small(self):
        x = [1,2,3,3,4]
        y = [3,2,6,1,6,1,4,1]
        W, pval = stats.ansari(x,y)
        assert_almost_equal(W,23.5,11)
        assert_almost_equal(pval,0.13499256881897437,11)

    def test_approx(self):
        ramsay = np.array((111, 107, 100, 99, 102, 106, 109, 108, 104, 99,
                           101, 96, 97, 102, 107, 113, 116, 113, 110, 98))
        parekh = np.array((107, 108, 106, 98, 105, 103, 110, 105, 104,
                           100, 96, 108, 103, 104, 114, 114, 113, 108, 106, 99))
        W, pval = stats.ansari(ramsay, parekh)
        assert_almost_equal(W,185.5,11)
        assert_almost_equal(pval,0.18145819972867083,11)

    def test_exact(self):
        W,pval = stats.ansari([1,2,3,4],[15,5,20,8,10,12])
        assert_almost_equal(W,10.0,11)
        assert_almost_equal(pval,0.533333333333333333,7)

    def test_bad_arg(self):
        assert_raises(ValueError, stats.ansari, [], [1])
        assert_raises(ValueError, stats.ansari, [1], [])

warnings.filterwarnings('ignore',
                        message="Ties preclude use of exact statistic.")


class TestBartlett(TestCase):

    def test_data(self):
        args = [g1, g2, g3, g4, g5, g6, g7, g8, g9, g10]
        T, pval = stats.bartlett(*args)
        assert_almost_equal(T,20.78587342806484,7)
        assert_almost_equal(pval,0.0136358632781,7)

    def test_bad_arg(self):
        """Too few args raises ValueError."""
        assert_raises(ValueError, stats.bartlett, [1])


class TestLevene(TestCase):

    def test_data(self):
        args = [g1, g2, g3, g4, g5, g6, g7, g8, g9, g10]
        W, pval = stats.levene(*args)
        assert_almost_equal(W,1.7059176930008939,7)
        assert_almost_equal(pval,0.0990829755522,7)

    def test_trimmed1(self):
        """Test that center='trimmed' gives the same result as center='mean' when proportiontocut=0."""
        W1, pval1 = stats.levene(g1, g2, g3, center='mean')
        W2, pval2 = stats.levene(g1, g2, g3, center='trimmed', proportiontocut=0.0)
        assert_almost_equal(W1, W2)
        assert_almost_equal(pval1, pval2)

    def test_trimmed2(self):
        x = [1.2, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 100.0]
        y = [0.0, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 200.0]
        # Use center='trimmed'
        W1, pval1 = stats.levene(x, y, center='trimmed', proportiontocut=0.125)
        # Trim the data here, and use center='mean'
        W2, pval2 = stats.levene(x[1:-1], y[1:-1], center='mean')
        # Result should be the same.
        assert_almost_equal(W1, W2)
        assert_almost_equal(pval1, pval2)

    def test_equal_mean_median(self):
        x = np.linspace(-1,1,21)
        y = x**3
        W1, pval1 = stats.levene(x, y, center='mean')
        W2, pval2 = stats.levene(x, y, center='median')
        assert_almost_equal(W1, W2)
        assert_almost_equal(pval1, pval2)

    def test_bad_keyword(self):
        x = np.linspace(-1,1,21)
        assert_raises(TypeError, stats.levene, x, x, portiontocut=0.1)

    def test_bad_center_value(self):
        x = np.linspace(-1,1,21)
        assert_raises(ValueError, stats.levene, x, x, center='trim')

    def test_too_few_args(self):
        assert_raises(ValueError, stats.levene, [1])


class TestBinomP(TestCase):

    def test_data(self):
        pval = stats.binom_test(100,250)
        assert_almost_equal(pval,0.0018833009350757682,11)
        pval = stats.binom_test(201,405)
        assert_almost_equal(pval,0.92085205962670713,11)
        pval = stats.binom_test([682,243],p=3.0/4)
        assert_almost_equal(pval,0.38249155957481695,11)

    def test_bad_len_x(self):
        """Length of x must be 1 or 2."""
        assert_raises(ValueError, stats.binom_test, [1,2,3])

    def test_bad_n(self):
        """len(x) is 1, but n is invalid."""
        # Missing n
        assert_raises(ValueError, stats.binom_test, [100])
        # n less than x[0]
        assert_raises(ValueError, stats.binom_test, [100], n=50)

    def test_bad_p(self):
        assert_raises(ValueError, stats.binom_test, [50, 50], p=2.0)


class TestFindRepeats(TestCase):
    def test_basic(self):
        a = [1,2,3,4,1,2,3,4,1,2,5]
        res,nums = stats.find_repeats(a)
        assert_array_equal(res,[1,2,3,4])
        assert_array_equal(nums,[3,3,2,2])


class TestFligner(TestCase):

    def test_data(self):
        # numbers from R: fligner.test in package stats
        x1 = np.arange(5)
        assert_array_almost_equal(stats.fligner(x1,x1**2),
                           (3.2282229927203536, 0.072379187848207877), 11)

    def test_trimmed1(self):
        """Test that center='trimmed' gives the same result as center='mean' when proportiontocut=0."""
        Xsq1, pval1 = stats.fligner(g1, g2, g3, center='mean')
        Xsq2, pval2 = stats.fligner(g1, g2, g3, center='trimmed', proportiontocut=0.0)
        assert_almost_equal(Xsq1, Xsq2)
        assert_almost_equal(pval1, pval2)

    def test_trimmed2(self):
        x = [1.2, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 100.0]
        y = [0.0, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 200.0]
        # Use center='trimmed'
        Xsq1, pval1 = stats.fligner(x, y, center='trimmed', proportiontocut=0.125)
        # Trim the data here, and use center='mean'
        Xsq2, pval2 = stats.fligner(x[1:-1], y[1:-1], center='mean')
        # Result should be the same.
        assert_almost_equal(Xsq1, Xsq2)
        assert_almost_equal(pval1, pval2)

    # The following test looks reasonable at first, but fligner() uses the
    # function stats.rankdata(), and in one of the cases in this test,
    # there are ties, while in the other (because of normal rounding
    # errors) there are not.  This difference leads to differences in the
    # third significant digit of W.
    #
    #def test_equal_mean_median(self):
    #    x = np.linspace(-1,1,21)
    #    y = x**3
    #    W1, pval1 = stats.fligner(x, y, center='mean')
    #    W2, pval2 = stats.fligner(x, y, center='median')
    #    assert_almost_equal(W1, W2)
    #    assert_almost_equal(pval1, pval2)

    def test_bad_keyword(self):
        x = np.linspace(-1,1,21)
        assert_raises(TypeError, stats.fligner, x, x, portiontocut=0.1)

    def test_bad_center_value(self):
        x = np.linspace(-1,1,21)
        assert_raises(ValueError, stats.fligner, x, x, center='trim')

    def test_bad_num_args(self):
        """Too few args raises ValueError."""
        assert_raises(ValueError, stats.fligner, [1])


def test_mood():
    # numbers from R: mood.test in package stats
    x1 = np.arange(5)
    assert_array_almost_equal(stats.mood(x1,x1**2),
            (-1.3830857299399906, 0.16663858066771478), 11)

def test_mood_bad_arg():
    """Raise ValueError when the sum of the lengths of the args is less than 3."""
    assert_raises(ValueError, stats.mood, [1], [])

def test_oneway_bad_arg():
    """Raise ValueError is fewer than two args are given."""
    assert_raises(ValueError, stats.oneway, [1])

def test_wilcoxon_bad_arg():
    """Raise ValueError when two args of different lengths are given."""
    assert_raises(ValueError, stats.wilcoxon, [1], [1,2])

def test_mvsdist_bad_arg():
    """Raise ValueError if fewer than two data points are given."""
    data = [1]
    assert_raises(ValueError, stats.mvsdist, data)

def test_kstat_bad_arg():
    """Raise ValueError if n > 4 or n > 1."""
    data = [1]
    n = 10
    assert_raises(ValueError, stats.kstat, data, n=n)

def test_kstatvar_bad_arg():
    """Raise ValueError is n is not 1 or 2."""
    data = [1]
    n = 10
    assert_raises(ValueError, stats.kstatvar, data, n=n)

def test_probplot_bad_arg():
    """Raise ValueError when given an invalid distribution."""
    data = [1]
    assert_raises(ValueError, stats.probplot, data, dist="plate_of_shrimp")

def test_ppcc_max_bad_arg():
    """Raise ValueError when given an invalid distribution."""
    data = [1]
    assert_raises(ValueError, stats.ppcc_max, data, dist="plate_of_shrimp")

def test_boxcox_bad_arg():
    """Raise ValueError if any data value is negative."""
    x = np.array([-1])
    assert_raises(ValueError, stats.boxcox, x)

if __name__ == "__main__":
    run_module_suite()
