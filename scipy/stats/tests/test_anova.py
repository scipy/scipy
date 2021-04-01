
import pytest
import numpy as np
from numpy.testing import assert_equal, assert_allclose
from scipy.stats import anova
from scipy.special import fdtrc
from .._anova_util import _nway_groups


class TestOneway:

    @pytest.mark.parametrize('use_labeled', [False, True])
    def test_zar_10_1(self, use_labeled):
        # This is Example 10.1 from Zar, Biostatistical Analysis (5th ed).
        # When use_labeled is True, the data is put into two 1D lists,
        # and the result is computed with anova.oneway_from_labels_values.

        f1 = [60.8, 67.0, 65.0, 68.6, 61.7]
        f2 = [68.7, 67.7, 75.0, 73.3, 71.8]
        f3 = [69.6, 77.1, 75.2, 71.5]
        f4 = [61.9, 64.2, 63.1, 66.7, 60.3]

        if use_labeled:
            labels = sum([[k]*len(f) for k, f in enumerate([f1, f2, f3, f4])],
                         [])
            values = f1 + f2 + f3 + f4
            result = anova.oneway_from_labels_values(labels, values)
        else:
            result = anova.oneway(f1, f2, f3, f4)

        # The expected results can be found in the conclusion of
        # Example 10.1 on page 198.
        assert_equal(result.DFb, 3)
        assert_equal(result.DFw, 15)
        assert_allclose(result.SSb, 338.9374, rtol=5e-7)
        assert_allclose(result.SSw, 140.7500, rtol=5e-7)
        assert_allclose(result.MSb, 112.9791, rtol=1e-6)
        assert_allclose(result.MSw, 9.3833, rtol=1e-5)
        assert_allclose(result.pvalue, 0.00029, rtol=0.025)
        assert_allclose(result.F, 12.04, rtol=1e-4)

    def test_zar_10_2(self):
        # This is Example 10.2 from Zar, Biostatistical Analysis (5th ed).
        # The completed worked example is on page 201.

        t1 = [34, 36, 34, 35, 34]
        t2 = [37, 36, 35, 37, 37]
        t3 = [34, 37, 35, 37, 36]
        t4 = [36, 34, 37, 34, 35]

        result = anova.oneway(t1, t2, t3, t4)

        assert_equal(result.DFb, 3)
        assert_equal(result.DFw, 16)
        assert_allclose(result.MSb, 3.00, rtol=1e-12)
        assert_allclose(result.MSw, 1.25, rtol=1e-12)
        assert_allclose(result.F, 2.40, rtol=1e-12)
        assert_allclose(result.pvalue, 0.11, rtol=0.05)

    def test_sokalrohlf_9_1(self):
        # This is the example shown in Table 9.1 of Sokal and Rohlf,
        # Biometry (4th ed.), page 219.

        ctrl = np.array([1.33, 1.49, 1.43, 1.33, 1.54, 1.41, 1.49, 1.49,
                         1.32, 1.47])
        glu2 = np.array([1.75, 1.72, 1.67, 1.69, 1.61, 1.67, 1.67, 1.75,
                         1.69, 1.64])
        fru2 = np.array([1.72, 1.64, 1.79, 1.72, 1.75, 1.79, 1.64, 1.67,
                         1.75, 1.72])
        gf11 = np.array([1.72, 1.69, 1.72, 1.64, 1.75, 1.79, 1.72, 1.75,
                         1.75, 1.69])
        suc2 = np.array([1.61, 1.52, 1.54, 1.59, 1.56, 1.61, 1.54, 1.54,
                         1.61, 1.49])

        result = anova.oneway(ctrl, glu2, fru2, gf11, suc2)

        assert_equal(result.DFb, 4)
        assert_equal(result.DFw, 45)
        assert_allclose(result.MSb, 0.16019, rtol=2e-5)
        assert_allclose(result.MSw, 0.00297, rtol=1e-3)
        assert_allclose(result.F, 53.88, rtol=1e-4)
        assert_allclose(result.pvalue, 1e-16, rtol=0.5)

    def test_r_example(self):
        # An example of one-way ANOVA with R.
        #
        # R code:
        #
        #   a = factor(c(1, 1, 1, 2, 2, 2, 2, 3, 3, 3))
        #   v = c(1.0, 1.1, 0.4, 5.5, 1.2, 4.5, 2.2, 3.1, 1.2, 1.4)
        #   av = aov(v ~ a)
        #   options(digits=10)
        #   summary(av)
        #
        # R output:
        #
        #               Df   Sum Sq  Mean Sq F value  Pr(>F)
        #   a            2 11.14733 5.573667  2.7176 0.13383
        #   Residuals    7 14.35667 2.050952

        a = [1, 1, 1, 2, 2, 2, 2, 3, 3, 3]
        v = [1.0, 1.1, 0.4, 5.5, 1.2, 4.5, 2.2, 3.1, 1.2, 1.4]
        result = anova.oneway_from_labels_values(labels=a, values=v)

        assert_equal(result.DFb, 2)
        assert_equal(result.DFw, 7)
        assert_allclose(result.SSb, 11.14733, rtol=5e-7)
        assert_allclose(result.SSw, 14.35667, rtol=5e-7)
        assert_allclose(result.MSb, 5.573667, rtol=5e-7)
        assert_allclose(result.MSw, 2.050952, rtol=5e-7)
        assert_allclose(result.F, 2.7176, rtol=5e-5)
        assert_allclose(result.pvalue, 0.13383, rtol=5e-5)


class TestTwowayBalanced:

    def test_zar_12_1(self):
        # Example 12.1 (continued in Example 12.2) from
        # Zar, Biostatistical Analysis (fifth ed.), pages 251-255.

        conc = np.array(
            [[[16.3, 20.4, 12.4, 15.8, 9.5],
              [15.3, 17.4, 10.9, 10.3, 6.7]],
             [[38.1, 26.2, 32.3, 35.8, 30.2],
              [34.0, 22.8, 27.8, 25.0, 29.3]]])

        result = anova.twoway_from_data_grid(conc)

        # The expected results are summarized in Example 12.2
        # on page 255.
        assert_equal(result.DFA, 1)
        assert_equal(result.DFB, 1)
        assert_equal(result.DFAB, 1)
        assert_equal(result.DFerror, 16)

        assert_allclose(result.SSA, 1386.1125, rtol=5e-8)
        assert_allclose(result.SSB, 70.3125, rtol=5e-6)
        assert_allclose(result.SSAB, 4.9005, rtol=5e-5)
        assert_allclose(result.SSerror, 301.3920, rtol=5e-7)

        assert_allclose(result.MSA, 1386.1125, rtol=5e-8)
        assert_allclose(result.MSB, 70.3125, rtol=5e-6)
        assert_allclose(result.MSAB, 4.9005, rtol=5e-5)
        assert_allclose(result.MSerror, 18.8370, rtol=5e-6)

        assert_allclose(result.FA, 73.6, rtol=5e-3)
        assert_allclose(result.FB, 3.73, rtol=5e-3)
        assert_allclose(result.FAB, 0.260, rtol=5e-3)

        assert_allclose(result.pA, 0.00000022, rtol=5e-2)
        assert_allclose(result.pB, 0.071, rtol=5e-2)
        assert_allclose(result.pAB, 0.62, rtol=5e-2)

    def test_sokalrohlf_11_1(self):
        # Example from Sokal & Rohlf, Biometry (fourth ed.),
        # Box 11.1, page 322.

        consumption = np.array(
            [[[709, 679, 699],
              [592, 538, 476]],
             [[657, 594, 677],
              [508, 505, 539]]])

        result = anova.twoway_from_data_grid(consumption)

        assert_equal(result.DFA, 1)
        assert_equal(result.DFB, 1)
        assert_equal(result.DFAB, 1)
        assert_equal(result.DFerror, 8)

        assert_allclose(result.SSA, 3780.75, rtol=5e-6)
        assert_allclose(result.SSB, 61204.08, rtol=5e-7)
        assert_allclose(result.SSAB, 918.75, rtol=5e-5)
        assert_allclose(result.SSerror, 11666.67, rtol=5e-7)

        assert_allclose(result.MSA, 3780.75, rtol=5e-6)
        assert_allclose(result.MSB, 61204.08, rtol=5e-7)
        assert_allclose(result.MSAB, 918.75, rtol=5e-5)
        assert_allclose(result.MSerror, 1458.33, rtol=5e-6)

        # The F statistics are not shown in Box 11.1.

        assert_allclose(result.pA, 0.1460, rtol=5e-4)
        assert_allclose(result.pB, 0.0002, rtol=0.5)
        assert_allclose(result.pAB, 0.4503, rtol=5e-4)

    def test_r_balanced(self):
        # A simple example of balanced two-way ANOVA in R, with
        # two replicates:
        #
        #   x = factor(c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3))
        #   y = factor(c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2))
        #   v = c(0.5, 0.5, 1.0, 2.0, 2.1, 1.7, 0.5, 0.8,
        #         0.2, 0.1, 0.2, 0.1)
        #   av = aov(v ~ x + y + x*y)
        #   options(digits=10)
        #   summary(av)
        #
        # R output:
        #
        #               Df    Sum Sq   Mean Sq F value  Pr(>F)
        #   x            2 2.7516667 1.3758333 2.81261 0.13748
        #   y            1 0.0408333 0.0408333 0.08348 0.78236
        #   x:y          2 0.2216667 0.1108333 0.22658 0.80378
        #   Residuals    6 2.9350000 0.4891667

        x = np.array([1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3])
        y = np.array([1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2])
        v = np.array([0.5, 0.5, 1.0, 2.0, 2.1, 1.7, 0.5, 0.8,
                      0.2, 0.1, 0.2, 0.1])
        result = anova.twoway_from_a_b_values(x, y, v, sstype=1)

        assert_equal(result.DFA, 2)
        assert_equal(result.DFB, 1)
        assert_equal(result.DFAB, 2)
        assert_equal(result.DFerror, 6)

        assert_allclose(result.SSA, 2.7516667, rtol=5e-8)
        assert_allclose(result.SSB, 0.0408333, rtol=5e-6)
        assert_allclose(result.SSAB, 0.2216667, rtol=5e-7)
        assert_allclose(result.SSerror, 2.9350000, rtol=5e-8)

        assert_allclose(result.MSA, 1.3758333, rtol=5e-8)
        assert_allclose(result.MSB, 0.0408333, rtol=5e-6)
        assert_allclose(result.MSAB, 0.1108333, rtol=5e-7)
        assert_allclose(result.MSerror, 0.4891667, rtol=5e-7)

        assert_allclose(result.FA, 2.81261, rtol=5e-6)
        assert_allclose(result.FB, 0.08348, rtol=5e-4)
        assert_allclose(result.FAB, 0.22658, rtol=5e-5)

        assert_allclose(result.pA, 0.13748, rtol=5e-5)
        assert_allclose(result.pB, 0.78236, rtol=5e-5)
        assert_allclose(result.pAB, 0.80378, rtol=5e-5)

    def test_draper_smith(self):
        # This example is from section 23.6 of
        #
        #   Draper and Smith, Applied Regression Analysis (third edition),
        #   Wiley, New York, USA (1998).
        #
        # The data is given in Table 23.4, and the results of the two-way
        # ANOVA (not including the p-values) are shown in Table 23.5.

        data = [[[4, 6], [11, 7], [5, 9]],
                [[6, 4], [13, 15], [9, 7]],
                [[13, 15], [15, 9], [13, 13]],
                [[12, 12], [12, 14], [7, 9]]]

        result = anova.twoway_from_data_grid(data)

        assert_equal(result.DFA, 3)
        assert_equal(result.DFB, 2)
        assert_equal(result.DFAB, 6)
        assert_equal(result.DFerror, 12)

        # With just 2 replicates per cell with integer values,
        # the calculation should be very accurate.
        assert_allclose(result.SSA, 120, rtol=1e-8)
        assert_allclose(result.SSB, 48, rtol=1e-8)
        assert_allclose(result.SSAB, 84, rtol=1e-8)
        assert_allclose(result.SSerror, 48, rtol=1e-8)

        assert_allclose(result.MSA, 40, rtol=1e-8)
        assert_allclose(result.MSB, 24, rtol=1e-8)
        assert_allclose(result.MSAB, 14, rtol=1e-8)
        assert_allclose(result.MSerror, 4, rtol=1e-8)

        assert_allclose(result.FA, 10.0, rtol=1e-8)
        assert_allclose(result.FB, 6.0, rtol=1e-8)
        assert_allclose(result.FAB, 3.5, rtol=1e-8)


class TestTwowayUnbalanced:

    def test_rcompanion_example(self):
        # From https://rcompanion.org/rcompanion/d_08.html
        # The data is unbalanced.
        acttext = """
          1 male    ff        1.884
          2 male    ff        2.283
          3 male    fs        2.396
          4 female  ff        2.838
          5 male    fs        2.956
          6 female  ff        4.216
          7 female  ss        3.620
          8 female  ff        2.889
          9 female  fs        3.550
         10 male    fs        3.105
         11 female  fs        4.556
         12 female  fs        3.087
         13 male    ff        4.939
         14 male    ff        3.486
         15 female  ss        3.079
         16 male    fs        2.649
         17 female  fs        1.943
         19 female  ff        4.198
         20 female  ff        2.473
         22 female  ff        2.033
         24 female  fs        2.200
         25 female  fs        2.157
         26 male    ss        2.801
         28 male    ss        3.421
         29 female  ff        1.811
         30 female  fs        4.281
         32 female  fs        4.772
         34 female  ss        3.586
         36 female  ff        3.944
         38 female  ss        2.669
         39 female  ss        3.050
         41 male    ss        4.275
         43 female  ss        2.963
         46 female  ss        3.236
         48 female  ss        3.673
         49 male    ss        3.110
        """

        actfields = [w.split()[1:] for w in acttext.strip().splitlines()]
        sex, genotype, activity = list(zip(*actfields))
        activity = [float(s) for s in activity]

        result = anova.twoway_from_a_b_values(a=sex, b=genotype,
                                              values=activity, sstype=1)

        # The expected results are shown on the web page in the R calculations
        # in the subsection titled "Fit the linear model and conduct ANOVA".
        assert_equal(result.DFA, 1)
        assert_equal(result.DFB, 2)
        assert_equal(result.DFAB, 2)
        assert_equal(result.DFerror, 30)

        assert_allclose(result.SSA, 0.0681, rtol=5e-3)
        assert_allclose(result.SSB, 0.2772, rtol=5e-4)
        assert_allclose(result.SSAB, 0.8146, rtol=5e-4)
        assert_allclose(result.SSerror, 23.7138, rtol=5e-6)

        assert_allclose(result.MSA, 0.06808, rtol=5e-4)
        assert_allclose(result.MSB, 0.13862, rtol=5e-5)
        assert_allclose(result.MSAB, 0.40732, rtol=5e-5)
        assert_allclose(result.MSerror, 0.79046, rtol=5e-5)

        assert_allclose(result.FA, 0.0861, rtol=5e-3)
        assert_allclose(result.FB, 0.1754, rtol=5e-4)
        assert_allclose(result.FAB, 0.5153, rtol=5e-4)

        assert_allclose(result.pA, 0.7712, rtol=5e-4)
        assert_allclose(result.pB, 0.8400, rtol=5e-4)
        assert_allclose(result.pAB, 0.6025, rtol=5e-4)

    @pytest.mark.parametrize('sstype', [1, 2, 3])
    def test_r_unbalanced(self, sstype):
        # A two-way unbalanced example.  Test with ANOVA types 1, 2, and 3.
        #
        # R code:
        #
        #   options(digits=10)
        #   a = as.factor(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,
        #                   1, 2, 2, 2, 2, 2, 2, 2))
        #   b = as.factor(c(0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 0, 0, 0, 0, 1, 2, 2,
        #                   3, 0, 0, 1, 1, 2, 3, 3))
        #   vals = c(1, 1, 2, 1, 1, 2, 4, 4, 4, 3, 9, 3, 1, 4, 4, 4, 4,
        #            7, 6, 5, 2, 3, 4, 5, 6)
        #   summary(aov(vals ~ a + b + a*b))
        #
        # Output:
        #
        #              Df   Sum Sq   Mean Sq F value   Pr(>F)
        #  a            2 28.18571 14.092857 4.81069 0.027306 *
        #  b            3 19.88442  6.628140 2.26256 0.129505
        #  a:b          6 13.84653  2.307755 0.78777 0.594972
        #  Residuals   13 38.08333  2.929487
        #
        # R code for Type II and Type III ANOVA (continued from above code)
        #
        #   library(car)
        #   Anova(lm(vals ~ a + b + a*b), type='II')
        #   Anova(lm(vals ~ a * b, contrasts=list(a=contr.sum, b=contr.sum)),
        #         type='III"')
        #
        # Output:
        #
        #   Anova Table (Type II tests)
        #
        #   Response: vals
        #                Sum Sq Df F value   Pr(>F)
        #   a         25.125691  2 4.28841 0.037129 *
        #   b         19.884421  3 2.26256 0.129505
        #   a:b       13.846531  6 0.78777 0.594972
        #   Residuals 38.083333 13
        #
        #   Anova Table (Type III tests)
        #
        #   Response: vals
        #                  Sum Sq Df   F value     Pr(>F)
        #   (Intercept) 318.24197  1 108.63402 1.1105e-07 ***
        #   a            22.99267  2   3.92435   0.046411 *
        #   b            17.73820  3   2.01835   0.161123
        #   a:b          13.84653  6   0.78777   0.594972
        #   Residuals    38.08333 13
        #
        a = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,
                      1, 2, 2, 2, 2, 2, 2, 2])
        b = np.array([0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 0, 0, 0, 0, 1, 2, 2,
                      3, 0, 0, 1, 1, 2, 3, 3])
        vals = np.array([1, 1, 2, 1, 1, 2, 4, 4, 4, 3, 9, 3, 1, 4, 4, 4, 4,
                         7, 6, 5, 2, 3, 4, 5, 6])

        result = anova.twoway_from_a_b_values(a, b, vals, sstype=sstype)

        assert_equal(result.DFA, 2)
        assert_equal(result.DFB, 3)
        assert_equal(result.DFAB, 6)
        assert_equal(result.DFerror, 13)
        assert_allclose(result.SSAB, 13.846531, rtol=5e-8)
        assert_allclose(result.SSerror, 38.083333, rtol=5e-8)
        assert_allclose(result.MSAB, 2.307755, rtol=5e-7)
        assert_allclose(result.MSerror, 2.929487, rtol=5e-7)
        assert_allclose(result.FAB, 0.78777, rtol=5e-5)
        assert_allclose(result.pAB, 0.594972, rtol=5e-5)

        if sstype == 1:
            assert_allclose(result.SSA, 28.18571, rtol=5e-7)
            assert_allclose(result.SSB, 19.88442, rtol=5e-7)
            assert_allclose(result.MSA, 14.092857, rtol=5e-8)
            assert_allclose(result.MSB, 6.628140, rtol=5e-7)
            assert_allclose(result.FA, 4.81069, rtol=5e-6)
            assert_allclose(result.FB, 2.26256, rtol=5e-6)
            assert_allclose(result.pA, 0.027306, rtol=5e-5)
            assert_allclose(result.pB, 0.129505, rtol=5e-5)
        elif sstype == 2:
            assert_allclose(result.SSA, 25.125691, rtol=5e-8)
            assert_allclose(result.SSB, 19.884421, rtol=5e-8)
            assert_allclose(result.MSA, result.SSA/result.DFA, rtol=5e-8)
            assert_allclose(result.MSB, result.SSB/result.DFB, rtol=5e-8)
            assert_allclose(result.FA, 4.28841, rtol=5e-6)
            assert_allclose(result.FB, 2.26256, rtol=5e-6)
            assert_allclose(result.pA, 0.037129, rtol=5e-5)
            assert_allclose(result.pB, 0.129505, rtol=5e-6)
        else:
            # sstype == 3:
            assert_allclose(result.SSA, 22.99267, rtol=5e-7)
            assert_allclose(result.SSB, 17.73820, rtol=5e-7)
            assert_allclose(result.MSA, result.SSA/result.DFA, rtol=5e-8)
            assert_allclose(result.MSB, result.SSB/result.DFB, rtol=5e-8)
            assert_allclose(result.FA, 3.92435, rtol=5e-6)
            assert_allclose(result.FB, 2.01835, rtol=5e-6)
            assert_allclose(result.pA, 0.046411, rtol=5e-5)
            assert_allclose(result.pB, 0.161123, rtol=5e-6)

    @pytest.mark.parametrize('sstype', [1, 2, 3])
    def test_smith_cribbie(self, sstype):
        # Apply ANOVA to the example used in
        #
        #   Smith and Cribbie, Factorial ANOVA with unbalanced data: A
        #   fresh look at the types of sums of squares, Journal of Data
        #   Science 12(2014), 385-404.
        #
        # The paper shows results for types 1, 2, and 3.

        data = [[[3.0, 3.0, 2.8],
                 [5.1, 5.2, 4.7, 4.9, 4.9, 5.0],
                 [2.1, 1.9, 2.0, 1.8]],
                [[2.3, 2.4, 2.1],
                 [3.9, 4.1, 3.8],
                 [1.2, 1.3, 1.1, 1.1, 1.0]]]

        result = anova.twoway_from_data_grid(data, sstype=sstype)

        assert_equal(result.DFA, 1)
        assert_equal(result.DFB, 2)
        assert_equal(result.DFAB, 2)
        assert_equal(result.DFerror, 18)
        assert_allclose(result.SSAB, 0.121, rtol=5e-3)
        assert_allclose(result.SSerror, 0.375, rtol=5e-3)
        assert_allclose(result.MSAB, 0.061, rtol=5e-2)
        assert_allclose(result.MSerror, 0.021, rtol=5e-2)
        assert_allclose(result.FAB, 2.906, rtol=5e-4)
        assert_allclose(result.pAB, 0.081, rtol=5e-2)

        # The paper reports the very small probabilities pA and pB
        # as "< 0.0001".  Here we'll recompute the p-value using fdtrc.
        assert_allclose(result.pA,
                        fdtrc(result.DFA, result.DFerror, result.FA),
                        rtol=1e-8)
        assert_allclose(result.pB,
                        fdtrc(result.DFB, result.DFerror, result.FB),
                        rtol=1e-8)

        if sstype == 1:
            assert_allclose(result.SSA, 11.023, rtol=5e-5)
            assert_allclose(result.SSB, 37.94, rtol=5e-4)
            assert_allclose(result.MSA, 11.023, rtol=5e-5)
            assert_allclose(result.MSB, 18.97, rtol=5e-4)
            assert_allclose(result.FA, 528.625, rtol=5e-6)
            assert_allclose(result.FB, 909.757, rtol=5e-6)
        elif sstype == 2:
            assert_allclose(result.SSA, 4.139, rtol=5e-4)
            assert_allclose(result.SSB, 37.94, rtol=5e-4)
            assert_allclose(result.MSA, 4.139, rtol=5e-4)
            assert_allclose(result.MSB, 18.97, rtol=5e-4)
            assert_allclose(result.FA, 198.497, rtol=5e-6)
            assert_allclose(result.FB, 909.757, rtol=5e-6)
        else:
            # sstype == 3:
            assert_allclose(result.SSA, 3.897, rtol=5e-4)
            assert_allclose(result.SSB, 35.989, rtol=5e-5)
            assert_allclose(result.MSA, 3.897, rtol=5e-4)
            assert_allclose(result.MSB, 17.995, rtol=5e-5)
            assert_allclose(result.FA, 186.888, rtol=5e-6)
            assert_allclose(result.FB, 862.971, rtol=5e-6)

    @pytest.mark.parametrize('sstype', [1, 2, 3])
    def test_sas_example(self, sstype):
        # Compare results to the SAS example at
        # https://support.sas.com/documentation/cdl/en/statug/63347/
        #         HTML/default/viewer.htm#statug_glm_sect049.htm
        #
        text = """
            1 1 42 44 36 13 19 22
            1 2 33  . 26  . 33 21
            1 3 31 -3  . 25 25 24
            2 1 28  . 23 34 42 13
            2 2  . 34 33 31  . 36
            2 3  3 26 28 32  4 16
            3 1  .  .  1 29  . 19
            3 2  . 11  9  7  1 -6
            3 3 21  1  .  9  3  .
            4 1 24  .  9 22 -2 15
            4 2 27 12 12 -5 16 15
            4 3 22  7 25  5 12  .
        """
        data = np.empty((4, 3), dtype=object)
        for row in text.strip().splitlines():
            fields = [int(field) for field in row.split() if field != '.']
            data[fields[0] - 1, fields[1] - 1] = fields[2:]

        result = anova.twoway_from_data_grid(data, sstype=sstype)

        assert_equal(result.DFA, 3)
        assert_equal(result.DFB, 2)
        assert_equal(result.DFAB, 6)
        assert_equal(result.DFerror, 46)

        assert_allclose(result.SSAB, 707.266259, rtol=5e-9)
        assert_allclose(result.SSerror, 5080.816667, rtol=5e-10)

        assert_allclose(result.MSAB, 117.877710, rtol=5e-9)
        assert_allclose(result.MSerror, 110.452536, rtol=5e-9)

        assert_allclose(result.FAB, 1.07, rtol=5e-3)
        assert_allclose(result.pAB, 0.3958, rtol=5e-4)

        if sstype == 1:
            assert_allclose(result.SSA, 3133.238506, rtol=5e-10)
            assert_allclose(result.SSB, 418.833741, rtol=5e-9)
            assert_allclose(result.MSA, 1044.412835, rtol=5e-10)
            assert_allclose(result.MSB, 209.416870, rtol=5e-9)
            assert_allclose(result.FA, 9.46, rtol=5e-3)
            assert_allclose(result.FB, 1.90, rtol=5e-3)
            assert result.pA < 0.0001
            assert_allclose(result.pB, 0.1617, rtol=5e-4)
        elif sstype == 2:
            assert_allclose(result.SSA, 3063.432863, rtol=5e-10)
            assert_allclose(result.SSB, 418.833741, rtol=5e-9)
            assert_allclose(result.MSA, 1021.144288, rtol=5e-10)
            assert_allclose(result.MSB, 209.416870, rtol=5e-9)
            assert_allclose(result.FA, 9.25, rtol=5e-3)
            assert_allclose(result.FB, 1.90, rtol=5e-3)
            assert result.pA < 0.0001
            assert_allclose(result.pB, 0.1617, rtol=5e-4)
        else:
            # sstype == 3:
            assert_allclose(result.SSA, 2997.471860, rtol=5e-10)
            assert_allclose(result.SSB, 415.873046, rtol=5e-9)
            assert_allclose(result.MSA, 999.157287, rtol=5e-9)
            assert_allclose(result.MSB, 207.936523, rtol=5e-9)
            assert_allclose(result.FA, 9.05, rtol=5e-3)
            assert_allclose(result.FB, 1.88, rtol=5e-3)
            assert result.pA < 0.0001
            assert_allclose(result.pB, 0.1637, rtol=5e-4)

    def test_neter_wasserman_kutner_22(self):
        # This example is from section 22.2 of
        #
        #   Neter, Wasserman and Kutner, Applied Linear Statistical
        #   Models (second edition), Irwin, Homewood, Illinois, USA (1985)
        #
        # The data is given in Table 22.1, and the results (without p-values)
        # are shown in Table 22.4.  The analysis in the text amounts to a
        # Type 3 ANOVA.

        data = [[[1.4, 2.4, 2.2], [2.1, 1.7], [0.7, 1.1]],
                [[2.4], [2.5, 1.8, 2.0], [0.5, 0.9, 1.3]]]

        result = anova.twoway_from_data_grid(data, sstype=3)

        assert_equal(result.DFA, 1)
        assert_equal(result.DFB, 2)
        assert_equal(result.DFAB, 2)
        assert_equal(result.DFerror, 8)

        assert_allclose(result.SSA, 0.1200, rtol=5e-4)
        assert_allclose(result.SSB, 4.1897, rtol=5e-5)
        assert_allclose(result.SSAB, 0.0754, rtol=5e-3)
        assert_allclose(result.SSerror, 1.3000, rtol=5e-5)

        assert_allclose(result.MSA, 0.1200, rtol=5e-4)
        assert_allclose(result.MSB, 2.0949, rtol=5e-5)
        assert_allclose(result.MSAB, 0.0377, rtol=5e-3)
        assert_allclose(result.MSerror, 0.1625, rtol=5e-4)

        assert_allclose(result.FA, 0.74, rtol=5e-2)
        assert_allclose(result.FB, 12.89, rtol=5e-4)
        assert_allclose(result.FAB, 0.23, rtol=5e-2)


class TestTwoway1:

    def test_biostathandbook_example(self):
        # From http://www.biostathandbook.com/twowayanova.html
        # Snake-ID    Day-1   Day-2   Day-3   Day-4
        # D1  85  58  15  57
        # D3  107 51  30  12
        # D5  61  60  68  36
        # D8  22  41  63  21
        # D11 40  45  28  10
        # D12 65  27  3   16

        data = np.array([[85, 58, 15, 57],
                         [107, 51, 30, 12],
                         [61, 60, 68, 36],
                         [22, 41, 63, 21],
                         [40, 45, 28, 10],
                         [65, 27, 3, 16]])

        # The web page reports:
        #   The effect of snake is not significant (F[5,15]=1.24, P=0.34),
        #   while the effect of day is significant (F[3,15]=3.32, P=0.049).
        # That agrees with the output of the following:
        result = anova.twoway_from_data_grid(data)

        # Test just the values explicitly shown in the web page.
        assert_equal(result.DFA, 5)
        assert_equal(result.DFB, 3)
        assert_equal(result.DFerror, 15)

        assert_allclose(result.FA, 1.24, rtol=5e-3)
        assert_allclose(result.FB, 3.32, rtol=5e-3)

        assert_allclose(result.pA, 0.34, rtol=0.05)
        assert_allclose(result.pB, 0.049, rtol=0.05)

    def test_sokalrohlf_11_4(self):
        # Example from Sokal & Rohlf, Biometry (fourth ed.), Box 11.4
        # (page 340 - 342).
        data = np.array([[23.8, 24.0, 24.6, 24.8],
                         [22.6, 22.4, 22.9, 23.2],
                         [22.2, 22.1, 22.1, 22.2],
                         [21.2, 21.8, 21.0, 21.2],
                         [18.4, 19.3, 19.0, 18.8],
                         [13.5, 14.4, 14.2, 13.8],
                         [9.80, 9.90, 10.4, 9.60],
                         [6.00, 6.00, 6.30, 6.30],
                         [5.80, 5.90, 6.00, 5.80],
                         [5.60, 5.60, 5.50, 5.60]])
        # Note that we pass in data.T to swap the roles of the
        # rows and columns.  This is only to match the convention
        # used in the text, where "A" refers to the "days" factor,
        # and "B" refers to the "depth" factor.
        result = anova.twoway_from_data_grid(data.T)

        assert_equal(result.DFA, 3)
        assert_equal(result.DFB, 9)
        assert_equal(result.DFerror, 27)

        assert_allclose(result.SSA, 0.5620, rtol=5e-4)
        assert_allclose(result.SSB, 2119.6510, rtol=5e-8)
        assert_allclose(result.SSerror, 2.2430, rtol=5e-5)

        assert_allclose(result.MSA, 0.1873, rtol=5e-4)
        assert_allclose(result.MSB, 235.5168, rtol=5e-7)
        assert_allclose(result.MSerror, 0.0831, rtol=5e-3)

        assert_allclose(result.FA, 2.255, rtol=5e-4)
        assert_allclose(result.FB, 2835.0, rtol=5e-5)

        assert_allclose(result.pA, 0.1048, rtol=5e-4)
        # The book reports pB < 1e-16. The computed value should be
        #
        #   >>> fdtrc(9, 27, 2835.0)
        #   8.829369277294019e-38
        #
        expected_pB = fdtrc(result.DFB, result.DFerror, result.FB)
        assert_allclose(result.pB, expected_pB, rtol=1e-12)

    def test_r_example_1_rep(self):
        # A simple example of two-way ANOVA with a single replicate;
        # an interaction term is not included in the formula.
        #
        # R code:
        #
        #   x = factor(c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3))
        #   y = factor(c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5))
        #   v = c(1.0, 2.0, 2.1, 1.7, 0.5, 0.8, 1.0, 4.0, 2.5, 2.5, 0.1,
        #         0.2, 0.1, 0.1, 0.2)
        #   av = aov(v ~ x + y)
        #   options(digits=10)
        #   summary(av)
        #
        # R output:
        #
        #               Df   Sum Sq  Mean Sq F value   Pr(>F)
        #   x            2 10.52133 5.260667 7.98280 0.012417 *
        #   y            4  3.44400 0.861000 1.30653 0.345417
        #   Residuals    8  5.27200 0.659000
        #   ---
        #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
        #

        a = np.array([1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3])
        b = np.array([1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5])
        v = np.array([1.0, 2.0, 2.1, 1.7, 0.5, 0.8, 1.0, 4.0, 2.5, 2.5, 0.1,
                      0.2, 0.1, 0.1, 0.2])
        result = anova.twoway_from_a_b_values(a, b, v, sstype=1)

        assert_equal(result.DFA, 2)
        assert_equal(result.DFB, 4)
        assert_equal(result.DFerror, 8)

        assert_allclose(result.SSA, 10.52133, rtol=5e-7)
        assert_allclose(result.SSB, 3.44400, rtol=5e-6)
        assert_allclose(result.SSerror, 5.27200, rtol=5e-6)

        assert_allclose(result.MSA, 5.260667, rtol=5e-7)
        assert_allclose(result.MSB, 0.861000, rtol=5e-6)
        assert_allclose(result.MSerror, 0.659000, rtol=5e-6)

        assert_allclose(result.FA, 7.98280, rtol=5e-3)
        assert_allclose(result.FB, 1.30653, rtol=5e-3)

        assert_allclose(result.pA, 0.012417, rtol=5e-5)
        assert_allclose(result.pB, 0.345417, rtol=5e-6)


class TestTwowayNestedBalanced:
    """
    Tests of two-way nested ANOVA in which the number of replicates
    is the same in all the cells.
    """

    def test_schools_teachers(self):
        # http://www.personal.psu.edu/mar36/stat_461/nested/nested_designs.html
        # Scroll down to the line that says "What if we do it correctly?"
        #
        # R code:
        #
        #   scores <- c(25, 29, 14, 11, 11, 6, 22, 18, 17, 20, 5, 2)
        #   school <- factor(c("A", "A", "A", "A", "B", "B", "B", "B",
        #                      "C", "C", "C", "C"))
        #   teacher <- factor(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6))
        #   options(digits=10)
        #   res1 <- lm(scores ~ school + school/teacher)
        #   anova(res1)
        #
        # Alternatively,
        #
        #   summary(aov(scores ~ school + school/teacher))
        #
        # R output:
        #
        #   Analysis of Variance Table
        #
        #   Response: scores
        #                  Df Sum Sq   Mean Sq  F value     Pr(>F)
        #   school          2  156.5  78.25000 11.17857 0.00947254 **
        #   school:teacher  3  567.5 189.16667 27.02381 0.00069701 ***
        #   Residuals       6   42.0   7.00000
        #

        scores = np.array([25, 29, 14, 11, 11, 6, 22, 18, 17, 20, 5, 2])
        data = scores.reshape(3, -1, 2)

        result = anova.twoway_nested_from_data_grid(data, btype='fixed')

        assert_equal(result.DFA, 2)
        assert_equal(result.DFBA, 3)
        assert_equal(result.DFerror, 6)

        assert_allclose(result.SSA, 156.5, rtol=5e-4)
        assert_allclose(result.SSBA, 567.5, rtol=5e-4)
        assert_allclose(result.SSerror, 42.0, rtol=5e-3)

        assert_allclose(result.MSA, 78.25000, rtol=5e-7)
        assert_allclose(result.MSBA, 189.16667, rtol=5e-8)
        assert_allclose(result.MSerror, 7.00000, rtol=5e-6)

        assert_allclose(result.FA, 11.17857, rtol=5e-7)
        assert_allclose(result.FBA, 27.02381, rtol=5e-7)

        assert_allclose(result.pA, 0.00947254, rtol=5e-6)
        assert_allclose(result.pBA, 0.00069701, rtol=5e-5)

    def test_sokalrohlf_box_10_1(self):
        # Example from Sokal & Rohlf, Biometry (fourth ed.), page 281-
        # Two-level nested ANOVA.
        # Factor A has 3 groups.
        # Factor B nested in A has 4 groups.
        # 2 measurements for each nested combination of factors.

        data = [[[58.5, 59.5],
                 [77.8, 80.9],
                 [84.0, 83.6],
                 [70.1, 68.3]],
                [[69.8, 69.8],
                 [56.0, 54.5],
                 [50.7, 49.3],
                 [63.8, 65.8]],
                [[56.6, 57.5],
                 [77.8, 79.2],
                 [69.9, 69.2],
                 [62.1, 64.5]]]

        result = anova.twoway_nested_from_data_grid(data)

        # The expected results can be found in the "Completed anova"
        # section of the table on page 282.
        assert_equal(result.DFA, 2)
        assert_equal(result.DFBA, 9)
        assert_equal(result.DFerror, 12)

        assert_allclose(result.SSA, 665.6758, rtol=5e-7)
        assert_allclose(result.SSBA, 1720.6775, rtol=5e-8)
        assert_allclose(result.SSerror, 15.6200, rtol=5e-6)

        assert_allclose(result.MSA, 332.8379, rtol=5e-7)
        assert_allclose(result.MSBA, 191.1864, rtol=5e-7)
        assert_allclose(result.MSerror, 1.3017, rtol=5e-5)

        assert_allclose(result.FA, 1.741, rtol=5e-4)
        assert_allclose(result.FBA, 146.88, rtol=5e-5)

        assert_allclose(result.pA, 0.2295, rtol=5e-4)
        # The p-value for the nested factor is reported to be 1.3467e-7
        # in the text, but this is not correct.  It appears they used
        # the same degrees of freedom as for the calculation of pA.  That
        # is, pA is computed as
        #
        #   >>> fdtrc(2, 9, 1.741)
        #   0.2295157484151867
        #
        # which gives the expected value 0.2295.  The value reported for pB
        # can be obtained with
        #
        #   >>> fdtrc(2, 9, 146.88)
        #   1.346319257591261e-07
        #
        # but those are the wrong degrees of free for FBA.  The calculation
        # should be
        #
        #   >>> fdtrc(9, 12, 146.88)
        #   6.980793718640759e-11
        #
        expected_pBA = fdtrc(result.DFBA, result.DFerror, result.FBA)
        assert_allclose(result.pBA, expected_pBA, rtol=1e-12)

    def test_zar_15_1(self):
        # This is Example 15.1 from Zar, Biostatistical Analysis (5th ed).

        g = np.array([[[102, 104],
                       [103, 104]],
                      [[108, 110],
                       [109, 108]],
                      [[104, 106],
                       [105, 107]]])

        result = anova.twoway_nested_from_data_grid(g)

        assert_equal(result.DFA, 2)
        assert_equal(result.DFBA, 3)
        assert_equal(result.DFerror, 6)

        assert_allclose(result.SSA, 61.17, rtol=5e-4)
        assert_allclose(result.SSBA, 1.50, rtol=5e-3)
        assert_allclose(result.SSerror, 9.00, rtol=5e-3)

        assert_allclose(result.MSA, 30.58, rtol=5e-4)
        assert_allclose(result.MSBA, 0.50, rtol=5e-3)
        assert_allclose(result.MSerror, 1.50, rtol=5e-3)

        assert_allclose(result.FA, 61.16, rtol=5e-4)
        assert_allclose(result.FBA, 0.33, rtol=5e-2)

        assert_allclose(result.pA, 0.0037, rtol=5e-2)
        assert_allclose(result.pBA, 0.80, rtol=5e-2)

    @pytest.mark.parametrize('btype', ['random', 'fixed'])
    def test_statsdirect_example(self, btype):
        # https://www.statsdirect.com/help/analysis_of_variance/nested.htm

        datastr = """
        3.28 3.52 2.88 2.46 1.87 2.19 2.77 3.74 2.55 3.78 4.07 3.31
        3.09 3.48 2.80 2.44 1.92 2.19 2.66 3.44 2.55 3.87 4.12 3.31
        """

        data = np.array([float(t) for t in datastr.split()]).reshape(2, 4, 3)
        data = np.moveaxis(data, 0, -1)
        # data.shape is now (4, 3, 2).

        result = anova.twoway_nested_from_data_grid(data, btype=btype)

        # The expected values are from the web page.
        # These parts are the same whether btype is 'fixed' or 'random'.
        assert_equal(result.DFA, 3)
        assert_equal(result.DFBA, 8)
        assert_equal(result.DFerror, 12)
        assert_allclose(result.SSA, 7.560346, rtol=5e-7)
        assert_allclose(result.SSBA, 2.6302, rtol=5e-5)
        assert_allclose(result.SSerror, 0.07985, rtol=5e-4)
        assert_allclose(result.MSA, 2.520115, rtol=5e-4)
        assert_allclose(result.MSBA, 0.328775, rtol=5e-3)
        assert_allclose(result.MSerror, 0.006654, rtol=5e-3)
        assert_allclose(result.FBA, 49.408892, rtol=5e-8)
        #
        # The web page reports only that result.pBA < 0.0001
        # Here we'll check the actual p-value.
        #
        #   >>> fdtrc(8, 12, 49.408892)
        #   5.0904479461488075e-08
        #
        expected_pBA = fdtrc(result.DFBA, result.DFerror, result.FBA)
        assert_allclose(result.pBA, expected_pBA, rtol=1e-12)

        if btype == 'random':
            assert_allclose(result.FA, 7.665167, rtol=5e-7)
            assert_allclose(result.pA, 0.0097, rtol=5e-2)
        else:
            # btype == 'fixed'
            assert_allclose(result.FA, 378.727406, rtol=5e-9)
            # The web page reports the corresponding p-value to be < 0.0001
            #
            #   >>> fdtrc(3, 12, 378.727406)
            #   3.804729913731049e-12
            #
            expected_pA_fixed = fdtrc(result.DFA, result.DFerror, result.FA)
            assert_allclose(result.pA, expected_pA_fixed, rtol=1e-12)

    def test_realstatistics_example(self):
        # https://www.real-statistics.com/anova-random-nested-factors/
        #     nested-anova/data-analysis-nested-anova/

        g = [[[34, 56, 78, 67, 61, 48, 72, 51, 52, 60],
              [45, 76, 65, 52, 23, 73, 68, 71, 59, 80],
              [40, 39, 78, 24, 59, 53, 51, 32, 41, 45]],
             [[46, 32, 56, 52, 79, 28, 36, 41, 64, 37],
              [23, 56, 61, 34, 28, 35, 30, 46, 41, 34],
              [75, 57, 63, 25, 60, 52, 45, 48, 38, 26]],
             [[23, 56, 61, 34, 28, 35, 30, 46, 41, 34],
              [75, 57, 63, 25, 60, 52, 45, 48, 38, 26],
              [55, 73, 69, 64, 59, 62, 67, 31, 48, 75]]]

        result = anova.twoway_nested_from_data_grid(g)

        assert_equal(result.DFA, 2)
        assert_equal(result.DFBA, 6)
        assert_equal(result.DFerror, 81)

        assert_allclose(result.SSA, 1559.756, rtol=5e-7)
        assert_allclose(result.SSBA, 4137.133, rtol=5e-7)
        assert_allclose(result.SSerror, 17326.1, rtol=5e-6)

        assert_allclose(result.MSA, 779.8778, rtol=5e-7)
        assert_allclose(result.MSBA, 689.5222, rtol=5e-7)
        assert_allclose(result.MSerror, 213.9025, rtol=5e-7)

        assert_allclose(result.FA, 1.131041, rtol=5e-7)
        assert_allclose(result.FBA, 3.223536, rtol=5e-7)

        assert_allclose(result.pA, 0.382988, rtol=5e-6)
        assert_allclose(result.pBA, 0.006857, rtol=5e-4)

    def test_nist_eng_stats_handbook_example(self):
        # Test the example from
        # https://www.itl.nist.gov/div898/handbook/ppc/section2/ppc233.htm
        # The web page does not show the p-values associated with the
        # F statistics in the ANOVA table.

        data = np.array([[[0.125, 0.118, 0.123, 0.126, 0.118],
                          [0.127, 0.122, 0.125, 0.128, 0.129],
                          [0.125, 0.120, 0.125, 0.126, 0.127],
                          [0.126, 0.124, 0.124, 0.127, 0.120],
                          [0.128, 0.119, 0.126, 0.129, 0.121]],
                         [[0.124, 0.116, 0.122, 0.126, 0.125],
                          [0.128, 0.125, 0.121, 0.129, 0.123],
                          [0.127, 0.119, 0.124, 0.125, 0.114],
                          [0.126, 0.125, 0.126, 0.130, 0.124],
                          [0.129, 0.120, 0.125, 0.124, 0.117]]])

        # Move the last axis to the beginning, so we have the correct
        # layout for the nested data.
        data = np.moveaxis(data, -1, 0)

        result = anova.twoway_nested_from_data_grid(data)

        assert_equal(result.DFA, 4)
        assert_equal(result.DFBA, 5)
        assert_equal(result.DFerror, 40)

        assert_allclose(result.SSA, 3.03e-4, rtol=5e-3)
        assert_allclose(result.SSBA, 1.86e-5, rtol=5e-3)
        assert_allclose(result.SSerror, 3.46e-4, rtol=5e-3)

        assert_allclose(result.MSA, 7.58e-5, rtol=5e-3)
        assert_allclose(result.MSBA, 3.72e-6, rtol=5e-3)
        assert_allclose(result.MSerror, 8.7e-6, rtol=5e-2)

        assert_allclose(result.FA, 20.38, rtol=5e-4)
        # The web page reports 0.428 for the following, but we calculate
        # 0.430, so use a looser tolerance.
        assert_allclose(result.FBA, 0.428, rtol=5e-2)

    @pytest.mark.parametrize('usegrid', [True, False])
    def test_r_balanced_nested_example(self, usegrid):
        # A example of nested ANOVA using R.
        # The data set has two replicates.
        #
        # R code:
        #
        #   x = factor(c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3))
        #   y = factor(c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2))
        #   v = c(2.1, 1.7, 0.5, 0.8, 1.0, 4.0, 2.5, 2.5, 0.2,
        #         0.1, 0.1, 0.2)
        #   av = aov(v ~ x + x:y)
        #   options(digits=10)
        #   summary(av)
        #
        # R output:
        #
        #               Df   Sum Sq  Mean Sq F value   Pr(>F)
        #   x            2 11.05167 5.525833 8.40431 0.018204 *
        #   x:y          3  2.25250 0.750833 1.14195 0.405200
        #   Residuals    6  3.94500 0.657500
        #   ---
        #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
        #

        x = np.array([1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3])
        y = np.array([1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2])
        v = np.array([2.1, 1.7, 0.5, 0.8, 1.0, 4.0, 2.5, 2.5, 0.2,
                      0.1, 0.1, 0.2])
        if usegrid:
            # Form the array with the utility function _nway_groups.
            levels, groups = _nway_groups(x, y, values=v)
            data = np.array(groups.tolist())
            result = anova.twoway_nested_from_data_grid(data, btype='fixed')
        else:
            result = anova.twoway_nested_from_a_b_values(x, y, values=v,
                                                         btype='fixed')

        assert_equal(result.DFA, 2)
        assert_equal(result.DFBA, 3)
        assert_equal(result.DFerror, 6)

        assert_allclose(result.SSA, 11.05167, rtol=5e-7)
        assert_allclose(result.SSBA, 2.25250, rtol=5e-6)
        assert_allclose(result.SSerror, 3.94500, rtol=5e-6)

        assert_allclose(result.MSA, 5.525833, rtol=5e-7)
        assert_allclose(result.MSBA, 0.750833, rtol=5e-6)
        assert_allclose(result.MSerror, 0.657500, rtol=5e-6)

        assert_allclose(result.FA, 8.40431, rtol=5e-6)
        assert_allclose(result.FBA, 1.14195, rtol=5e-6)

        assert_allclose(result.pA, 0.018204, rtol=5e-5)
        assert_allclose(result.pBA, 0.405200, rtol=5e-6)

    def test_rcompanion_nested_random(self):
        # From https://rcompanion.org/rcompanion/d_07.html
        text = """
             Tech  Rat Protein
             Janet 1   1.119
             Janet 1   1.2996
             Janet 1   1.5407
             Janet 1   1.5084
             Janet 1   1.6181
             Janet 1   1.5962
             Janet 1   1.2617
             Janet 1   1.2288
             Janet 1   1.3471
             Janet 1   1.0206
             Janet 2   1.045
             Janet 2   1.1418
             Janet 2   1.2569
             Janet 2   0.6191
             Janet 2   1.4823
             Janet 2   0.8991
             Janet 2   0.8365
             Janet 2   1.2898
             Janet 2   1.1821
             Janet 2   0.9177
             Janet 3   0.9873
             Janet 3   0.9873
             Janet 3   0.8714
             Janet 3   0.9452
             Janet 3   1.1186
             Janet 3   1.2909
             Janet 3   1.1502
             Janet 3   1.1635
             Janet 3   1.151
             Janet 3   0.9367
             Brad  5   1.3883
             Brad  5   1.104
             Brad  5   1.1581
             Brad  5   1.319
             Brad  5   1.1803
             Brad  5   0.8738
             Brad  5   1.387
             Brad  5   1.301
             Brad  5   1.3925
             Brad  5   1.0832
             Brad  6   1.3952
             Brad  6   0.9714
             Brad  6   1.3972
             Brad  6   1.5369
             Brad  6   1.3727
             Brad  6   1.2909
             Brad  6   1.1874
             Brad  6   1.1374
             Brad  6   1.0647
             Brad  6   0.9486
             Brad  7   1.2574
             Brad  7   1.0295
             Brad  7   1.1941
             Brad  7   1.0759
             Brad  7   1.3249
             Brad  7   0.9494
             Brad  7   1.1041
             Brad  7   1.1575
             Brad  7   1.294
             Brad  7   1.4543
        """
        data = np.genfromtxt(text.splitlines(), dtype=None, names=True,
                             encoding='utf8')
        tech = data['Tech']
        rat = data['Rat']
        protein = data['Protein']
        result = anova.twoway_nested_from_a_b_values(tech, rat, protein,
                                                     btype='random')
        # The web page shows how `lme` and `anova.lme` from the nlme library
        # are used to compute the F-value and p-value for `Tech` factor
        # (i.e. the A factor) in the result.  We get the expected results for
        # those terms from the web page.  The other results are the same for
        # btype 'fixed' or 'random', so we get them from this R code:
        #
        #   > summary(aov(Protein ~ Tech + Tech/Rat, data=Data))
        #               Df    Sum Sq    Mean Sq F value    Pr(>F)
        #   Tech         1 0.0384105 0.03841046 1.06585 0.3064875
        #   Tech:Rat     4 0.5739754 0.14349386 3.98179 0.0066601 **
        #   Residuals   54 1.9460277 0.03603755

        assert_equal(result.DFA, 1)
        assert_equal(result.DFBA, 4)
        assert_equal(result.DFerror, 54)

        assert_allclose(result.SSA, 0.03841046, rtol=5e-7)
        assert_allclose(result.SSBA, 0.5739754, rtol=5e-7)
        assert_allclose(result.SSerror, 1.9460277, rtol=5e-7)

        assert_allclose(result.MSA, 0.03841046, rtol=5e-7)
        assert_allclose(result.MSBA, 0.14349386, rtol=5e-8)
        assert_allclose(result.MSerror, 0.03603755, rtol=5e-7)

        assert_allclose(result.FA, 0.2677, rtol=5e-4)  # from web page
        assert_allclose(result.FBA, 3.98179, rtol=5e-6)

        assert_allclose(result.pA, 0.6322, rtol=5e-4)  # from web page
        assert_allclose(result.pBA, 0.0066601, rtol=5e-5)


class TestTwowayNestedUnbalanced:
    """
    Test of two-way nested ANOVA where the number of replicates
    per cell are not all the same.
    """

    # Notes on comparing to results from R:
    #
    # For balanced nested data, with data such as
    #
    #   a <- as.factor(c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3))
    #   b <- as.factor(c(1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3))
    #   y <- c(11, 12, 11, 15, 12, 20,
    #          14, 14, 15, 16, 14, 14,
    #          13, 12, 11, 11, 13, 10)
    #
    # where b is nested in a and b is a fixed effect, one can use
    # `aov(y ~ a + a:b)` for the nested ANOVA, e.g.
    #
    #  > summary(aov(y ~ a + a:b))
    #              Df        Sum Sq        Mean Sq F value  Pr(>F)
    #  a            2 24.7777777778 12.38888888889 2.42391 0.14384
    #  a:b          6 26.3333333333  4.38888888889 0.85870 0.55804
    #  Residuals    9 46.0000000000  5.11111111111
    #
    # and the result will agree with twoway_nested_from_a_b_values, with
    # btype "fixed".  E.g.
    #
    # >>> a = np.array([1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3])
    # >>> b = np.array([1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3])
    # >>> y = np.array([11, 12, 11, 15, 12, 20,
    # ...               14, 14, 15, 16, 14, 14,
    # ...               13, 12, 11, 11, 13, 10])
    # >>> result = anova.twoway_nested_from_a_b_values(a, b, values=y,
    # ...                                              btype='fixed')
    # >>> print(result)
    # ANOVA two-way, factor B nested in factor A, B is fixed
    # Source               SS DF          MS          F          p
    # Factor A    24.77777778  2 12.38888889 2.42391304 0.14383815
    # Factor B(A) 26.33333333  6  4.38888889 0.85869565 0.55804209
    # Error               46.  9  5.11111111
    # Total       97.11111111 17
    #
    # See `TestTwowayNestedBalanced.test_r_balanced_nested_example` for a
    # test that uses a comparison like this.
    #
    # One can find different recommendations for how to do nested ANOVA in R
    # when there are different numbers of replicates per "cell". The
    # computation in anova.two_way_nested_from_a_b_values agrees with this
    # version (as shown at https://www.flutterbys.com.au/stats/tut/tut9.2a.html
    # in Section 4):
    #
    #   contrasts(b) <- contr.sum
    #   library(car)
    #   Anova(aov(y ~ a/b), type='III')
    #
    # See `test_r_unbalanced_nested` in this class for a test where the
    # result of anova.twoway_nested_from_a_b_values is tested against the R
    # function Anova from the car package.

    def test_neter_wasserman_kutner_29_6(self):
        # Example from section 29.6 of
        #
        #   Neter, Wasserman and Kutner, Applied Linear Statistical
        #   Models (second edition), Irwin, Homewood, Illinois, USA (1985)
        #
        # The data is in Table 29.8(a) (page 982), and the results are
        # in Table 29.9 (page 983).  The p-values are not shown in the table.

        # The integers 1, 2, ... are used to encode the levels of factors
        # A and B (instead of "city" and "instructor").
        a = np.array([1, 1, 1, 1, 1, 2, 2, 2, 2])
        b = np.array([1, 1, 2, 3, 3, 1, 1, 2, 2])
        values = [20, 22, 8, 9, 13, 4, 8, 16, 20]

        result = anova.twoway_nested_from_a_b_values(a, b, values,
                                                     btype='fixed')

        assert_equal(result.DFA, 1)
        assert_equal(result.DFBA, 3)
        assert_equal(result.DFerror, 4)

        assert_allclose(result.SSA, 3.76, rtol=5e-3)
        assert_allclose(result.SSBA, 295.20, rtol=5e-5)
        assert_allclose(result.SSerror, 26.0, rtol=5e-3)

        assert_allclose(result.MSA, 3.76, rtol=5e-3)
        assert_allclose(result.MSBA, 98.4, rtol=5e-3)
        assert_allclose(result.MSerror, 6.50, rtol=5e-3)

        assert_allclose(result.FA, 0.58, rtol=5e-2)
        assert_allclose(result.FBA, 15.1, rtol=5e-3)

    def test_r_unbalanced_nested(self):
        # R code (from https://www.flutterbys.com.au/stats/tut/tut9.2a.html,
        # Section 4):
        #
        #   library(car)
        #   a <- as.factor(c(1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3))
        #   b <- as.factor(c(1, 1, 2, 2, 3, 1, 2, 2, 3, 1, 2, 2, 3, 3))
        #   values <- c(10, 10, 11, 11, 10, 12, 13, 10, 15, 11, 12, 14, 15, 15)
        #   options(digits=10)
        #   contrasts(b) <- contr.sum
        #   Anova(aov(values ~ a/b), type='III')
        #
        # R output:
        #
        #   Anova Table (Type III tests)
        #
        #   Response: values
        #              Sum Sq Df   F value     Pr(>F)
        #   (Intercept)  480.5  1 369.61538 7.0214e-06 ***
        #   a             19.5  2   7.50000    0.03125 *
        #   a:b           20.9  6   2.67949    0.14945
        #   Residuals      6.5  5
        #
        a = np.array([1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3])
        b = np.array([1, 1, 2, 2, 3, 1, 2, 2, 3, 1, 2, 2, 3, 3])
        y = np.array([10, 10, 11, 11, 10, 12, 13, 10, 15, 11, 12, 14, 15, 15])

        result = anova.twoway_nested_from_a_b_values(a, b, values=y,
                                                     btype='fixed')

        assert_equal(result.DFA, 2)
        assert_equal(result.DFBA, 6)
        assert_equal(result.DFerror, 5)

        assert_allclose(result.SSA, 19.5, rtol=1e-8)
        assert_allclose(result.SSBA, 20.9, rtol=1e-8)
        assert_allclose(result.SSerror, 6.5, rtol=1e-8)

        assert_allclose(result.MSA, result.SSA/result.DFA, rtol=1e-8)
        assert_allclose(result.MSBA, result.SSBA/result.DFBA, rtol=1e-8)
        assert_allclose(result.MSerror, result.SSerror/result.DFerror,
                        rtol=1e-8)

        assert_allclose(result.FA, 7.5, rtol=1e-8)
        assert_allclose(result.FBA, 2.67949, rtol=5e-6)

        assert_allclose(result.pA, 0.03125, rtol=5e-4)
        assert_allclose(result.pBA, 0.14945, rtol=5e-5)


class TestTwowayFromABValuesValidation:

    @pytest.mark.parametrize('b', [[1, 1, 2, 2, 2], [1, 1, 2, 2]])
    def test_bad_lengths(self, b):
        a = [1, 2, 2, 3]
        values = [10, 11, 12, 13, 14]
        with pytest.raises(ValueError, match='must have the same length'):
            anova.twoway_from_a_b_values(a, b, values, sstype=1)

    def test_bad_type(self):
        a = [1, 1, 1, 2, 2, 2, 3, 3, 3]
        b = [1, 1, 2, 1, 2, 2, 1, 1, 2]
        values = np.arange(len(a))
        with pytest.raises(ValueError, match='`sstype` must be'):
            anova.twoway_from_a_b_values(a, b, values, sstype=99)

    def test_missing_type(self):
        a = [1, 1, 1, 2, 2, 2, 3, 3, 3]
        b = [1, 1, 2, 1, 2, 2, 1, 1, 2]
        values = np.arange(len(a))
        with pytest.raises(ValueError, match='data is unbalanced'):
            anova.twoway_from_a_b_values(a, b, values)

    def test_empty_cell(self):
        # This data has no values for a = 2 and b = 3.
        a = [1, 1, 1, 2, 2, 2, 3, 3, 3]
        b = [1, 2, 3, 1, 2, 2, 1, 2, 3]
        values = np.arange(len(a))
        with pytest.raises(ValueError, match='has no data'):
            anova.twoway_from_a_b_values(a, b, values, sstype=1)


class TestTwowayFromDataGridValidation:

    def test_bad_shape(self):
        data = [[[1, 2], [1, 2, 3], [1, 3]],
                [[1, 1], [2, 2, 3]]]
        with pytest.raises(ValueError, match='must have the same length'):
            anova.twoway_from_data_grid(data, sstype=3)

    def test_bad_type(self):
        data = [[[1, 2], [3, 5]],
                [[2, 9], [1, 2, 3]]]
        with pytest.raises(ValueError, match='`sstype` must be'):
            anova.twoway_from_data_grid(data, sstype=99)

    def test_missing_type(self):
        data = [[[1, 2], [3, 5]],
                [[2, 9], [1, 2, 3]]]
        with pytest.raises(ValueError, match='data is unbalanced'):
            anova.twoway_from_data_grid(data)

    def test_empty_cell(self):
        data = [[[1, 2], [3, 5]],
                [[2, 9], []]]
        with pytest.raises(ValueError, match='has no data'):
            anova.twoway_from_data_grid(data, sstype=1)
