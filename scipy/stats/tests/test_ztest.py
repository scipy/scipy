from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (assert_, assert_equal, assert_allclose,
                           assert_almost_equal)  # avoid new uses

import pytest
from pytest import raises as assert_raises
from scipy.stats._hypotests import (epps_singleton_2samp, cramervonmises,
                                    _cdf_cvm, cramervonmises_2samp,
                                    _pval_cvm_2samp_exact, barnard_exact)
import scipy.stats as stats
from scipy.stats import distributions
from .common_tests import check_named_results
import ztest.py as ztest


class TestZTests:
    def test_one_sample_z(self):
        # Test one Sample z-test

        # 1,500 women followed the Atkin’s diet for a month. 
        # A random sample of 29 women gained an average of 6.7 pounds. 
        # Test the hypothesis that the average weight gain per woman for the month was over 5 pounds. 
        # The standard deviation for all women in the group was 7.1.
        # The alpha will be 0.10
        N = 29 
        sample_mean = 6.7 
        target_mean = 5.0
        stdev_sample = 7.1
        op = 'Right'
        alpha = 0.10

        # _one_sample_z(N, sample_mean, target_mean, stdev_sample, op, alpha):
        test_result = ztest._one_sample_z(N, sample_mean, target_mean, stdev_sample, op, alpha)
        assert test_result == 0

    def test_two_independent_sample_z(self):
        # Test two independent sample Z test
        # _two_independent_sample_z(P1_positive, P1_total, P2_positive, P2_total, alpha):
        P1_positive = 20
        P1_total = 100
        P2_positive = 9
        P2_total = 36
        alpha = 10
        #Z statistic = (0.05)
        #overall sample proportion = 0.2132
        # z denominator = 0.00633684
        # z stat = 7.89

        test_result = ztest._two_independent_sample_z(P1_positive, P1_total, P2_positive, P2_total, alpha)
        assert test_result == 0


    def test_two_dependent_sample_z(self):
        # Test two independent sample Z test
        # _two_dependent_sample_z(N, sample_mean, sample_SD, op, alpha):
        N = 100
        sample_mean = -219
        sample_SD = 725
        op = 0
        alpha = 0.05

        test_result = ztest._two_dependent_sample_z(N, sample_mean, sample_SD, op, alpha)
        assert test_result == 1