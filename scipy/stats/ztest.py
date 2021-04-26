# Copyright 2002 Gary Strangman.  All rights reserved
# Copyright 2002-2016 The SciPy Developers
#
# The original code from Gary Strangman was heavily adapted for
# use in SciPy by Travis Oliphant.  The original code came with the
# following disclaimer:
#
# This software is provided "as-is".  There are no expressed or implied
# warranties of any kind, including, but not limited to, the warranties
# of merchantability and fitness for a given application.  In no event
# shall Gary Strangman be liable for any direct, indirect, incidental,
# special, exemplary or consequential damages (including, but not limited
# to, loss of use, data or profits, or business interruption) however
# caused and on any theory of liability, whether in contract, strict
# liability or tort (including negligence or otherwise) arising in any way
# out of the use of this software, even if advised of the possibility of
# such damage.

"""
A collection of basic statistical tests for Python to run a Z-test.

References
----------
.. [CRCProbStat2000] Zwillinger, D. and Kokoska, S. (2000). CRC Standard
   Probability and Statistics Tables and Formulae. Chapman & Hall: New
   York. 2000.

"""
import warnings
import math
from math import gcd
from collections import namedtuple
from itertools import permutations

import numpy as np
from numpy import array, asarray, ma

from scipy.spatial.distance import cdist
from scipy.ndimage import measurements
from scipy._lib._util import (check_random_state, MapWrapper,
                              rng_integers, float_factorial)
import scipy.special as special
from scipy import linalg
from . import distributions
from . import mstats_basic
from ._stats_mstats_common import (_find_repeats, linregress, theilslopes,
                                   siegelslopes)
from ._stats import (_kendall_dis, _toint64, _weightedrankedtau,
                     _local_correlations)
from dataclasses import make_dataclass


# Functions/classes in other files should be added in `__init__.py`, not here
__all__ = ['find_repeats', 'gmean', 'hmean', 'mode', 'tmean', 'tvar',
           'tmin', 'tmax', 'tstd', 'tsem', 'moment', 'variation',
           'skew', 'kurtosis', 'describe', 'skewtest', 'kurtosistest',
           'normaltest', 'jarque_bera', 'itemfreq',
           'scoreatpercentile', 'percentileofscore',
           'cumfreq', 'relfreq', 'obrientransform',
           'sem', 'zmap', 'zscore', 'iqr', 'gstd', 'median_absolute_deviation',
           'median_abs_deviation',
           'sigmaclip', 'trimboth', 'trim1', 'trim_mean',
           'f_oneway', 'F_onewayConstantInputWarning',
           'F_onewayBadInputSizesWarning',
           'PearsonRConstantInputWarning', 'PearsonRNearConstantInputWarning',
           'pearsonr', 'fisher_exact',
           'SpearmanRConstantInputWarning', 'spearmanr', 'pointbiserialr',
           'kendalltau', 'weightedtau', 'multiscale_graphcorr',
           'linregress', 'siegelslopes', 'theilslopes', 'ttest_1samp',
           'ttest_ind', 'ttest_ind_from_stats', 'ttest_rel',
           'kstest', 'ks_1samp', 'ks_2samp',
           'chisquare', 'power_divergence', 'mannwhitneyu',
           'tiecorrect', 'ranksums', 'kruskal', 'friedmanchisquare',
           'rankdata',
           'combine_pvalues', 'wasserstein_distance', 'energy_distance',
           'brunnermunzel', 'alexandergovern']
_one_sample_z_table = {
    '0.00': 0.5000, '0.01': 0.5040, '0.02': 0.5080, '0.03': 0.5120, '0.04': 0.5160, '0.05': 0.5199, '0.06': 0.5239, '0.07': 0.5279, '0.08': 0.5319, '0.09': 0.5359,
    '0.10': 0.5398, '0.11': 0.5438, '0.12': 0.5478, '0.13': 0.5517, '0.14': 0.5557, '0.15': 0.5596, '0.16': 0.5636, '0.17': 0.5675, '0.18': 0.5714, '0.19': 0.5753,
    '0.20': 0.5793, '0.21': 0.5832, '0.22': 0.5871, '0.23': 0.5910, '0.24': 0.5948, '0.25': 0.5987, '0.26': 0.6026, '0.27': 0.6064, '0.28': 0.6103, '0.29': 0.6141,
    '0.30': 0.6179, '0.31': 0.6217, '0.32': 0.6255, '0.33': 0.6293, '0.34': 0.6331, '0.35': 0.6368, '0.36': 0.6406, '0.37': 0.6443, '0.38': 0.6480, '0.39': 0.6517,
    '0.40': 0.6554, '0.41': 0.6591, '0.42': 0.6628, '0.43': 0.6664, '0.44': 0.6700, '0.45': 0.6736, '0.46': 0.6772, '0.47': 0.6808, '0.48': 0.6844, '0.49': 0.6879,
    '0.50': 0.6915, '0.51': 0.6950, '0.52': 0.6985, '0.53': 0.7019, '0.54': 0.7054, '0.55': 0.7088, '0.56': 0.7123, '0.57': 0.7157, '0.58': 0.7190, '0.59': 0.7224,
    '0.60': 0.7257, '0.61': 0.7291, '0.62': 0.7324, '0.63': 0.7357, '0.64': 0.7389, '0.65': 0.7422, '0.66': 0.7454, '0.67': 0.7486, '0.68': 0.7517, '0.69': 0.7549,
    '0.70': 0.7580, '0.71': 0.7611, '0.72': 0.7642, '0.73': 0.7673, '0.74': 0.7704, '0.75': 0.7734, '0.76': 0.7764, '0.77': 0.7794, '0.78': 0.7823, '0.79': 0.7852,
    '0.80': 0.7881, '0.81': 0.7910, '0.82': 0.7939, '0.83': 0.7967, '0.84': 0.7995, '0.85': 0.8023, '0.86': 0.8051, '0.87': 0.8078, '0.88': 0.8106, '0.89': 0.8133,
    '0.90': 0.8159, '0.91': 0.8186, '0.92': 0.8212, '0.93': 0.8238, '0.94': 0.8364, '0.95': 0.8289, '0.96': 0.8315, '0.97': 0.8340, '0.98': 0.8365, '0.99': 0.8389,
    '1.00': 0.8413, '1.01': 0.8438, '1.02': 0.8461, '1.03': 0.8485, '1.04': 0.8508, '1.05': 0.8531, '1.06': 0.8554, '1.07': 0.8577, '1.08': 0.8599, '1.09': 0.8621,
    '1.10': 0.8643, '1.11': 0.8665, '1.12': 0.8686, '1.13': 0.8708, '1.14': 0.8729, '1.15': 0.8749, '1.16': 0.8770, '1.17': 0.8790, '1.18': 0.8810, '1.19': 0.8830,
    '1.20': 0.8849, '1.21': 0.8869, '1.22': 0.8888, '1.23': 0.8907, '1.24': 0.8925, '1.25': 0.8944, '1.26': 0.8962, '1.27': 0.8980, '1.28': 0.8997, '1.29': 0.9015,
    '1.30': 0.9032, '1.31': 0.9049, '1.32': 0.9066, '1.33': 0.9082, '1.34': 0.9099, '1.35': 0.9115, '1.36': 0.9131, '1.37': 0.9147, '1.38': 0.9162, '1.39': 0.9177,
    '1.40': 0.9192, '1.41': 0.9207, '1.42': 0.9222, '1.43': 0.9236, '1.44': 0.9251, '1.45': 0.9265, '1.46': 0.9279, '1.47': 0.9292, '1.48': 0.9306, '1.49': 0.9319,
    '1.50': 0.9332, '1.51': 0.9345, '1.52': 0.9357, '1.53': 0.9370, '1.54': 0.9382, '1.55': 0.9394, '1.56': 0.9406, '1.57': 0.9418, '1.58': 0.9429, '1.59': 0.9441,
    '1.60': 0.9452, '1.61': 0.9463, '1.62': 0.9474, '1.63': 0.9484, '1.64': 0.9495, '1.65': 0.9505, '1.66': 0.9515, '1.67': 0.9525, '1.68': 0.9539, '1.69': 0.9545,
    '1.70': 0.9554, '1.71': 0.9564, '1.72': 0.9573, '1.73': 0.9582, '1.74': 0.9591, '1.75': 0.9599, '1.76': 0.9608, '1.77': 0.9616, '1.78': 0.9625, '1.79': 0.9633,
    '1.80': 0.9641, '1.81': 0.9649, '1.82': 0.9656, '1.83': 0.9664, '1.84': 0.9671, '1.85': 0.9678, '1.86': 0.9696, '1.87': 0.9693, '1.88': 0.9699, '1.89': 0.9706,
    '1.90': 0.9713, '1.91': 0.9719, '1.92': 0.9726, '1.93': 0.9732, '1.94': 0.9738, '1.95': 0.9744, '1.96': 0.9750, '1.97': 0.9756, '1.98': 0.9761, '1.99': 0.9767,
    '2.00': 0.9772, '2.01': 0.9778, '2.02': 0.9783, '2.03': 0.9788, '2.04': 0.9793, '2.05': 0.9798, '2.06': 0.9803, '2.07': 0.9808, '2.08': 0.9812, '2.09': 0.9817,
    '2.10': 0.9821, '2.11': 0.9826, '2.12': 0.9830, '2.13': 0.9834, '2.14': 0.9838, '2.15': 0.9842, '2.16': 0.9846, '2.17': 0.9850, '2.18': 0.9854, '2.19': 0.9857,
    '2.20': 0.9861, '2.21': 0.9864, '2.22': 0.9868, '2.23': 0.9871, '2.24': 0.9875, '2.25': 0.9878, '2.26': 0.9881, '2.27': 0.9884, '2.28': 0.9887, '2.29': 0.9890,
    '2.30': 0.9893, '2.31': 0.9896, '2.32': 0.9898, '2.33': 0.9901, '2.34': 0.9904, '2.35': 0.9906, '2.36': 0.9909, '2.37': 0.9911, '2.38': 0.9913, '2.39': 0.9916,
    '2.40': 0.9918, '2.41': 0.9920, '2.42': 0.9922, '2.43': 0.9925, '2.44': 0.9927, '2.45': 0.9929, '2.46': 0.9931, '2.47': 0.9932, '2.48': 0.9934, '2.49': 0.9936,
    '2.50': 0.9938, '2.51': 0.9940, '2.52': 0.9941, '2.53': 0.9943, '2.54': 0.9945, '2.55': 0.9946, '2.56': 0.9948, '2.57': 0.9949, '2.58': 0.9951, '2.59': 0.9952,
    '2.60': 0.9953, '2.61': 0.9955, '2.62': 0.9956, '2.63': 0.9957, '2.64': 0.9959, '2.65': 0.9960, '2.66': 0.9961, '2.67': 0.9962, '2.68': 0.9963, '2.69': 0.9964,
    '2.70': 0.9965, '2.71': 0.9966, '2.72': 0.9967, '2.73': 0.9968, '2.74': 0.9969, '2.75': 0.9970, '2.76': 0.9971, '2.77': 0.9972, '2.78': 0.9973, '2.79': 0.9974,
    '2.80': 0.9974, '2.81': 0.9975, '2.82': 0.9976, '2.83': 0.9977, '2.84': 0.9977, '2.85': 0.9978, '2.86': 0.9979, '2.87': 0.9979, '2.88': 0.9980, '2.89': 0.9981,
    '2.90': 0.9981, '2.91': 0.9982, '2.92': 0.9982, '2.93': 0.9983, '2.94': 0.9984, '2.95': 0.9984, '2.96': 0.9985, '2.97': 0.9985, '2.98': 0.9986, '2.99': 0.9986,
    '3.00': 0.9987, '3.01': 0.9987, '3.02': 0.9987, '3.03': 0.9988, '3.04': 0.9988, '3.05': 0.9989, '3.06': 0.9989, '3.07': 0.9989, '3.08': 0.9990, '3.09': 0.9990,
    '3.10': 0.9990, '3.11': 0.9991, '3.12': 0.9991, '3.13': 0.9991, '3.14': 0.9992, '3.15': 0.9992, '3.16': 0.9992, '3.17': 0.9992, '3.18': 0.9993, '3.19': 0.9993,
    '3.20': 0.9993, '3.21': 0.9993, '3.22': 0.9994, '3.23': 0.9994, '3.24': 0.9994, '3.25': 0.9994, '3.26': 0.9994, '3.27': 0.9995, '3.28': 0.9995, '3.29': 0.9995,
    '3.30': 0.9995, '3.31': 0.9995, '3.32': 0.9995, '3.33': 0.9996, '3.34': 0.9996, '3.35': 0.9996, '3.36': 0.9996, '3.37': 0.9996, '3.38': 0.9996, '3.39': 0.9997,
    '3.40': 0.9997, '3.41': 0.9997, '3.42': 0.9997, '3.43': 0.9997, '3.44': 0.9997, '3.45': 0.9997, '3.46': 0.9997, '3.47': 0.9997, '3.48': 0.9997, '3.49': 0.9998
}

_two_sample_z_values = {
    '0.10': 1.645,
    '0.05': 1.96,
    '0.02': 2.326,
    '0.01': 2.576
}



def _one_sample_z(N, sample_mean, target_mean, stdev_sample, op, alpha):
    """One sample Z test return 1 if reject null."""
    # How to run a one sample Z test: 
    # State the null hypothesis and alternative hypothesis 
    # Find the critical value of Z in the Z table 
    # Calculate the Z test statistic 
    # Compare the test statistic to the critical z value and decide to support or reject the null hypothesis

    # Z-FORMULA to find a z-score
    # Z = X(Sample mean)  - Uo (mean for hypothesis test) / (O (standard deviation) /radical(N) (# items in sample))
    # need to take in: 
        # num in sample 
        # Sample Mean 
        # Test Mean 
        # Sample Standard Deviation
    # Need an interface to do this

    # Op comes from drop down options: 
        # H0: μ≥μ0 --> op = 'Right'
        # H0: μ≤μ0 --> op = 'Left'
        # H0: μ=μ0 --> op = 'Two'

    Z_value_numerator = sample_mean -  target_mean
    Z_value_denominator = stdev_sample / (math.sqrt(N))
    Z_test_val = Z_value_numerator / Z_value_denominator

    # After Calculating this Z-value we need to look up the score that we get in the Z-table 
    rounded_z_val = round(Z_test_val, 2)
    if rounded_z_val < 0: 
        rounded_z_val = rounded_z_val* -1
    if rounded_z_val > 3.49: 
        rounded_z_val = 3.49
    
    area_val = _one_sample_z_table[rounded_z_val]
    tail_value = 1 - area_val
    final_val = 0
    # Now need to see if right tailed test, left tailed test, or two tailed test 
    # OP should take in a value selected from dropdown indicating: 
        # Two Tailed = "Two"
        # Right Tailed = "Right"
        # Left Tailed = "Left"
    if op == 'Two'
        final_val = tail_value * 2
    else: 
        final_val = tail_value
    
    if final_val < alpha:
        return 1
    else: 
        return 0

def _two_independent_sample_z(P1_positive, P1_total, P2_positive, P2_total, alpha):
    """Two independent sample Z test return 1 if reject null."""
    # The null hypothesis is always that the population means are equal --> No options for null hypothesis here
    # Variables that need to be input: 
        # Sample 1 Positive outcomes 
        # Sample 1 population size
        # Sample 2 Positive Outcomes
        # Sample 2 Population Size
        # Confidence Level
        # Alpha
    population1_proportion = P1_positive/P1_total
    population2_proportion = P2_positive/P2_total
    overall_sample_proportion = (P1_positive + P2_positive)/(P1_total + P2_total)

    # Breaking down formula into more readable features
    Z_stat_numerator = abs(population1_proportion - population2_proportion)
    Z_stat_rad1 = overall_sample_proportion * (1 - overall_sample_proportion)
    Z_stat_rad2 = (1/P1_total) + (1/P2_total)
    Z_stat_denominator = sqrt(Z_stat_rad1 * Z_stat_rad2)

    calculated_Z_stat = Z_stat_numerator/Z_stat_denominator
    table_Z_stat = _two_sample_z_values[alpha]
    if calculated_Z_stat > table_Z_stat: 
        return 1
    else: 
        return 0


def _two_dependent_sample_z(N, sample_mean, sample_SD, op, alpha):
    """Two dependent sample Z test return 1 if reject null."""
    # The null hypothesis is always either the population means are not equal, greater than, or less than
    # Variables that Need to be input for this test: 
        # Total sample population size: N 
        # Sample Mean of differences: sample_mean
        # Sample Standard Deviation: sample_SD
        # Null hypothesis Value: (This will be op)
            # µ = op --> op value represents the difference in the two means
    standard_error_mean = sample_SD / sqrt(N)

    # Compute the test statistic
    test_statistic_numerator = sample_mean - op
    test_statistic = test_statistic_numerator / standard_error_mean

    upper_bound_CI = _two_sample_z_values[alpha]
    lower_bound_CI = -1 * upper_bound_CI

    # If the value that we calculated for test stat is within interval, accept Null
    if test_statistic >= lower_bound_CI and test_statistic <= upper_bound_CI: 
        return 0
    else: 
        return 1