"""
===================================================================
Statistical functions for masked arrays (:mod:`scipy.stats.mstats`)
===================================================================

.. currentmodule:: scipy.stats.mstats

This module contains a large number of statistical functions that can
be used with masked arrays.

Most of these functions are similar to those in scipy.stats but might
have small differences in the API or in the algorithm used. Since this
is a relatively new package, some API changes are still possible.

.. autosummary::
   :toctree: generated/

   argstoarray
   betai
   chisquare
   count_tied_groups
   describe
   f_oneway
   f_value_wilks_lambda
   find_repeats
   friedmanchisquare
   kendalltau
   kendalltau_seasonal
   kruskalwallis
   ks_twosamp
   kurtosis
   kurtosistest
   linregress
   mannwhitneyu
   plotting_positions
   mode
   moment
   mquantiles
   msign
   normaltest
   obrientransform
   pearsonr
   plotting_positions
   pointbiserialr
   rankdata
   scoreatpercentile
   sem
   signaltonoise
   skew
   skewtest
   spearmanr
   theilslopes
   threshold
   tmax
   tmean
   tmin
   trim
   trima
   trimboth
   trimmed_stde
   trimr
   trimtail
   tsem
   ttest_onesamp
   ttest_ind
   ttest_onesamp
   ttest_rel
   tvar
   variation
   winsorize
   zmap
   zscore

"""
from __future__ import division, print_function, absolute_import

from .mstats_basic import *
from .mstats_extras import *
from scipy.stats import gmean, hmean
