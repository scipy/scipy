"""
===================================================================
Statistical functions for masked arrays (:mod:`scipy.stats.mstats`)
===================================================================

.. currentmodule:: scipy.stats.mstats

This module contains a large number of statistical functions that can
be used with masked arrays.

Most of these functions are similar to those in `scipy.stats` but might
have small differences in the API or in the algorithm used. Since this
is a relatively new package, some API changes are still possible.

Summary statistics
==================

.. autosummary::
   :toctree: generated/

   describe
   gmean
   hmean
   kurtosis
   mode
   mquantiles
   hdmedian
   hdquantiles
   hdquantiles_sd
   idealfourths
   plotting_positions
   meppf
   moment
   skew
   tmean
   tvar
   tmin
   tmax
   tsem
   variation
   find_repeats
   sem
   trimmed_mean
   trimmed_mean_ci
   trimmed_std
   trimmed_var

Frequency statistics
====================

.. autosummary::
   :toctree: generated/

   scoreatpercentile

Correlation functions
=====================

.. autosummary::
   :toctree: generated/

   f_oneway
   pearsonr
   spearmanr
   pointbiserialr
   kendalltau
   kendalltau_seasonal
   linregress
   siegelslopes
   theilslopes
   sen_seasonal_slopes

Statistical tests
=================

.. autosummary::
   :toctree: generated/

   ttest_1samp
   ttest_onesamp
   ttest_ind
   ttest_rel
   chisquare
   kstest
   ks_2samp
   ks_1samp
   ks_twosamp
   mannwhitneyu
   rankdata
   kruskal
   kruskalwallis
   friedmanchisquare
   brunnermunzel
   skewtest
   kurtosistest
   normaltest

Transformations
===============

.. autosummary::
   :toctree: generated/

   obrientransform
   trim
   trima
   trimmed_stde
   trimr
   trimtail
   trimboth
   winsorize
   zmap
   zscore

Other
=====

.. autosummary::
   :toctree: generated/

   argstoarray
   count_tied_groups
   msign
   compare_medians_ms
   median_cihs
   mjci
   mquantiles_cimj
   rsh

"""
import warnings
from . import _mstats_basic
from . import _mstats_extras
from ._mstats_basic import *  # noqa: F403
from ._mstats_extras import *  # noqa: F403
from ._stats_mstats_common import (_mstats_deprecation_table,
                                   _generate_deprecation_message)
# Functions that support masked array input in stats but need to be kept in the
# mstats namespace for backwards compatibility:
from scipy.stats import gmean, hmean, zmap, zscore, chisquare

__all__ = _mstats_basic.__all__ + _mstats_extras.__all__
__all__ += ['gmean', 'hmean', 'zmap', 'zscore', 'chisquare']

# process the deprecation messages now that `stats` is initialized
for fun in _mstats_deprecation_table:
   _generate_deprecation_message(fun, globals()[fun.__name__])
del _mstats_deprecation_table
del _generate_deprecation_message

msg = ("`scipy.stats.mstats` is deprecated as of SciPy 2.0.0 and will be removed "
       "in SciPy 2.4.0. See function documentation for alternatives.")
warnings.warn(msg, DeprecationWarning, stacklevel=2)
del warnings
