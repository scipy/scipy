"""
===================================================================
Statistical functions for masked arrays (:mod:`scipy.stats.mstats`)
===================================================================

.. currentmodule:: scipy.stats.mstats

.. deprecated:: 2.0.0

   This module is deprecated and will be removed in SciPy 2.4.0.
   See function documentation for alternatives.

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
from scipy import stats
from scipy._lib._array_api import xp_capabilities_table
from scipy._lib._docscrape import FunctionDoc
from . import _mstats_basic
from . import _mstats_extras
from ._mstats_basic import *  # noqa: F403
from ._mstats_extras import *  # noqa: F403
# Functions that support masked array input in stats but need to be kept in the
# mstats namespace for backwards compatibility:
from scipy.stats import gmean, hmean, zmap, zscore, chisquare

__all__ = _mstats_basic.__all__ + _mstats_extras.__all__
__all__ += ['gmean', 'hmean', 'zmap', 'zscore', 'chisquare']


_trim_transition_guide = ("See the Trimming and winsorization transition guide "
                          "for alternatives.")


_mstats_deprecation_table = {
    'hdquantiles': dict(replacement='quantile'),
    'hdmedian': dict(replacement='quantile'),
    'mquantiles': dict(replacement='quantile'),
    'scoreatpercentile': dict(replacement='quantile'),
    'plotting_positions': dict(replacement='estimated_cdf'),
    'trima': dict(notes=_trim_transition_guide),
    'trimr': dict(notes=_trim_transition_guide),
    'trim': dict(notes=_trim_transition_guide),
    'trimboth': dict(notes=_trim_transition_guide),
    'trimtail': dict(notes=_trim_transition_guide),
    'trimmed_mean': dict(notes=_trim_transition_guide),
    'trimmed_mean_ci': dict(notes=_trim_transition_guide),
    'trimmed_std': dict(notes=_trim_transition_guide),
    'trimmed_var': dict(notes=_trim_transition_guide),
    'trimmed_stde': dict(notes=_trim_transition_guide),
    'winsorize': dict(notes=_trim_transition_guide),
}


def _get_deprecation_message(fun, _mstats_deprecation_table, stats):
    # generate a deprecation message based on whether there is a stats replacement,
    # and whether that has MArray or only nan_policy='omit' support
    # TODO: clean this up a bit
    if fun.__name__ not in _mstats_deprecation_table:
        _mstats_deprecation_table[fun.__name__]= {}
    _deprecation_entry = _mstats_deprecation_table[fun.__name__]
    notes = _deprecation_entry.get("notes", "")
    replacement = _deprecation_entry.get('replacement', None)
    replacement = fun.__name__ if replacement is None else replacement
    replacement_fun = getattr(stats, str(replacement), False)
    replaced = bool(replacement_fun)
    marray = xp_capabilities_table[replacement_fun]['marray'] if replaced else False

    message = (f"`scipy.stats.mstats.{fun.__name__}` is deprecated as of "
               "SciPy 2.0.0 and will be removed, along with the "
               "`scipy.stats.mstats` namespace, in SciPy 2.4.0. ")
    legacy_message = ("If needed, the legacy `scipy.stats.mstats` from SciPy "
                      "1.18.0 can be installed as the package `mstats`.")

    if not replaced:
        message += "SciPy offers no replacement for this function. "
    else:
        message += "For similar functionality, "
        if marray:
            message += ("use MArray(s) instead of NumPy masked array(s) with "
                        f"`scipy.stats.{replacement}`. ")
        else:
            message += (f"use `scipy.stats.{replacement}` with regular NumPy "
                        "arrays, replacing masked values with NaNs and using "
                        "the `nan_policy='omit'` option. ")

    message = message + f"{notes} " if notes else message

    message += legacy_message
    return message


def _document_deprecation(wrapper, message):
    # append deprecation admonition to extended summary of wrapper
    doc = FunctionDoc(wrapper)
    doc['Extended Summary'].append("")
    doc['Extended Summary'].append(".. deprecated:: 2.0.0")
    doc['Extended Summary'].append("")
    doc['Extended Summary'].append(f"    {message}")
    doc = str(doc).split("\n", 1)[1].lstrip(" \n")  # remove signature
    wrapper.__doc__ = str(doc)


# modify all public `mstats` functions to include deprecation
# admonition in documentation and warn when executed
for fun_name in __all__:
   fun = globals()[fun_name]
   message =_get_deprecation_message(fun, _mstats_deprecation_table, stats)
   wrapper = warnings.deprecated(message)(fun)
   _document_deprecation(wrapper, message)
   globals()[fun_name] = wrapper

# warn when `scipy.stats.mstats` is intentionally imported
# (`scipy.stats.__getattr__` is overridden to enable lazy import)
msg = ("`scipy.stats.mstats` is deprecated as of SciPy 2.0.0 and will be removed "
       "in SciPy 2.4.0. See function documentation for alternatives.")
warnings.warn(msg, DeprecationWarning, stacklevel=4)

# clean up
del warnings
del stats
del _trim_transition_guide
del _mstats_deprecation_table
del _document_deprecation
del _get_deprecation_message
