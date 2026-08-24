import functools
import warnings
import numpy as np
from .._lib._bunch import _make_tuple_bunch
from scipy._lib._array_api import xp_capabilities_table
from scipy._lib._docscrape import FunctionDoc
from scipy import stats

__all__ = ['_find_repeats']

_mstats_deprecation_table = {}

# This is not a namedtuple for backwards compatibility. See PR #12983
TheilslopesResult = _make_tuple_bunch('TheilslopesResult',
                                      ['slope', 'intercept',
                                       'low_slope', 'high_slope'])
SiegelslopesResult = _make_tuple_bunch('SiegelslopesResult',
                                       ['slope', 'intercept'])


def _n_samples_optional_x(kwargs):
    return 2 if kwargs.get('x', None) is not None else 1


def _find_repeats(arr):
    # This function assumes it may clobber its input.
    if len(arr) == 0:
        return np.array(0, np.float64), np.array(0, np.intp)

    # XXX This cast was previously needed for the Fortran implementation,
    # should we ditch it?
    arr = np.asarray(arr, np.float64).ravel()
    arr.sort()

    # Taken from NumPy 1.9's np.unique.
    change = np.concatenate(([True], arr[1:] != arr[:-1]))
    unique = arr[change]
    change_idx = np.concatenate(np.nonzero(change) + ([arr.size],))
    freq = np.diff(change_idx)
    atleast2 = freq > 1
    return unique[atleast2], freq[atleast2]


# mstats deprecation
def _deprecate_mstats(replacement=None, notes=""):
    """Decorator to deprecate mstats function"""
    def decorator(fun):
        # Can't produce the deprecation message until stats module is initialized,
        # so store the args in a table to be processed after initialization.
        _mstats_deprecation_table[fun] = {'replacement': replacement, 'notes': notes}

        @functools.wraps(fun)
        def wrapper(*args, **kwargs):
            warnings.warn(_mstats_deprecation_table[fun]['message'],
                          DeprecationWarning, stacklevel=2)
            return fun(*args, **kwargs)
        return wrapper
    return decorator


def _generate_deprecation_message(fun, wrapper):
    notes = _mstats_deprecation_table[fun]['notes']
    replacement = _mstats_deprecation_table[fun]['replacement']
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
                        "the `nan_policy='omit' option. ")

    message = message + f"{notes} " if notes else message

    message += legacy_message

    doc = FunctionDoc(wrapper)
    doc['Extended Summary'].append("")
    doc['Extended Summary'].append(".. deprecated:: 2.0.0")
    doc['Extended Summary'].append("")
    doc['Extended Summary'].append(f"    {message}")
    doc = str(doc).split("\n", 1)[1].lstrip(" \n")  # remove signature
    wrapper.__doc__ = str(doc)
    _mstats_deprecation_table[fun]['message'] = message
