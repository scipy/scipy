import functools
import warnings
import numpy as np
from .._lib._bunch import _make_tuple_bunch
from scipy._lib._array_api import xp_capabilities_table
from scipy import stats

__all__ = ['_find_repeats']

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
    """Decorator to deprecate mstats function

    """
    def decorator(fun):
        replacement_ = fun.__name__ if replacement is None else replacement

        @functools.wraps(fun)
        def wrapper(*args, **kwargs):
            replacement_fun = getattr(stats, str(replacement_), False)
            replaced = bool(replacement_fun)
            marray = xp_capabilities_table[replacement_fun]['marray'] if replaced else False

            message = (f"`scipy.stats.mstats.{fun.__name__}` is deprecated as of "
                       "SciPy 2.0.0 and will be removed, along with the "
                       "`scipy.stats.mstats` namespace, in SciPy 2.2.0. ") 
            legacy_message = ("If needed, the legacy `scipy.stats.mstats` from SciPy "
                              "1.18.0 can be installed as the package `scipy_mstats`.")
            
            if not replaced:
                message += "SciPy offers no replacement for this function. "
            else:
                message += "For similar functionality, "
                if marray:
                    message += ("use MArray(s) instead of NumPy masked array(s) with "
                                f"`scipy.stats.{replacement_}`. ")
                else:
                    message += (f"use `scipy.stats.{replacement_}` with regular NumPy "
                                "arrays, replacing masked values with NaNs and using "
                                "the `nan_policy='omit' option. ")  
                    
            message = message + f"{notes} " if notes else message

            message += legacy_message
            warnings.warn(message, DeprecationWarning, stacklevel=2)
            return fun(*args, **kwargs)
        return wrapper

        # if replace_doc:
        #     doc = FunctionDoc(wrapper)
        #     parameter_names = [param.name for param in doc['Parameters']]
        #     if 'rng' in parameter_names:
        #         _type = "{None, int, `numpy.random.Generator`}, optional"
        #         _desc = _rng_desc.replace("{old_name}", old_name)
        #         old_doc = doc['Parameters'][parameter_names.index('rng')].desc
        #         old_doc_keep = old_doc[old_doc.index("") + 1:] if "" in old_doc else []
        #         new_doc = [_desc] + old_doc_keep
        #         _rng_parameter_doc = Parameter('rng', _type, new_doc)
        #         doc['Parameters'][parameter_names.index('rng')] = _rng_parameter_doc
        #         doc = str(doc).split("\n", 1)[1].lstrip(" \n")  # remove signature
        #         wrapper.__doc__ = str(doc)
        # return wrapper

    return decorator
