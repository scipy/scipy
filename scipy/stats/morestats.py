# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.stats` namespace for importing the functions
# included below.

import scipy
import warnings
from . import _morestats


__all__ = [  # noqa: F822
    'mvsdist',
    'bayes_mvs', 'kstat', 'kstatvar', 'probplot', 'ppcc_max', 'ppcc_plot',
    'boxcox_llf', 'boxcox', 'boxcox_normmax', 'boxcox_normplot',
    'shapiro', 'anderson', 'ansari', 'bartlett', 'levene',
    'fligner', 'mood', 'wilcoxon', 'median_test',
    'circmean', 'circvar', 'circstd', 'anderson_ksamp',
    'yeojohnson_llf', 'yeojohnson', 'yeojohnson_normmax',
    'yeojohnson_normplot', 'annotations', 'namedtuple', 'isscalar', 'log',
    'around', 'unique', 'arange', 'sort', 'amin', 'amax', 'atleast_1d',
    'array', 'compress', 'exp', 'ravel', 'count_nonzero', 'arctan2',
    'hypot', 'optimize', 'find_repeats',
    'chi2_contingency', 'distributions', 'rv_generic', 'Mean',
    'Variance', 'Std_dev', 'ShapiroResult', 'AndersonResult',
    'Anderson_ksampResult', 'AnsariResult', 'BartlettResult',
    'LeveneResult', 'FlignerResult', 'WilcoxonResult'
]


def __dir__():
    return __all__


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(
            f"`scipy.stats.morestats` has no attribute `{name}`; furthermore, "
            "`scipy.stats.morestats` is deprecated and will be removed in "
            "SciPy 2.0.0.")

    attr = getattr(scipy.stats, name, None)

    if attr is not None:
        message = (f"Please import `{name}` from the `scipy.stats` namespace; "
                   "the `scipy.stats.morestats` namespace is deprecated and "
                   "will be removed in SciPy 2.0.0.")
    else:
        message = (f"`scipy.stats.morestats.{name}` is deprecated along with "
                   "the `scipy.stats.morestats` namespace. "
                   f"`scipy.stats.morestats.{name}` will be removed in SciPy 1.13.0, and "
                   "the `scipy.stats.morestats` namespace will be removed in SciPy 2.0.0.")

    warnings.warn(message, category=DeprecationWarning, stacklevel=2)

    return getattr(_morestats, name)
