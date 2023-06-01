# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.stats` namespace for importing the functions
# included below.

import warnings
from . import _morestats


__all__ = [  # noqa: F822
    'mvsdist',
    'bayes_mvs', 'kstat', 'kstatvar', 'probplot', 'ppcc_max', 'ppcc_plot',
    'boxcox_llf', 'boxcox', 'boxcox_normmax', 'boxcox_normplot',
    'shapiro', 'anderson', 'ansari', 'bartlett', 'levene', 'binom_test',
    'fligner', 'mood', 'wilcoxon', 'median_test',
    'circmean', 'circvar', 'circstd', 'anderson_ksamp',
    'yeojohnson_llf', 'yeojohnson', 'yeojohnson_normmax',
    'yeojohnson_normplot', 'find_repeats',
    'chi2_contingency', 'distributions'
]

_deprecated = {'annotations', 'namedtuple', 'isscalar', 'log', 'around',
               'unique', 'arange', 'sort', 'amin', 'amax', 'atleast_1d',
               'array', 'compress', 'exp', 'ravel', 'count_nonzero', 'arctan2',
               'hypot', 'optimize', 'Mean', 'Variance', 'Std_dev',
               'ShapiroResult', 'AndersonResult', 'Anderson_ksampResult',
               'AnsariResult', 'BartlettResult', 'LeveneResult',
               'FlignerResult', 'WilcoxonResult', 'rv_generic'}

def __dir__():
    return __all__


def __getattr__(name):
    if name not in __all__ and name not in _deprecated:
        raise AttributeError(
            "`scipy.stats.morestats` is deprecated and has no attribute "
            f"`{name}`.")

    if name in __all__:
        message = (f"Please import `{name}` from the `scipy.stats` namespace; "
                   "the `scipy.stats.morestats` namespace is deprecated.")
    else:
        message = (f"`scipy.stats.morestats.{name}` is deprecated along with "
                   "the `scipy.stats.morestats` namespace and will be removed "
                   "in SciPy 1.13.")

    warnings.warn(message, category=DeprecationWarning, stacklevel=2)

    return getattr(_morestats, name)
