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
    'yeojohnson_normplot'
]


def __dir__():
    return __all__


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(
            "scipy.stats.morestats is deprecated and has no attribute "
            f"{name}. Try looking in scipy.stats instead.")

    warnings.warn(f"Please use `{name}` from the `scipy.stats` namespace, "
                  "the `scipy.stats.morestats` namespace is deprecated.",
                  category=DeprecationWarning, stacklevel=2)

    return getattr(_morestats, name)
