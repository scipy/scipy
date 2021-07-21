# This file is not meant for public use and will be removed in the future
# versions of SciPy. Use the `scipy.ndimage` namespace for importing the
# functions included below.

import warnings
from . import _filters


__all__ = ['correlate1d', 'convolve1d', 'gaussian_filter1d',  # noqa: F822
           'gaussian_filter', 'prewitt', 'sobel', 'generic_laplace',
           'laplace', 'gaussian_laplace', 'generic_gradient_magnitude',
           'gaussian_gradient_magnitude', 'correlate', 'convolve',
           'uniform_filter1d', 'uniform_filter', 'minimum_filter1d',
           'maximum_filter1d', 'minimum_filter', 'maximum_filter',
           'rank_filter', 'median_filter', 'percentile_filter',
           'generic_filter1d', 'generic_filter']


def __dir__():
    return __all__


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(
            f"scipy.ndimage.filters has no attribute {name}.")

    warnings.warn(f"Please use `{name}` from the `scipy.ndimage` namespace, "
                  "the `scipy.ndimage.filters` namespace is deprecated.",
                  category=DeprecationWarning, stacklevel=2)

    return getattr(_filters, name)
