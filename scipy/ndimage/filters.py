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

    return getattr(_filters, name)
