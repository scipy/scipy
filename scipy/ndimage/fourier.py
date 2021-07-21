from . import _fourier


__all__ = ['fourier_gaussian', 'fourier_uniform',  # noqa: F822
           'fourier_ellipsoid', 'fourier_shift']


def __dir__():
    return __all__


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(
            f"scipy.ndimage.fourier has no attribute {name}.")

    return getattr(_fourier, name)
