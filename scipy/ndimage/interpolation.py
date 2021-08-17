# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.ndimage` namespace for importing the functions
# included below.

import warnings
from . import _interpolation


__all__ = ['spline_filter1d', 'spline_filter',  # noqa: F822
           'geometric_transform', 'map_coordinates',
           'affine_transform', 'shift', 'zoom', 'rotate']


def __dir__():
    return __all__


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(
            "scipy.ndimage.interpolation is deprecated and has no attribute "
            f"{name}. Try looking in scipy.ndimage instead.")

    warnings.warn(f"Please use `{name}` from the `scipy.ndimage` namespace, "
                  "the `scipy.ndimage.interpolation` namespace is deprecated.",
                  category=DeprecationWarning, stacklevel=2)

    return getattr(_interpolation, name)
