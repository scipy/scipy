from . import _interpolation


__all__ = ['spline_filter1d', 'spline_filter',  # noqa: F822
           'geometric_transform', 'map_coordinates',
           'affine_transform', 'shift', 'zoom', 'rotate']


def __dir__():
    return __all__


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(
            f"scipy.ndimage.interpolation has no attribute {name}.")

    return getattr(_interpolation, name)
