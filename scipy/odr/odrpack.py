# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.odr` namespace for importing the functions
# included below.


__all__ = [  # noqa: F822
    'odr', 'OdrWarning', 'OdrError', 'OdrStop',
    'Data', 'RealData', 'Model', 'Output', 'ODR',
    'odr_error', 'odr_stop'
]


def __dir__():
    return __all__


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(
            "`scipy.odr.odrpack` is deprecated and will be removed in SciPy 1.19.0 and "
            f"has no attribute {name}.")

    import warnings
    from . import _models
    warnings.warn("`scipy.odr` is deprecated and will be removed in SciPy 1.19.0",
                  category=DeprecationWarning, stacklevel=2)

    return getattr(_models, name)
