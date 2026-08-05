# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.special` namespace for importing the functions
# included below.

from scipy._lib.deprecation import _sub_module_deprecation

from . import _orthogonal


__all__ = []
__all__ += _orthogonal.__all__
__all__ += ['airy']  # noqa: F822


def __dir__():
    return __all__


def __getattr__(name):
    return _sub_module_deprecation(sub_package="special", module="orthogonal",
                                   private_modules=["_orthogonal"], all=__all__,
                                   attribute=name)
