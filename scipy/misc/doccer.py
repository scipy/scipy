# This file is not meant for public use and will be removed in SciPy v2.0.0.


import warnings
from . import _doccer


__all__ = [  # noqa: F822
    'docformat', 'inherit_docstring_from', 'indentcount_lines',
    'filldoc', 'unindent_dict', 'unindent_string'
]


def __dir__():
    return __all__


def __getattr__(name):
    if name not in __all__:
        raise AttributeError(
            "scipy.misc.doccer is deprecated and has no attribute "
            f"{name}.")

    warnings.warn("The `scipy.misc.doccer` namespace is deprecated and "
                  "will be removed in SciPy v2.0.0.",
                  category=DeprecationWarning, stacklevel=2)

    return getattr(_doccer, name)
