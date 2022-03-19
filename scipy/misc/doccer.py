# This file is not meant for public use and will be removed in SciPy v2.0.0.


from scipy._lib import doccer


__all__ = [  # noqa: F822
    'docformat', 'inherit_docstring_from', 'indentcount_lines',
    'filldoc', 'unindent_dict', 'unindent_string', 'extend_notes_in_docstring',
    'replace_notes_in_docstring'
]


def __dir__():
    return __all__


def __getattr__(name):
    return getattr(doccer, name)
