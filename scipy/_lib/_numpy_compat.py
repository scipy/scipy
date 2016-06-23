"""Functions copypasted from newer versions of numpy.

"""
from __future__ import division, print_function, absolute_import

import warnings
import sys

import numpy as np
from numpy.testing.nosetester import import_nose

from scipy._lib._version import NumpyVersion

if NumpyVersion(np.__version__) > '1.7.0.dev':
    _assert_warns = np.testing.assert_warns
else:
    def _assert_warns(warning_class, func, *args, **kw):
        r"""
        Fail unless the given callable throws the specified warning.

        This definition is copypasted from numpy 1.9.0.dev.
        The version in earlier numpy returns None.

        Parameters
        ----------
        warning_class : class
            The class defining the warning that `func` is expected to throw.
        func : callable
            The callable to test.
        *args : Arguments
            Arguments passed to `func`.
        **kwargs : Kwargs
            Keyword arguments passed to `func`.

        Returns
        -------
        The value returned by `func`.

        """
        with warnings.catch_warnings(record=True) as l:
            warnings.simplefilter('always')
            result = func(*args, **kw)
            if not len(l) > 0:
                raise AssertionError("No warning raised when calling %s"
                        % func.__name__)
            if not l[0].category is warning_class:
                raise AssertionError("First warning for %s is not a "
                        "%s( is %s)" % (func.__name__, warning_class, l[0]))
        return result


def assert_raises_regex(exception_class, expected_regexp,
                        callable_obj=None, *args, **kwargs):
    """
    Fail unless an exception of class exception_class and with message that
    matches expected_regexp is thrown by callable when invoked with arguments
    args and keyword arguments kwargs.
    Name of this function adheres to Python 3.2+ reference, but should work in
    all versions down to 2.6.
    Notes
    -----
    .. versionadded:: 1.8.0
    """
    __tracebackhide__ = True  # Hide traceback for py.test
    nose = import_nose()

    if sys.version_info.major >= 3:
        funcname = nose.tools.assert_raises_regex
    else:
        # Only present in Python 2.7, missing from unittest in 2.6
            funcname = nose.tools.assert_raises_regexp

    return funcname(exception_class, expected_regexp, callable_obj,
                    *args, **kwargs)
