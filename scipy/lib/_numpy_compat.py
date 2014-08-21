"""Functions copypasted from newer versions of numpy.

"""
from __future__ import division, print_function, absolute_import

import warnings

import numpy as np

from scipy.lib._version import NumpyVersion

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


if NumpyVersion(np.__version__) >= '1.6.0':
    count_nonzero = np.count_nonzero
else:
    def count_nonzero(a):
        return (a != 0).sum()
