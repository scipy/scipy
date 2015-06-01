from __future__ import division, print_function, absolute_import

import functools
import operator
import sys
import warnings
import numbers
import inspect

import numpy as np


def _aligned_zeros(shape, dtype=float, order="C", align=None):
    """Allocate a new ndarray with aligned memory.

    Primary use case for this currently is working around a f2py issue
    in Numpy 1.9.1, where dtype.alignment is such that np.zeros() does
    not necessarily create arrays aligned up to it.

    """
    dtype = np.dtype(dtype)
    if align is None:
        align = dtype.alignment
    if not hasattr(shape, '__len__'):
        shape = (shape,)
    size = functools.reduce(operator.mul, shape) * dtype.itemsize
    buf = np.empty(size + align + 1, np.uint8)
    offset = buf.__array_interface__['data'][0] % align
    if offset != 0:
        offset = align - offset
    # Note: slices producing 0-size arrays do not necessarily change
    # data pointer --- so we use and allocate size+1
    buf = buf[offset:offset+size+1][:-1]
    data = np.ndarray(shape, dtype, buf, order=order)
    data.fill(0)
    return data


class DeprecatedImport(object):
    """
    Deprecated import, with redirection + warning.

    Examples
    --------
    Suppose you previously had in some module::

        from foo import spam

    If this has to be deprecated, do::

        spam = DeprecatedImport("foo.spam", "baz")

    to redirect users to use "baz" module instead.

    """

    def __init__(self, old_module_name, new_module_name):
        self._old_name = old_module_name
        self._new_name = new_module_name
        __import__(self._new_name)
        self._mod = sys.modules[self._new_name]

    def __dir__(self):
        return dir(self._mod)

    def __getattr__(self, name):
        warnings.warn("Module %s is deprecated, use %s instead"
                      % (self._old_name, self._new_name),
                      DeprecationWarning)
        return getattr(self._mod, name)


# copy-pasted from scikit-learn utils/validation.py
def check_random_state(seed):
    """Turn seed into a np.random.RandomState instance

    If seed is None (or np.random), return the RandomState singleton used
    by np.random.
    If seed is an int, return a new RandomState instance seeded with seed.
    If seed is already a RandomState instance, return it.
    Otherwise raise ValueError.
    """
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, (numbers.Integral, np.integer)):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError('%r cannot be used to seed a numpy.random.RandomState'
                     ' instance' % seed)


def _asarray_validated(a, check_finite=True,
                       sparse_ok=False, objects_ok=False, mask_ok=False):
    """
    Helper function for scipy argument validation.

    Many scipy linear algebra functions do support arbitrary array-like
    input arguments.  Examples of commonly unsupported inputs include
    matrices containing inf/nan, sparse matrix representations, and
    matrices with complicated elements.

    Parameters
    ----------
    a : array_like
        The array-like input.
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.
        Default: True
    sparse_ok : bool, optional
        True if scipy sparse matrices are allowed.
    objects_ok : bool, optional
        True if arrays with dype('O') are allowed.
    mask_ok : bool, optional
        True if masked arrays are allowed.

    Returns
    -------
    ret : ndarray
        The converted validated array.

    """
    if not sparse_ok:
        import scipy.sparse
        if scipy.sparse.issparse(a):
            msg = ('Sparse matrices are not supported by this function. '
                   'Perhaps one of the scipy.linalg.sparse functions '
                   'would work instead.')
            raise ValueError(msg)
    if not mask_ok:
        if np.ma.isMaskedArray(a):
            raise ValueError('masked arrays are not supported')
    if check_finite:
        a = np.asarray_chkfinite(a)
    else:
        a = np.asarray(a)
    if not objects_ok:
        if a.dtype is np.dtype('O'):
            raise ValueError('object arrays are not supported')
    return a


def inspect_args(func):
    """
    Args/defaults from inspect.getargspec (py2k) or inspect.signature (py3k).

    `inspect.signature` does not exist in Python 2.x, and `inspect.getargspec`
    is deprecated in Python 3.x.  This function is a simple wrapper around both
    that returns only the parameters and defaults lists (typically all that is
    needed within Scipy).

    Parameters
    ----------
    func : callable
        The function to inspect.

    Returns
    -------
    args : sequence of str
        The argument names.  For Python 2.x, a list of strings.  For Python
        3.x, an OrderedDict.
    defaults : list
        A list of defaults for keyword arguments.

    """
    if sys.version_info[0] >= 3:
        sig = inspect.signature(func)
        args = sig.parameters
        tmp_defaults = [d.default for d in sig.parameters.values()]
        defaults = [d for d in tmp_defaults if d is not inspect._empty]
    else:
        args, _, _, defaults = inspect.getargspec(func)

    return args, defaults
