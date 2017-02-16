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


if NumpyVersion(np.__version__) >= '1.10.0':
    from numpy import broadcast_to
else:
    # Definition of `broadcast_to` from numpy 1.10.0.

    def _maybe_view_as_subclass(original_array, new_array):
        if type(original_array) is not type(new_array):
            # if input was an ndarray subclass and subclasses were OK,
            # then view the result as that subclass.
            new_array = new_array.view(type=type(original_array))
            # Since we have done something akin to a view from original_array, we
            # should let the subclass finalize (if it has it implemented, i.e., is
            # not None).
            if new_array.__array_finalize__:
                new_array.__array_finalize__(original_array)
        return new_array

    def _broadcast_to(array, shape, subok, readonly):
        shape = tuple(shape) if np.iterable(shape) else (shape,)
        array = np.array(array, copy=False, subok=subok)
        if not shape and array.shape:
            raise ValueError('cannot broadcast a non-scalar to a scalar array')
        if any(size < 0 for size in shape):
            raise ValueError('all elements of broadcast shape must be non-'
                             'negative')
        broadcast = np.nditer(
            (array,), flags=['multi_index', 'refs_ok', 'zerosize_ok'],
            op_flags=['readonly'], itershape=shape, order='C').itviews[0]
        result = _maybe_view_as_subclass(array, broadcast)
        if not readonly and array.flags.writeable:
            result.flags.writeable = True
        return result

    def broadcast_to(array, shape, subok=False):
        return _broadcast_to(array, shape, subok=subok, readonly=True)


if NumpyVersion(np.__version__) >= '1.9.0':
    from numpy import unique
else:
    # the return_counts keyword was added in 1.9.0
    def unique(ar, return_index=False, return_inverse=False, return_counts=False):
        """
        Find the unique elements of an array.

        Returns the sorted unique elements of an array. There are three optional
        outputs in addition to the unique elements: the indices of the input array
        that give the unique values, the indices of the unique array that
        reconstruct the input array, and the number of times each unique value
        comes up in the input array.

        Parameters
        ----------
        ar : array_like
            Input array. This will be flattened if it is not already 1-D.
        return_index : bool, optional
            If True, also return the indices of `ar` that result in the unique
            array.
        return_inverse : bool, optional
            If True, also return the indices of the unique array that can be used
            to reconstruct `ar`.
        return_counts : bool, optional
            If True, also return the number of times each unique value comes up
            in `ar`.

            .. versionadded:: 1.9.0

        Returns
        -------
        unique : ndarray
            The sorted unique values.
        unique_indices : ndarray, optional
            The indices of the first occurrences of the unique values in the
            (flattened) original array. Only provided if `return_index` is True.
        unique_inverse : ndarray, optional
            The indices to reconstruct the (flattened) original array from the
            unique array. Only provided if `return_inverse` is True.
        unique_counts : ndarray, optional
            The number of times each of the unique values comes up in the
            original array. Only provided if `return_counts` is True.

            .. versionadded:: 1.9.0

        Notes
        -----
        Taken over from numpy 1.12.0-dev (c8408bf9c).  Omitted examples,
        see numpy documentation for those.

        """
        ar = np.asanyarray(ar).flatten()

        optional_indices = return_index or return_inverse
        optional_returns = optional_indices or return_counts

        if ar.size == 0:
            if not optional_returns:
                ret = ar
            else:
                ret = (ar,)
                if return_index:
                    ret += (np.empty(0, np.bool),)
                if return_inverse:
                    ret += (np.empty(0, np.bool),)
                if return_counts:
                    ret += (np.empty(0, np.intp),)
            return ret

        if optional_indices:
            perm = ar.argsort(kind='mergesort' if return_index else 'quicksort')
            aux = ar[perm]
        else:
            ar.sort()
            aux = ar
        flag = np.concatenate(([True], aux[1:] != aux[:-1]))

        if not optional_returns:
            ret = aux[flag]
        else:
            ret = (aux[flag],)
            if return_index:
                ret += (perm[flag],)
            if return_inverse:
                iflag = np.cumsum(flag) - 1
                inv_idx = np.empty(ar.shape, dtype=np.intp)
                inv_idx[perm] = iflag
                ret += (inv_idx,)
            if return_counts:
                idx = np.concatenate(np.nonzero(flag) + ([ar.size],))
                ret += (np.diff(idx),)
        return ret


if NumpyVersion(np.__version__) > '1.12.0.dev':
    polyvalfromroots = np.polynomial.polynomial.polyvalfromroots
else:
    def polyvalfromroots(x, r, tensor=True):
        """
        Evaluate a polynomial specified by its roots at points x.

        This function is copypasted from numpy 1.12.0.dev.

        If `r` is of length `N`, this function returns the value

        .. math:: p(x) = \prod_{n=1}^{N} (x - r_n)

        The parameter `x` is converted to an array only if it is a tuple or a
        list, otherwise it is treated as a scalar. In either case, either `x`
        or its elements must support multiplication and addition both with
        themselves and with the elements of `r`.

        If `r` is a 1-D array, then `p(x)` will have the same shape as `x`.  If
        `r` is multidimensional, then the shape of the result depends on the
        value of `tensor`. If `tensor is ``True`` the shape will be r.shape[1:]
        + x.shape; that is, each polynomial is evaluated at every value of `x`.
        If `tensor` is ``False``, the shape will be r.shape[1:]; that is, each
        polynomial is evaluated only for the corresponding broadcast value of
        `x`. Note that scalars have shape (,).

        Parameters
        ----------
        x : array_like, compatible object
            If `x` is a list or tuple, it is converted to an ndarray, otherwise
            it is left unchanged and treated as a scalar. In either case, `x`
            or its elements must support addition and multiplication with with
            themselves and with the elements of `r`.
        r : array_like
            Array of roots. If `r` is multidimensional the first index is the
            root index, while the remaining indices enumerate multiple
            polynomials. For instance, in the two dimensional case the roots of
            each polynomial may be thought of as stored in the columns of `r`.
        tensor : boolean, optional
            If True, the shape of the roots array is extended with ones on the
            right, one for each dimension of `x`. Scalars have dimension 0 for
            this action. The result is that every column of coefficients in `r`
            is evaluated for every element of `x`. If False, `x` is broadcast
            over the columns of `r` for the evaluation.  This keyword is useful
            when `r` is multidimensional. The default value is True.

        Returns
        -------
        values : ndarray, compatible object
            The shape of the returned array is described above.

        See Also
        --------
        polyroots, polyfromroots, polyval

        Examples
        --------
        >>> from numpy.polynomial.polynomial import polyvalfromroots
        >>> polyvalfromroots(1, [1,2,3])
        0.0
        >>> a = np.arange(4).reshape(2,2)
        >>> a
        array([[0, 1],
               [2, 3]])
        >>> polyvalfromroots(a, [-1, 0, 1])
        array([[ -0.,   0.],
               [  6.,  24.]])
        >>> r = np.arange(-2, 2).reshape(2,2) # multidimensional coefficients
        >>> r # each column of r defines one polynomial
        array([[-2, -1],
               [ 0,  1]])
        >>> b = [-2, 1]
        >>> polyvalfromroots(b, r, tensor=True)
        array([[-0.,  3.],
               [ 3., 0.]])
        >>> polyvalfromroots(b, r, tensor=False)
        array([-0.,  0.])
        """
        r = np.array(r, ndmin=1, copy=0)
        if r.dtype.char in '?bBhHiIlLqQpP':
            r = r.astype(np.double)
        if isinstance(x, (tuple, list)):
            x = np.asarray(x)
        if isinstance(x, np.ndarray):
            if tensor:
                r = r.reshape(r.shape + (1,)*x.ndim)
            elif x.ndim >= r.ndim:
                raise ValueError("x.ndim must be < r.ndim when tensor == "
                                 "False")
        return np.prod(x - r, axis=0)
