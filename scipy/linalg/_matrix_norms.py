""" Compute matrix norms.

Note that numpy.linalg.norm allows the 'axis' keyword argument in numpy 1.8+
and allows the 'keepdims' keyword argument in numpy 1.10+.

"""
from __future__ import division, print_function, absolute_import

from functools import partial

import numpy as np
from . import blas
from scipy.lib._version import NumpyVersion

from .decomp import _asarray_validated
from .decomp_svd import svdvals

__all__ = ['entrywise_norm', 'frobenius_norm', 'nuclear_norm',
           'spectral_norm', 'schatten_norm', 'induced_norm', 'ky_fan_norm']

_npv = NumpyVersion(np.__version__)


# This function is copied from numpy, with some modifications.
def _multi_svd_norm(x, row_axis, col_axis, op):
    r"""Compute functions of singular values of the 2-D matrices in `x`.

    Parameters
    ----------
    x : ndarray
    row_axis, col_axis : int
        The axes of `x` that hold the 2-D matrices.
    op : callable
        This should reduce a 1d array of singular values to one number.
    Returns
    -------
    result : float or ndarray
        If `x` is 2-D, the return values is a float.
        Otherwise, it is an array with ``x.ndim - 2`` dimensions.
    """
    if row_axis > col_axis:
        row_axis -= 1
    y = np.rollaxis(np.rollaxis(x, col_axis, x.ndim), row_axis, -1)
    # svdvals is not yet flexible enough to replace svd here.
    return np.apply_along_axis(op, -1, np.linalg.svd(y, compute_uv=0))


def _simple_vector_p_norm(p, v):
    return np.power(np.power(np.absolute(v), p).sum(), 1/p)


def _sum_of_first_k(k, v):
    return v[:k].sum()


def _asarray_atleast2d_validated(A, check_finite):
    A = _asarray_validated(A, check_finite)
    if A.ndim < 2:
        raise ValueError('Expected the ndarray to be at least 2d')
    return A


def _as_int_validated(k):
    try:
        if int(k) != k:
            raise ValueError
    except (OverflowError, ValueError):
        raise ValueError('expected an integer')
    return int(k)


def _numpy_kwarg(value, name, version):
    if value is not None:
        if _npv < version:
            raise NotImplementedError(
                    "'%s' requires numpy version >= %s" % (name, version))
    return value


def _checked_numpy_kwargs(axis, keepdims):
    kwargs = {}
    if axis is not None:
        kwargs['axis'] = _numpy_kwarg(axis, 'axis', '1.8.0.dev')
    if keepdims is not None:
        kwargs['keepdims'] = _numpy_kwarg(keepdims, 'keepdims', '1.10.0.dev')
    return kwargs


def _checked_nd_axis(axis):
    r"""
    Parameters
    ----------
    axis : {int, tuple of ints, None}, optional
        Defines the axis or axes of reduction.

    """
    if axis is not None and not isinstance(axis, tuple):
        try:
            axis = _as_int_validated(axis)
        except TypeError:
            raise TypeError("'axis' must be None, an integer, "
                            "or a tuple of integers")
    return axis


def _checked_2d_axis(A, axis):
    r"""
    Parameters
    ----------
    A : ndarray
        The array to be reduced.
    axis : {length-2 tuple of ints, None}, optional
        Defines the axes of reduction.

    Returns
    -------
    axis : length-2 tuple of ints or None
        Defines the axes of reduction.

    """
    if axis is None:
        if A.ndim > 2:
            raise ValueError('To compute matrix norms given an array '
                             'with shape larger than 2d, a pair of axes '
                             'must be provided.')
        return None
    elif isinstance(axis, tuple) and len(axis) == 2:
        if (0 <= axis[0] < A.ndim) and (0 <= axis[1] < A.ndim):
            return axis
        else:
            raise ValueError('The array has ndim %s which is incompatible '
                             'with the axis %s.' % (A.ndim, axis))
    else:
        raise ValueError("'axis' must be None or a length-2 tuple of integers")


def _restore_dims(shape, ret, axis):
    # Restore dimensions for broadcasting.
    if axis is None:
        axis = (0, 1)
    ret_shape = list(shape)
    ret_shape[axis[0]] = 1
    ret_shape[axis[1]] = 1
    return np.reshape(ret, ret_shape)


def entrywise_norm(A, p, axis=None, keepdims=None, check_finite=True):
    r"""
    Compute the entrywise p-norm of an array.

    Compute the vector p-norm of the flattened array.

    .. math::

       \lVert A \rVert_p
           = \lVert \text{vec} \left( A \right) \rVert_p
           = \left(
                \sum_{i=1}^m
                \sum_{j=1}^n
                \lvert a_{ij} \rvert^p
             \right)^{\frac{1}{p}}
       

    Parameters
    ----------
    A : array_like
        Input array.
    p : float
        A floating point number greater than or equal to 1.
        Infinity is allowed, and is represented by the numpy.inf object.
        This is the value p that parameterizes the vector p-norm.
    axis : {int, tuple of ints, None}, optional
        If `axis` is an integer, it specifies the axis of `A` along which to
        compute the vector norms. If `axis` is a tuple, it specifies the
        axes that hold ndarrays, and the entrywise norms of these arrays
        are computed. If `axis` is None then the entrywise norm of
        the entire array is returned.
        Requires numpy 1.8+
    keepdims : bool, optional
        If this is set to True, the axes which are normed over are left in the
        result as dimensions with size one.  With this option the result will
        broadcast correctly against the original `A`.
        Requires numpy 1.10+
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.
        Default: True

    Returns
    -------
    n : float or ndarray
        Norm of the array.

    """
    if p < 1:
        raise ValueError('p must be at least 1')
    A = _asarray_validated(A, check_finite=check_finite)
    axis = _checked_nd_axis(axis)
    kwargs = _checked_numpy_kwargs(axis, keepdims)
    if p == 1:
        return np.absolute(A).sum(**kwargs)
    elif p == 2:
        return np.sqrt((A.conj() * A).real.sum(**kwargs))
    elif p == np.inf:
        return np.absolute(A).max(**kwargs)
    else:
        return np.power(np.power(np.absolute(A), p).sum(**kwargs), 1/p)


def frobenius_norm(A, axis=None, keepdims=None, check_finite=True):
    r"""
    Frobenius matrix norm.

    Compute the vector 2-norm of the flattened matrix.
    This is the square root of the sum of the squares
    of the elements of the matrix.

    .. math::

       \lVert A \rVert_2
           = \lVert \text{vec} \left( A \right) \rVert_2
           = \left(
                \sum_{i=1}^m
                \sum_{j=1}^n
                \lvert a_{ij} \rvert^2
             \right)^{\frac{1}{2}}

    Parameters
    ----------
    A : array_like
        Input array.
    axis : pair of ints, optional
        If `axis` is a pair of ints, it specifies the axes that hold 2-D
        matrices, and the norms of these matrices are computed.
        Requires numpy 1.8+
    keepdims : bool, optional
        If this is set to True, the axes which are normed over are left in the
        result as dimensions with size one.  With this option the result will
        broadcast correctly against the original `A`.
        Requires numpy 1.10+
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.
        Default: True

    Returns
    -------
    n : float or ndarray
        Norm of the array.

    """
    A = _asarray_atleast2d_validated(A, check_finite)
    axis = _checked_2d_axis(A, axis)
    kwargs = _checked_numpy_kwargs(axis, keepdims)
    original_shape = A.shape
    A = np.sqrt((A.conj() * A).real.sum(**kwargs))
    if keepdims:
        A = _restore_dims(original_shape, A, axis)
    return A


def _unitarily_invariant_norm(A, op, axis, keepdims):
    # This helper function does very little conversion or validation.
    original_shape = A.shape
    if axis is None:
        A = op(svdvals(A, check_finite=False))
    else:
        row_axis, col_axis = axis
        A = _multi_svd_norm(A, row_axis, col_axis, op)
    if keepdims:
        A = _restore_dims(original_shape, A, axis)
    return A


def nuclear_norm(A, axis=None, keepdims=None, check_finite=True):
    r"""
    Compute the nuclear norm of a matrix.

    This is the sum of the singular values of the matrix.

    .. math::

       \lVert A \rVert_{*}
           = \text{trace} \left( \sqrt{ A^{*} A } \right)
           = \sum_{i=1}^{\min \left\{ m, n \right\}} \sigma_i

    Parameters
    ----------
    A : array_like
        Input array.
    axis : pair of ints, optional
        If `axis` is a pair of ints, it specifies the axes that hold 2-D
        matrices, and the norms of these matrices are computed.
        Requires numpy 1.8+
    keepdims : bool, optional
        If this is set to True, the axes which are normed over are left in the
        result as dimensions with size one.  With this option the result will
        broadcast correctly against the original `A`.
        Requires numpy 1.10+
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.
        Default: True

    Returns
    -------
    n : float or ndarray
        Norm of the array.

    """
    A = _asarray_atleast2d_validated(A, check_finite)
    axis = _checked_2d_axis(A, axis)
    _checked_numpy_kwargs(axis, keepdims)
    return _unitarily_invariant_norm(A, np.sum, axis, keepdims)


def spectral_norm(A, axis=None, keepdims=None, check_finite=True):
    r"""
    Compute the spectral norm of a matrix.

    This is the maximum singular value of the matrix.

    .. math::

       \lVert A \rVert_{2}
           = \sqrt{ \lambda_{\max} \left( A^{*} A \right) }
           = \sigma_{\max} \left( A \right)

    Parameters
    ----------
    A : array_like
        Input array.
    axis : pair of ints, optional
        If `axis` is a pair of ints, it specifies the axes that hold 2-D
        matrices, and the norms of these matrices are computed.
        Requires numpy 1.8+
    keepdims : bool, optional
        If this is set to True, the axes which are normed over are left in the
        result as dimensions with size one.  With this option the result will
        broadcast correctly against the original `A`.
        Requires numpy 1.10+
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.
        Default: True

    Returns
    -------
    n : float or ndarray
        Norm of the array.

    """
    A = _asarray_atleast2d_validated(A, check_finite)
    axis = _checked_2d_axis(A, axis)
    _checked_numpy_kwargs(axis, keepdims)
    return _unitarily_invariant_norm(A, np.max, axis, keepdims)


def schatten_norm(A, p, axis=None, keepdims=None, check_finite=True):
    r"""
    Compute the Schatten p-norm of a matrix.

    This is the p-norm of the vector of singular values of the matrix.

    .. math::

       \lVert A \rVert_p
           = \left(
                \sum_{i=1}^{\min \left\{ m, n \right\}}
                \sigma_i^p
             \right)^{\frac{1}{p}}

    Parameters
    ----------
    A : array_like
        Input array.
    p : float
        The Schatten norm corresponds to a vector p-norm
        of the singular values of the matrix.
    axis : pair of ints, optional
        If `axis` is a pair of ints, it specifies the axes that hold 2-D
        matrices, and the norms of these matrices are computed.
        Requires numpy 1.8+
    keepdims : bool, optional
        If this is set to True, the axes which are normed over are left in the
        result as dimensions with size one.  With this option the result will
        broadcast correctly against the original `A`.
        Requires numpy 1.10+
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.
        Default: True

    Returns
    -------
    n : float or ndarray
        Norm of the array.

    """
    if p < 1:
        raise ValueError('p must be at least 1')
    if p == 1:
        return nuclear_norm(A, axis, keepdims, check_finite)
    elif p == 2:
        return frobenius_norm(A, axis, keepdims, check_finite)
    elif p == np.inf:
        return spectral_norm(A, axis, keepdims, check_finite)
    else:
        op = partial(_simple_vector_p_norm, p)
        A = _asarray_atleast2d_validated(A, check_finite)
        axis = _checked_2d_axis(A, axis)
        _checked_numpy_kwargs(axis, keepdims)
        return _unitarily_invariant_norm(A, op, axis, keepdims)


def induced_norm(A, p, axis=None, keepdims=None, check_finite=True):
    r"""
    Compute the induced p-norm of a matrix.

    .. math::

        \lVert A \rVert_p = \sup_{x \neq 0} \dfrac
            {\lVert A x \rVert_p}
            {\lVert x \rVert_p}

    Parameters
    ----------
    A : array_like
        Input array.
    p : float
        The induced norm is computed with respect to a p-norm on vectors.
    axis : pair of ints, optional
        If `axis` is a pair of ints, it specifies the axes that hold 2-D
        matrices, and the norms of these matrices are computed.
        Requires numpy 1.8+
    keepdims : bool, optional
        If this is set to True, the axes which are normed over are left in the
        result as dimensions with size one.  With this option the result will
        broadcast correctly against the original `A`.
        Requires numpy 1.10+
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.
        Default: True

    Returns
    -------
    n : float or ndarray
        Norm of the array.

    """
    if p < 1:
        raise ValueError('p must be at least 1')
    elif p == 2:
        return spectral_norm(A, axis, keepdims, check_finite)
    elif p in (1, np.inf):
        A = _asarray_atleast2d_validated(A, check_finite)
        axis = _checked_2d_axis(A, axis)
        _checked_numpy_kwargs(axis, keepdims)
        original_shape = A.shape
        if axis is None:
            axis = (0, 1)
        row_axis, col_axis = axis
        if p == 1:
            if col_axis > row_axis:
                col_axis -= 1
            A = np.absolute(A).sum(axis=row_axis).max(axis=col_axis)
        elif p == np.inf:
            if row_axis > col_axis:
                row_axis -= 1
            A = np.absolute(A).sum(axis=col_axis).max(axis=row_axis)
        if keepdims:
            A = _restore_dims(original_shape, A, axis)
        return A
    else:
        raise ValueError('The induced norm has been implemented only '
                         'for p in {1, 2, inf} where inf means '
                         'the numpy.inf object.')


def ky_fan_norm(A, k, axis=None, keepdims=None, check_finite=True):
    r"""
    Compute the Ky-Fan k-norm of a matrix.

    This is the sum of the k largest singular values of the matrix.
    So given singular values 
    :math:`\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_n \geq 0`,
    the Ky-Fan k-norm is defined as follows.

    .. math::

       \lVert A \rVert_{k} = \sum_{i=1}^k \sigma_i

    Parameters
    ----------
    A : array_like
        Input array.
    k : int
        The number of leading singular values to sum.
    axis : pair of ints, optional
        If `axis` is a pair of ints, it specifies the axes that hold 2-D
        matrices, and the norms of these matrices are computed.
        Requires numpy 1.8+
    keepdims : bool, optional
        If this is set to True, the axes which are normed over are left in the
        result as dimensions with size one.  With this option the result will
        broadcast correctly against the original `A`.
        Requires numpy 1.10+
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.
        Default: True

    Returns
    -------
    n : float or ndarray
        Norm of the array.

    """
    k = _as_int_validated(k)
    if k < 1:
        raise ValueError('k must be at least 1')
    elif k == 1:
        return spectral_norm(A, axis, keepdims, check_finite)
    else:
        A = _asarray_atleast2d_validated(A, check_finite)
        axis = _checked_2d_axis(A, axis)
        _checked_numpy_kwargs(axis, keepdims)
        if axis is None:
            spectrum_length = min(A.shape)
        else:
            spectrum_length = min(A.shape[axis[0]], A.shape[axis[1]])
        if k == spectrum_length:
            return nuclear_norm(A, axis, keepdims, check_finite=False)
        elif k < spectrum_length:
            # The singular values will be provided in decreasing order.
            op = partial(_sum_of_first_k, k)
            return _unitarily_invariant_norm(A, op, axis, keepdims)
        else:
            raise ValueError('The integer k must be no greater '
                             'than the minimum of the width and height '
                             'of the input array.')
