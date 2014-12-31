from __future__ import division, print_function, absolute_import

from functools import partial

import numpy as np
from . import blas
from scipy.lib._version import NumpyVersion

from .decomp import _asarray_validated
from .decomp_svd import svdvals

__all__ = ['elementwise_norm', 'frobenius_norm', 'nuclear_norm',
           'spectral_norm', 'schatten_norm', 'induced_norm', 'ky_fan_norm']

_npv = NumpyVersion(np.__version__)

_nrm2_prefix = {'f': 's', 'F': 'sc', 'D': 'dz'}

#TODO
# strategy: do not special-case the euclidean norm here.
# if > 1.8.0.dev is available then allow axes
# if > 1.10.0.dev is available then allow keepdims
# numpy multi-svd may be required
# Except for elementwise norm, all norms here require exactly two axes.
# The elementwise norm is flexible to support any number of axes,
# and by default it uses all available axes.


def _parameterized_vector_norm(p, A, axis=None):
    """
    Parameterized vector norm.

    Parameters
    ----------
    p : float
        The value of p defining the p-norm of the vector.
    A : ndarray
        The ndarray to be reduced.
    axis : int
        The axis to which the reduction should be applied.

    Returns
    -------
    n : float or ndarray
        Norm of the ndarray.

    Notes
    -----
    This is used functionally in association with the _mult_svd_norm
    to compute Schatten norms with p != 2.
    The Schatten norm with p == 2 is the Frobenius norm,
    in which case the svd is not required.

    """
    if p == 1:
        return np.absolute(A).sum(axis=axis)
    elif p == 2:
        return np.sqrt((A.conj() * A).real.sum(axis=axis)
    elif p == np.inf:
        return np.absolute(A).max(axis=axis)
    else:
        return np.power(np.power(np.absolute(A), p).sum(axis=axis), 1/p)


def _asarray_atleast2d_validated(A, check_finite):
    A = _asarray_validated(A, check_finite)
    if A.ndim < 2:
        raise ValueError('Expected the ndarray to be at least 2d')
    return A


def _euclidean_vector_norm(v):
    if v.dtype.char in 'fdFD':
        func_name = _nrm2_prefix.get(a.dtype.char, 'd') + 'nrm2'
        nrm2 = getattr(blas, func_name)
        return nrm2(v)
    else:
        return np.linalg.norm(v)


def _as_int_validated(k):
    if int(k) != k:
        raise TypeError('expected an integer')
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
    """
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


def _checked_2d_axis(nd, axis):
    """
    Parameters
    ----------
    nd : integer
        ndim of the ndarray to be reduced
    axis : {length-2 tuple of ints, None}, optional
        Defines the axes of reduction.

    Returns
    -------
    axis : length-2 tuple of ints or None
        Defines the axes of reduction.

    """
    if axis is None:
        if nd > 2:
            raise ValueError('To compute matrix norms given an array '
                             'with shape larger than 2d, a pair of axes '
                             'must be provided.')
    elif not isinstance(axis, tuple) or len(axis) != 2:
        raise TypeError("'axis' must be None or a length-2 tuple of integers")
    return axis



def _unsafe_elementwise_norm(A, p, axis=None, keepdims=None):
    """
    Elementwise norm.

    Parameters
    ----------
    A : ndarray
        The ndarray is assumed to have finite elements.
    p : float
        A floating point number greater than or equal to 1.
        Infinity is allowed, and is represented by the numpy.inf object.
        This is the value p that parameterizes the vector p-norm.
    axis : {int, tuple of ints, None}, optional
        If `axis` is an integer, it specifies the axis of `A` along which to
        compute the vector norms. If `axis` is a tuple, it specifies the
        axes that hold ndarrays, and the elementwise norms of these arrays
        are computed. If `axis` is None then the elementwise norm of
        the entire array is returned.
        Requires numpy 1.8+
    keepdims : bool, optional
        If this is set to True, the axes which are normed over are left in the
        result as dimensions with size one.  With this option the result will
        broadcast correctly against the original `A`.
        Requires numpy 1.10+

    Returns
    -------
    n : float or ndarray
        Norm of the ndarray.

    """
    if p < 1:
        raise ValueError('p must be at least 1')
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


def elementwise_norm(A, p, axis=None, keepdims=None, check_finite=True):
    A = _asarray_validated(A, check_finite=check_finite)
    return _unsafe_elementwise_norm(A, p, axis=axis, keepdims=keepdims)


def frobenius_norm(A, axis=None, keepdims=None, check_finite=True):
    A = _asarray_atleast2d_validated(A, check_finite)
    axis = _checked_2d_axis(A.ndim, axis)
    return _unsafe_elementwise_norm(A, 2, axis=axis, keepdims=keepdims)


def nuclear_norm(A, axis=None, keepdims=None, check_finite=True):
    A = _asarray_atleast2d_validated(A, check_finite)
    axis = _checked_2d_axis(A.ndim, axis)
    if axis is None and A.ndim == 2:
        return svdvals(A, check_finite=False).sum()
    else:



def spectral_norm(A, axis=None, keepdims=None, check_finite=True):
    A = _asarray_atleast2d_validated(A, check_finite)
    axis = _checked_2d_axis(axis)
    return svdvals(A, check_finite=False).max()


def schatten_norm(A, p, axis=None, keepdims=None, check_finite=True):
    if p < 1:
        raise ValueError('p must be at least 1')
    if p == 1:
        return nuclear_norm(A, check_finite=check_finite)
    elif p == 2:
        return frobenius_norm(A, check_finite=check_finite)
    elif p == np.inf:
        return spectral_norm(A, check_finite=check_finite)
    else:
        s = svdvals(A, check_finite=check_finite)
        return elementwise_norm(s, p)


def induced_norm(A, p, axis=None, keepdims=None, check_finite=True):
    _check_axis_and_keepdims(axis, keepdims)
    if p < 1:
        raise ValueError('p must be at least 1')
    if p == 1:
        A = _asarray_2d_validated(A, check_finite=check_finite)
        return np.absolute(A).sum(axis=0).max()
    elif p == 2:
        return spectral_norm(A, check_finite=check_finite)
    elif p == np.inf:
        return np.absolute(A).sum(axis=1).max()
    else:
        raise NotImplementedError('The induced norm has been implemented only '
                                  'for p in {1, 2, inf} where inf means '
                                  'the numpy.inf object.')


def ky_fan_norm(A, k, axis=None, keepdims=None, check_finite=True):
    _check_axis_and_keepdims(axis, keepdims)
    k = _as_int_validated(k)
    if k < 1:
        raise ValueError('k must be at least 1')
    A = _asarray_2d_validated(A, check_finite=check_finite)
    # The singular values are already sorted from largest to smallest.
    s = svdvals(A, check_finite=False)
    if k >= s.shape[0]:
        raise ValueError('The value of k must be less than the minimum '
                         'of the width and height of the 2d input array.')
    return s[:k].sum()
