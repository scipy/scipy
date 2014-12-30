from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.linalg import LinAlgError
from . import blas

__all__ = ['LinAlgError', 'norm',
           'elementwise_norm', 'frobenius_norm', 'nuclear_norm',
           'spectral_norm', 'schatten_norm', 'induced_norm', 'ky_fan_norm']

_nrm2_prefix = {'f': 's', 'F': 'sc', 'D': 'dz'}


def norm(a, ord=None):
    """
    Matrix or vector norm.

    This function is able to return one of seven different matrix norms,
    or one of an infinite number of vector norms (described below), depending
    on the value of the ``ord`` parameter.

    Parameters
    ----------
    x : (M,) or (M, N) array_like
        Input array.
    ord : {non-zero int, inf, -inf, 'fro'}, optional
        Order of the norm (see table under ``Notes``). inf means numpy's
        `inf` object.

    Returns
    -------
    norm : float
        Norm of the matrix or vector.

    Notes
    -----
    For values of ``ord <= 0``, the result is, strictly speaking, not a
    mathematical 'norm', but it may still be useful for various numerical
    purposes.

    The following norms can be calculated:

    =====  ============================  ==========================
    ord    norm for matrices             norm for vectors
    =====  ============================  ==========================
    None   Frobenius norm                2-norm
    'fro'  Frobenius norm                --
    inf    max(sum(abs(x), axis=1))      max(abs(x))
    -inf   min(sum(abs(x), axis=1))      min(abs(x))
    0      --                            sum(x != 0)
    1      max(sum(abs(x), axis=0))      as below
    -1     min(sum(abs(x), axis=0))      as below
    2      2-norm (largest sing. value)  as below
    -2     smallest singular value       as below
    other  --                            sum(abs(x)**ord)**(1./ord)
    =====  ============================  ==========================

    The Frobenius norm is given by [1]_:

        :math:`||A||_F = [\\sum_{i,j} abs(a_{i,j})^2]^{1/2}`

    References
    ----------
    .. [1] G. H. Golub and C. F. Van Loan, *Matrix Computations*,
           Baltimore, MD, Johns Hopkins University Press, 1985, pg. 15

    Examples
    --------
    >>> from scipy.linalg import norm
    >>> a = np.arange(9) - 4
    >>> a
    array([-4, -3, -2, -1,  0,  1,  2,  3,  4])
    >>> b = a.reshape((3, 3))
    >>> b
    array([[-4, -3, -2],
           [-1,  0,  1],
           [ 2,  3,  4]])

    >>> norm(a)
    7.745966692414834
    >>> norm(b)
    7.745966692414834
    >>> norm(b, 'fro')
    7.745966692414834
    >>> norm(a, np.inf)
    4
    >>> norm(b, np.inf)
    9
    >>> norm(a, -np.inf)
    0
    >>> norm(b, -np.inf)
    2

    >>> norm(a, 1)
    20
    >>> norm(b, 1)
    7
    >>> norm(a, -1)
    -4.6566128774142013e-010
    >>> norm(b, -1)
    6
    >>> norm(a, 2)
    7.745966692414834
    >>> norm(b, 2)
    7.3484692283495345

    >>> norm(a, -2)
    nan
    >>> norm(b, -2)
    1.8570331885190563e-016
    >>> norm(a, 3)
    5.8480354764257312
    >>> norm(a, -3)
    nan

    """
    # Differs from numpy only in non-finite handling and the use of
    # blas
    a = np.asarray_chkfinite(a)
    if ord in (None, 2) and (a.ndim == 1) and (a.dtype.char in 'fdFD'):
        # use blas for fast and stable euclidean norm
        func_name = _nrm2_prefix.get(a.dtype.char, 'd') + 'nrm2'
        nrm2 = getattr(blas, func_name)
        return nrm2(a)
    return np.linalg.norm(a, ord=ord)


def _datacopied(arr, original):
    """
    Strict check for `arr` not sharing any data with `original`,
    under the assumption that arr = asarray(original)

    """
    if arr is original:
        return False
    if not isinstance(original, np.ndarray) and hasattr(original, '__array__'):
        return False
    return arr.base is None


from .decomp import _asarray_validated
from .decomp_svd import svdvals


def _euclidean_vector_norm(v):
    if v.dtype.char in 'fdFD':
        func_name = _nrm2_prefix.get(a.dtype.char, 'd') + 'nrm2'
        nrm2 = getattr(blas, func_name)
        return nrm2(v)
    else:
        return np.linalg.norm(v)


def _asarray_2d_validated(A, check_finite=True):
    A = _asarray_validated(A, check_finite=check_finite)
    if A.ndim != 2:
        raise ValueError('A 2d input array is required.')
    return A


def _as_int_validated(k):
    if int(k) != k:
        raise TypeError('expected an integer')
    return int(k)


def _check_axis_and_keepdims(axis, keepdims):
    if axis is not None or keepdims is not None:
        raise NotImplementedError('The axis and keepdims arguments '
                                  'are not yet supported.')


def elementwise_norm(A, p, axis=None, keepdims=None, check_finite=True):
    _check_axis_and_keepdims(axis, keepdims)
    if p < 1:
        raise ValueError('p must be at least 1')
    A = _asarray_validated(A, check_finite=check_finite)
    if p == 1:
        return np.absolute(A).sum()
    elif p == 2:
        return _euclidean_vector_norm(A.ravel())
    elif p == np.inf:
        return np.absolute(A).max()
    else:
        return np.power(A, p).sum() ** (1 / p)


def frobenius_norm(A, axis=None, keepdims=None, check_finite=True):
    _check_axis_and_keepdims(axis, keepdims)
    A = _asarray_2d_validated(A, check_finite=check_finite)
    return _euclidean_vector_norm(A.ravel())


def nuclear_norm(A, axis=None, keepdims=None, check_finite=True):
    _check_axis_and_keepdims(axis, keepdims)
    A = _asarray_2d_validated(A, check_finite=check_finite)
    return svdvals(A, check_finite=False).sum()


def spectral_norm(A, axis=None, keepdims=None, check_finite=True):
    _check_axis_and_keepdims(axis, keepdims)
    A = _asarray_2d_validated(A, check_finite=check_finite)
    return svdvals(A, check_finite=False).max()


def schatten_norm(A, p, axis=None, keepdims=None, check_finite=True):
    _check_axis_and_keepdims(axis, keepdims)
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
