from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.linalg import LinAlgError
from . import blas

__all__ = ['LinAlgError', 'norm']

_nrm2_prefix = {'f': 's', 'F': 'sc', 'D': 'dz'}


def norm(a, ord=None, axis=None):
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
    axis : {int, 2-tuple of ints, None}, optional
        If `axis` is an integer, it specifies the axis of `x` along which to
        compute the vector norms.  If `axis` is a 2-tuple, it specifies the
        axes that hold 2-D matrices, and the matrix norms of these matrices
        are computed.  If `axis` is None then either a vector norm (when `x`
        is 1-D) or a matrix norm (when `x` is 2-D) is returned.

    Returns
    -------
    norm : float or ndarray
        Norm of the matrix or vector(s).

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

    Using the `axis` argument to compute vector norms:

    >>> c = np.array([[ 1, 2, 3],
    ...               [-1, 1, 4]])
    >>> norm(c, axis=0)
    array([ 1.41421356,  2.23606798,  5.        ])
    >>> norm(c, axis=1)
    array([ 3.74165739,  4.24264069])
    >>> norm(c, ord=1, axis=1)
    array([6, 6])

    Using the `axis` argument to compute matrix norms:

    >>> m = np.arange(8).reshape(2,2,2)
    >>> norm(m, axis=(1,2))
    array([  3.74165739,  11.22497216])
    >>> norm(m[0, :, :]), norm(m[1, :, :])
    (3.7416573867739413, 11.224972160321824)

    """
    # Differs from numpy only in non-finite handling and the use of blas.
    a = np.asarray_chkfinite(a)

    # Under simple conditions, use blas for fast and stable euclidean norm.
    if (axis is None and ord in (None, 2) and
            (a.ndim == 1) and (a.dtype.char in 'fdFD')):
        func_name = _nrm2_prefix.get(a.dtype.char, 'd') + 'nrm2'
        nrm2 = getattr(blas, func_name)
        return nrm2(a)
    return np.linalg.norm(a, ord=ord, axis=axis)


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
