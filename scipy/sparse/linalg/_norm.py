"""Sparse matrix norms.

"""
import numpy as np
from scipy.sparse import issparse
from scipy.sparse.linalg.eigen.lobpcg import lobpcg

from numpy import Inf, sqrt, abs

__all__ = ['norm']


def _sparse_frobenius_norm(x):
    if np.issubdtype(x.dtype, np.complexfloating):
        sqnorm = abs(x).power(2).sum()
    else:
        sqnorm = x.power(2).sum()
    return sqrt(sqnorm)


# More suitable for larger-scale sparse matrics
def _sparse_spectral_norm(x, ord):
    # Check the shape of 'x'
    diff = x.shape[1] - x.shape[0]
    idx = diff if diff > 0 else 0

    xHx = x.H @ x
    eigvec = np.zeros((xHx.shape[0], (ord != 2) * idx + 1))
    if ord == 2:
        # Use LOBPCG to obtain the largest sing. value of 'x' efficiently
        return sqrt(lobpcg(xHx, eigvec)[0][0])
    else:
        # Use LOBPCG to obtain the smallest sing. value of 'x' efficiently
        smallestEig = abs(lobpcg(xHx, eigvec, largest=False)[0][idx])
        # For singular matrices, smallestEig equals zero in the sense of 1e-14
        if smallestEig < 1e-14:
            smallestEig = 0.
        return sqrt(smallestEig)


def norm(x, ord=None, axis=None):
    """
    Norm of a sparse matrix

    This function is able to return one of seven different matrix norms,
    depending on the value of the ``ord`` parameter.

    Parameters
    ----------
    x : a sparse matrix
        Input sparse matrix.
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
    n : float or ndarray

    Notes
    -----
    This docstring is modified based on numpy.linalg.norm.
    https://github.com/numpy/numpy/blob/master/numpy/linalg/linalg.py

    The following norms can be calculated:

    =====  ============================
    ord    norm for sparse matrices
    =====  ============================
    None   Frobenius norm
    'fro'  Frobenius norm
    inf    max(sum(abs(x), axis=1))
    -inf   min(sum(abs(x), axis=1))
    0      abs(x).sum(axis=axis)
    1      max(sum(abs(x), axis=0))
    -1     min(sum(abs(x), axis=0))
    2      Spectral norm (largest singular value)
    -2     Smallest singular value
    other  Not implemented
    =====  ============================

    The Frobenius norm is given by [1]_, pg. 71:

        :math:`||A||_F = [\\sum_{i,j} abs(a_{i,j})^2]^{1/2}`

    The Spectral norm is given by [1]_, pg. 73:

        :math:`||A||_2 = \\sqrt(\\lambda_{\\max}(A^HA))`

    References
    ----------
    .. [1] G. H. Golub and C. F. Van Loan, *Matrix Computations*,
        4th Edition, Baltimore, Johns Hopkins University Press, 2013.

    Examples
    --------
    >>> from scipy.sparse import *
    >>> import numpy as np
    >>> from scipy.sparse.linalg import norm
    >>> a = np.arange(9) - 4
    >>> a
    array([-4, -3, -2, -1, 0, 1, 2, 3, 4])
    >>> b = a.reshape((3, 3))
    >>> b
    array([[-4, -3, -2],
           [-1, 0, 1],
           [ 2, 3, 4]])

    >>> b = csr_matrix(b)
    >>> norm(b)
    7.745966692414834
    >>> norm(b, 'fro')
    7.745966692414834
    >>> norm(b, np.inf)
    9
    >>> norm(b, -np.inf)
    2
    >>> norm(b, 1)
    7
    >>> norm(b, -1)
    6
    >>> norm(b, 2)
    7.3484692283495345
    >>> norm(b, -2)
    0.

    """
    if not issparse(x):
        raise TypeError("input is not sparse. use numpy.linalg.norm")

    # Check the default case first and handle it immediately.
    if axis is None and ord in (None, 'fro', 'f'):
        return _sparse_frobenius_norm(x)

    # Some norms require functions that are not implemented for all types.
    x = x.tocsr()

    if axis is None:
        axis = (0, 1)
    elif not isinstance(axis, tuple):
        msg = "'axis' must be None, an integer or a tuple of integers"
        try:
            int_axis = int(axis)
        except TypeError as e:
            raise TypeError(msg) from e
        if axis != int_axis:
            raise TypeError(msg)
        axis = (int_axis,)

    nd = 2
    if len(axis) == 2:
        row_axis, col_axis = axis
        if not (-nd <= row_axis < nd and -nd <= col_axis < nd):
            raise ValueError('Invalid axis %r for an array with shape %r' %
                             (axis, x.shape))
        if row_axis % nd == col_axis % nd:
            raise ValueError('Duplicate axes given.')
        if ord == 2:
            return _sparse_spectral_norm(x, ord)
        elif ord == -2:
            return _sparse_spectral_norm(x, ord)
        elif ord == 1:
            return abs(x).sum(axis=row_axis).max(axis=col_axis)[0,0]
        elif ord == Inf:
            return abs(x).sum(axis=col_axis).max(axis=row_axis)[0,0]
        elif ord == -1:
            return abs(x).sum(axis=row_axis).min(axis=col_axis)[0,0]
        elif ord == -Inf:
            return abs(x).sum(axis=col_axis).min(axis=row_axis)[0,0]
        elif ord in (None, 'f', 'fro'):
            # The axis order does not matter for this norm.
            return _sparse_frobenius_norm(x)
        else:
            raise ValueError("Invalid norm order for matrices.")
    elif len(axis) == 1:
        a, = axis
        if not (-nd <= a < nd):
            raise ValueError('Invalid axis %r for an array with shape %r' %
                             (axis, x.shape))
        if ord == Inf:
            M = abs(x).max(axis=a)
        elif ord == -Inf:
            M = abs(x).min(axis=a)
        elif ord == 0:
            # Zero norm
            M = (x != 0).sum(axis=a)
        elif ord == 1:
            # special case for speedup
            M = abs(x).sum(axis=a)
        elif ord in (2, None):
            M = sqrt(abs(x).power(2).sum(axis=a))
        else:
            try:
                ord + 1
            except TypeError as e:
                raise ValueError('Invalid norm order for vectors.') from e
            M = np.power(abs(x).power(ord).sum(axis=a), 1 / ord)
        return M.A.ravel()
    else:
        raise ValueError("Improper number of dimensions to norm.")
