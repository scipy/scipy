import numpy as np
from scipy.sparse import issparse

from numpy.core import (
    array, asarray, zeros, empty, empty_like, transpose, intc, single, double,
    csingle, cdouble, inexact, complexfloating, newaxis, ravel, all, Inf, dot,
    add, multiply, sqrt, maximum, fastCopyAndTranspose, sum, isfinite, size,
    finfo, errstate, geterrobj, longdouble, rollaxis, amin, amax, product, abs,
    broadcast
    )


def isComplexType(t):
    return issubclass(t, complexfloating)

def norm(x, ord=None, axis=None):
    """
    Norm of a sparse matrix

    This function is able to return one of seven different matrix norms,
    or one of an infinite number of vector norms (described below), depending
    on the value of the ``ord`` parameter.
    
    Parameters
    ----------
    x : array_like
    Input array. If `axis` is None, `x` must be 1-D or 2-D.
    ord : {non-zero int, inf, -inf, 'fro'}, optional
    Order of the norm (see table under ``Notes``). inf means numpy's
    `inf` object.
    axis : {int, None}, optional
    If `axis` is an integer, it specifies the axis of `x` along which to
    compute the vector norms. 
    
    Returns
    -------
    n : float or matrix
    
    Notes
    -----
    Some of the ord is not implemented because some associated functions like, _multi_svd_norm, is not yet available for sparse matrix. 
    
    This docstring is modified based on numpy.linalg.norm. https://github.com/numpy/numpy/blob/master/numpy/linalg/linalg.py 
    
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
    2      Not implemented  
    -2     Not implemented      
    other  Not implemented                               
    =====  ============================  
    
    The Frobenius norm is given by [1]_:
    
    :math:`||A||_F = [\\sum_{i,j} abs(a_{i,j})^2]^{1/2}`
    
    References
    ----------
    .. [1] G. H. Golub and C. F. Van Loan, *Matrix Computations*,
    Baltimore, MD, Johns Hopkins University Press, 1985, pg. 15
    
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
    
    Using the `axis` argument to compute vector norms:
    
    >>> c = np.array([[ 1, 2, 3],
    ... [-1, 1, 4]])
    >>> c = csr_matrix(c)
    >>> norm(c, axis=0)
    matrix[[ 1.41421356, 2.23606798, 5. ]]
    >>> norm(c, axis=1)
    matrix[[ 3.74165739, 4.24264069]]
    >>> norm(c, ord=1, axis=1)
    matrix[[6]
    [6]]

"""
    if not issparse(x):
        raise TypeError("input is not sparse. use numpy.linalg.norm")

    # Check the default case first and handle it immediately.
    if ord in [None, 'fro', 'f'] and axis is None:
        if isComplexType(x.dtype.type):
            sqnorm = dot(x.real, x.real) + dot(x.imag, x.imag)
        else:
            sqnorm = x.power(2).sum()
        return sqrt(sqnorm)

    # Normalize the `axis` argument to a tuple.
    nd = x.ndim
    if axis is None:
        axis = tuple(range(nd))
    
    if np.isscalar(axis):
        if ord == Inf:
            return max(abs(x).sum(axis=axis))
        elif ord == -Inf:
            return min(abs(x).sum(axis=axis))
        elif ord == 0:
            # Zero norm
            return (x != 0).sum(axis=axis)
        elif ord == 1:
            # special case for speedup
            return abs(x).sum(axis=axis)
        elif ord == -1:
            return min(abs(x).sum(axis=axis))             
        elif ord is None:            
            return sqrt(x.power(2).sum(axis=axis))        
        else:
            raise NotImplementedError
    elif len(axis) == 2:
        row_axis, col_axis = axis
        if not (-nd <= row_axis < nd and -nd <= col_axis < nd):
            raise ValueError('Invalid axis %r for an array with shape %r' %
                             (axis, x.shape))
        if row_axis % nd == col_axis % nd:
            raise ValueError('Duplicate axes given.')
        if ord == 2:
            raise NotImplementedError
            #return _multi_svd_norm(x, row_axis, col_axis, amax)
        elif ord == -2:
            raise NotImplementedError
            #return _multi_svd_norm(x, row_axis, col_axis, amin)
        elif ord == 1:
            return abs(x).sum(axis=row_axis).max(axis=col_axis)[0,0]
        elif ord == Inf:
            return abs(x).sum(axis=col_axis).max(axis=row_axis)[0,0]
        elif ord == -1:
            return abs(x).sum(axis=row_axis).min(axis=col_axis)[0,0]
        elif ord == -Inf:
            return abs(x).sum(axis=col_axis).min(axis=row_axis)[0,0]
        elif ord in [None, 'fro', 'f']:
            return sqrt(x.power(2).sum(axis=axis))
        else:
            raise ValueError("Invalid norm order for matrices.")
    else:
        raise ValueError("Improper number of dimensions to norm.")

