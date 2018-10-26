"""
Functions for identifying crossings in signals.
"""
from __future__ import division, print_function, absolute_import

import numpy as np


__all__ = ['argcross', 'argupcross', 'argdowncross']


def _select_cross_comparator(cross_type):
    """
    Select comparators used in finding crossings.

    Parameters
    ----------
    cross_type : str
        Accepts 'up' or 'down'.

    Return
    ------
    comparator_1 : callable
        Function to use to compare two data points.
        Should take two arrays as arguments.
    comparator_2 : callable
        Function to use to compare two data points.
        Should take two arrays as arguments.

    """
    if cross_type == 'up':
        return np.less_equal, np.greater
    elif cross_type == 'down':
        return np.greater_equal, np.less
    else:
        raise ValueError('cross_type must be "up" or "down"')


def _boolcross(data, threshold, cross_type, axis=0, mode='clip'):
    """
    Calculate the crossings of `data`.

    Crossings are calculated by finding locations where
    `data - threshold` changes sign.

    Parameters
    ----------
    data : ndarray
        Array in which to find the crossings.
    threshold : float, optional
        Value to check for upcrossings. If None is passed, the mean is used as
        threshold.
    cross_type : str
        The type of crossing to find. 'up' finds up-crossings (i.e. sign change
        from minus to pluss) and 'down' finds down-corssings (i.e. sign change
        from plus to minus).
    axis : int, optional
        Axis over which to select from `data`.  Default is 0.
    mode : str, optional
        How the edges of the vector are treated.  'wrap' (wrap around) or
        'clip' (treat overflow as the same as the last (or first) element).
        Default 'clip'.  See numpy.take

    Returns
    -------
    crossing : ndarray
        Boolean array of the same shape as `data` that is True at an crossing,
        False otherwise. The position before an crossing is returned.

    """
    datalen = data.shape[axis]
    locs = np.arange(0, datalen)
    if threshold is None:
        threshold = data.mean(axis=axis, keepdims=True)

    shape = list(data.shape)
    shape[axis] = 1
    padding = np.zeros(shape, dtype=bool)

    comparator_1, comparator_2 = _select_cross_comparator(cross_type)

    arr = data.take(locs[:-1], axis=axis, mode=mode)
    cond = comparator_1(arr, threshold)
    data.take(locs[1:], axis=axis, mode=mode, out=arr)
    cond &= comparator_2(arr, threshold)
    return np.concatenate([cond, padding], axis=axis)


def argupcross(x, threshold=None, axis=0, mode='clip'):
    """
    Calculate the up-crossings of `data`.

    The index before an up-crossing, i.e. sign change
    from minus to pluss, is returned.

    Parameters
    ----------
    data : ndarray
        Array in which to find the up-crossings.
    threshold : float, optional
        Value to check for up-crossings. If None is passed, the mean is used as
        threshold.
    axis : int, optional
        Axis over which to select from `data`.  Default is 0.
    mode : str, optional
        How the edges of the vector are treated.  'wrap' (wrap around) or
        'clip' (treat overflow as the same as the last (or first) element).
        Default 'clip'.  See numpy.take

    Returns
    -------
    upcrossing : tuple of ndarrays
        Indices of the up-crossings in arrays of integers. ``upcrossing[k]``
        is the array of indices of axis `k` of `data`. Note that the
        return value is a tuple even when `data` is one-dimensional.

    See Also
    --------
    argcross, argdowncross

    """
    results = _boolcross(x, threshold, 'up', axis=axis, mode=mode)
    return np.nonzero(results)


def argdowncross(x, threshold=None, axis=0, mode='clip'):
    """
    Calculate the down-crossings of `data`.

    The index before an down-crossing is returned. i.e. sign change
    from pluss to change, is returned.


    Parameters
    ----------
    data : ndarray
        Array in which to find the down-crossings.
    threshold : float, optional
        Value to check for down-crossings. If None is passed, the mean is
        used as threshold.
    axis : int, optional
        Axis over which to select from `data`. Default is 0.
    mode : str, optional
        How the edges of the vector are treated. 'wrap' (wrap around) or
        'clip' (treat overflow as the same as the last (or first) element).
        Default 'clip'. See numpy.take

    Returns
    -------
    downcrossing : tuple of ndarrays
        Indices of the down-crossings in arrays of integers.
        ``downcrossing[k]`` is the array of indices of axis `k` of `data`.
        Note that the return value is a tuple even when `data` is
        one-dimensional.

    See Also
    --------
    argcross, argupcross

    """
    results = _boolcross(x, threshold, 'down', axis=axis, mode=mode)
    return np.nonzero(results)


def argcross(data, threshold=None, axis=0, mode='clip'):
    """
    Calculate the crossings of `data`.

    The index before an crossing is returned.

    Parameters
    ----------
    data : ndarray
        Array in which to find the crossings.
    threshold : float, optional
        Value to check for crossings. If None is passed, the mean is used as
        threshold.
    axis : int, optional
        Axis over which to select from `data`.  Default is 0.
    mode : str, optional
        How the edges of the vector are treated. 'wrap' (wrap around) or
        'clip' (treat overflow as the same as the last (or first) element).
        Default 'clip'. See numpy.take

    Returns
    -------
    crossing : tuple of ndarrays
        Indices of the crossings in arrays of integers. ``crossing[k]`` is
        the array of indices of axis `k` of `data`. Note that the return value
        is a tuple even when `data` is one-dimensional.

    See Also
    --------
    argupcross, argdowncross

    """
    results = _boolcross(data, threshold, 'up', axis=axis, mode=mode)
    results |= _boolcross(data, threshold, 'down', axis=axis, mode=mode)
    return np.nonzero(results)
