"""
Functions for identifying peaks in signals.
"""

import numpy as np


def argrelextrema(data, comparator,
                  axis=0, order=1, mode='clip'):
    """
    Calculate the relative extrema of `data`.

    Relative extrema are calculated by finding locations where
    comparator(data[n],data[n+1:n+order+1]) = True.

    Parameters
    -----------
    data: array-like
    comparator: function
        function to use to compare two data points.
        Should take 2 numbers as arguments
    axis: integer, optional
        axis over which to select data
    order: integer, optional
        How many points on each side to require
        a comparator(n,n+x) = True.
    mode: string, optional
        How the edges of the vector are treated.
        'wrap' (wrap around) or 'clip' (treat overflow
        as the same as the last (or first) element).
        Default 'clip'. See numpy.take

    Returns
    ----------
    extrema: array-like
        Indices of the extrema, as boolean array
        of same shape as data. True for an extrema,
        False else.

    See also
    --------
    argrelmax,argrelmin
    """

    if((int(order) != order) or (order < 1)):
        raise ValueError('Order must be an integer >= 1')

    datalen = data.shape[axis]
    locs = np.arange(0, datalen)

    results = np.ones(data.shape, dtype=bool)
    main = data.take(locs, axis=axis, mode=mode)
    for shift in xrange(1, order + 1):
        plus = data.take(locs + shift, axis=axis, mode=mode)
        minus = data.take(locs - shift, axis=axis, mode=mode)
        results &= comparator(main, plus)
        results &= comparator(main, minus)
        if(~results.any()):
            return results
            #return (np.array([]),)*2
    #arglocs = np.where(results)
    return results


def argrelmin(data, axis=0, order=1, mode='clip'):
    """
    Calculate the relative minima of `data`.

    See also
    --------
    argrelextrema,argrelmax
    """
    return argrelextrema(data, np.less, axis, order, mode)


def argrelmax(data, axis=0, order=1, mode='clip'):
    """
    Calculate the relative maxima of `data`.

    See also
    --------
    argrelextrema,argrelmin
    """
    return argrelextrema(data, np.greater, axis, order, mode)