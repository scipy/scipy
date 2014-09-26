from __future__ import division, print_function, absolute_import

import warnings

import numpy as np
from scipy.lib.six import callable, string_types


def binned_statistic(x, values, statistic='mean',
                     bins=10, range=None):
    """
    Compute a binned statistic for a set of data.

    This is a generalization of a histogram function.  A histogram divides
    the space into bins, and returns the count of the number of points in
    each bin.  This function allows the computation of the sum, mean, median,
    or other statistic of the values within each bin.

    Parameters
    ----------
    x : array_like
        A sequence of values to be binned.
    values : array_like
        The values on which the statistic will be computed.  This must be
        the same shape as `x`.
    statistic : string or callable, optional
        The statistic to compute (default is 'mean').
        The following statistics are available:

          * 'mean' : compute the mean of values for points within each bin.
            Empty bins will be represented by NaN.
          * 'median' : compute the median of values for points within each
            bin. Empty bins will be represented by NaN.
          * 'count' : compute the count of points within each bin.  This is
            identical to an unweighted histogram.  `values` array is not
            referenced.
          * 'sum' : compute the sum of values for points within each bin.
            This is identical to a weighted histogram.
          * function : a user-defined function which takes a 1D array of
            values, and outputs a single numerical statistic. This function
            will be called on the values in each bin.  Empty bins will be
            represented by function([]), or NaN if this returns an error.

    bins : int or sequence of scalars, optional
        If `bins` is an int, it defines the number of equal-width
        bins in the given range (10, by default). If `bins` is a sequence,
        it defines the bin edges, including the rightmost edge, allowing
        for non-uniform bin widths.
    range : (float, float) or [(float, float)], optional
        The lower and upper range of the bins.  If not provided, range
        is simply ``(x.min(), x.max())``.  Values outside the range are
        ignored.

    Returns
    -------
    statistic : array
        The values of the selected statistic in each bin.
    bin_edges : array of dtype float
        Return the bin edges ``(length(statistic)+1)``.
    binnumber : 1-D ndarray of ints
        This assigns to each observation an integer that represents the bin
        in which this observation falls. Array has the same length as values.

    See Also
    --------
    numpy.histogram, binned_statistic_2d, binned_statistic_dd

    Notes
    -----
    All but the last (righthand-most) bin is half-open.  In other words, if
    `bins` is::

      [1, 2, 3, 4]

    then the first bin is ``[1, 2)`` (including 1, but excluding 2) and the
    second ``[2, 3)``.  The last bin, however, is ``[3, 4]``, which *includes*
    4.

    .. versionadded:: 0.11.0

    Examples
    --------
    >>> stats.binned_statistic([1, 2, 1, 2, 4], np.arange(5), statistic='mean',
    ... bins=3)
    (array([ 1.,  2.,  4.]), array([ 1.,  2.,  3.,  4.]), array([1, 2, 1, 2, 3]))

    >>> stats.binned_statistic([1, 2, 1, 2, 4], np.arange(5), statistic='mean', bins=3)
    (array([ 1.,  2.,  4.]), array([ 1.,  2.,  3.,  4.]), array([1, 2, 1, 2, 3]))

    """
    try:
        N = len(bins)
    except TypeError:
        N = 1

    if N != 1:
        bins = [np.asarray(bins, float)]

    if range is not None:
        if len(range) == 2:
            range = [range]

    medians, edges, xy = binned_statistic_dd([x], values, statistic,
                                             bins, range)

    return medians, edges[0], xy


def binned_statistic_2d(x, y, values, statistic='mean',
                        bins=10, range=None):
    """
    Compute a bidimensional binned statistic for a set of data.

    This is a generalization of a histogram2d function.  A histogram divides
    the space into bins, and returns the count of the number of points in
    each bin.  This function allows the computation of the sum, mean, median,
    or other statistic of the values within each bin.

    Parameters
    ----------
    x : (N,) array_like
        A sequence of values to be binned along the first dimension.
    y : (M,) array_like
        A sequence of values to be binned along the second dimension.
    values : (N,) array_like
        The values on which the statistic will be computed.  This must be
        the same shape as `x`.
    statistic : string or callable, optional
        The statistic to compute (default is 'mean').
        The following statistics are available:

          * 'mean' : compute the mean of values for points within each bin.
            Empty bins will be represented by NaN.
          * 'median' : compute the median of values for points within each
            bin. Empty bins will be represented by NaN.
          * 'count' : compute the count of points within each bin.  This is
            identical to an unweighted histogram.  `values` array is not
            referenced.
          * 'sum' : compute the sum of values for points within each bin.
            This is identical to a weighted histogram.
          * function : a user-defined function which takes a 1D array of
            values, and outputs a single numerical statistic. This function
            will be called on the values in each bin.  Empty bins will be
            represented by function([]), or NaN if this returns an error.

    bins : int or [int, int] or array-like or [array, array], optional
        The bin specification:

          * the number of bins for the two dimensions (nx=ny=bins),
          * the number of bins in each dimension (nx, ny = bins),
          * the bin edges for the two dimensions (x_edges = y_edges = bins),
          * the bin edges in each dimension (x_edges, y_edges = bins).

    range : (2,2) array_like, optional
        The leftmost and rightmost edges of the bins along each dimension
        (if not specified explicitly in the `bins` parameters):
        [[xmin, xmax], [ymin, ymax]]. All values outside of this range will be
        considered outliers and not tallied in the histogram.

    Returns
    -------
    statistic : (nx, ny) ndarray
        The values of the selected statistic in each two-dimensional bin
    xedges : (nx + 1) ndarray
        The bin edges along the first dimension.
    yedges : (ny + 1) ndarray
        The bin edges along the second dimension.
    binnumber : 1-D ndarray of ints
        This assigns to each observation an integer that represents the bin
        in which this observation falls. Array has the same length as `values`.

    See Also
    --------
    numpy.histogram2d, binned_statistic, binned_statistic_dd

    Notes
    -----

    .. versionadded:: 0.11.0

    """

    # This code is based on np.histogram2d
    try:
        N = len(bins)
    except TypeError:
        N = 1

    if N != 1 and N != 2:
        xedges = yedges = np.asarray(bins, float)
        bins = [xedges, yedges]

    medians, edges, xy = binned_statistic_dd([x, y], values, statistic,
                                             bins, range)

    return medians, edges[0], edges[1], xy


def binned_statistic_dd(sample, values, statistic='mean',
                        bins=10, range=None):
    """
    Compute a multidimensional binned statistic for a set of data.

    This is a generalization of a histogramdd function.  A histogram divides
    the space into bins, and returns the count of the number of points in
    each bin.  This function allows the computation of the sum, mean, median,
    or other statistic of the values within each bin.

    Parameters
    ----------
    sample : array_like
        Data to histogram passed as a sequence of D arrays of length N, or
        as an (N,D) array.
    values : array_like
        The values on which the statistic will be computed.  This must be
        the same shape as x.
    statistic : string or callable, optional
        The statistic to compute (default is 'mean').
        The following statistics are available:

          * 'mean' : compute the mean of values for points within each bin.
            Empty bins will be represented by NaN.
          * 'median' : compute the median of values for points within each
            bin. Empty bins will be represented by NaN.
          * 'count' : compute the count of points within each bin.  This is
            identical to an unweighted histogram.  `values` array is not
            referenced.
          * 'sum' : compute the sum of values for points within each bin.
            This is identical to a weighted histogram.
          * function : a user-defined function which takes a 1D array of
            values, and outputs a single numerical statistic. This function
            will be called on the values in each bin.  Empty bins will be
            represented by function([]), or NaN if this returns an error.

    bins : sequence or int, optional
        The bin specification:

          * A sequence of arrays describing the bin edges along each dimension.
          * The number of bins for each dimension (nx, ny, ... =bins)
          * The number of bins for all dimensions (nx=ny=...=bins).

    range : sequence, optional
        A sequence of lower and upper bin edges to be used if the edges are
        not given explicitely in `bins`. Defaults to the minimum and maximum
        values along each dimension.

    Returns
    -------
    statistic : ndarray, shape(nx1, nx2, nx3,...)
        The values of the selected statistic in each two-dimensional bin
    edges : list of ndarrays
        A list of D arrays describing the (nxi + 1) bin edges for each
        dimension
    binnumber : 1-D ndarray of ints
        This assigns to each observation an integer that represents the bin
        in which this observation falls. Array has the same length as values.

    See Also
    --------
    np.histogramdd, binned_statistic, binned_statistic_2d

    Notes
    -----

    .. versionadded:: 0.11.0

    """
    if isinstance(statistic, string_types):
        if statistic not in ['mean', 'median', 'count', 'sum', 'std']:
            raise ValueError('unrecognized statistic "%s"' % statistic)
    elif callable(statistic):
        pass
    else:
        raise ValueError("statistic not understood")

    # This code is based on np.histogramdd
    try:
        # Sample is an ND-array.
        N, D = sample.shape
    except (AttributeError, ValueError):
        # Sample is a sequence of 1D arrays.
        sample = np.atleast_2d(sample).T
        N, D = sample.shape

    nbin = np.empty(D, int)
    edges = D * [None]
    dedges = D * [None]

    try:
        M = len(bins)
        if M != D:
            raise AttributeError('The dimension of bins must be equal '
                                 'to the dimension of the sample x.')
    except TypeError:
        bins = D * [bins]

    # Select range for each dimension
    # Used only if number of bins is given.
    if range is None:
        smin = np.atleast_1d(np.array(sample.min(0), float))
        smax = np.atleast_1d(np.array(sample.max(0), float))
    else:
        smin = np.zeros(D)
        smax = np.zeros(D)
        for i in np.arange(D):
            smin[i], smax[i] = range[i]

    # Make sure the bins have a finite width.
    for i in np.arange(len(smin)):
        if smin[i] == smax[i]:
            smin[i] = smin[i] - .5
            smax[i] = smax[i] + .5

    # Create edge arrays
    for i in np.arange(D):
        if np.isscalar(bins[i]):
            nbin[i] = bins[i] + 2  # +2 for outlier bins
            edges[i] = np.linspace(smin[i], smax[i], nbin[i] - 1)
        else:
            edges[i] = np.asarray(bins[i], float)
            nbin[i] = len(edges[i]) + 1  # +1 for outlier bins
        dedges[i] = np.diff(edges[i])

    nbin = np.asarray(nbin)

    # Compute the bin number each sample falls into.
    Ncount = {}
    for i in np.arange(D):
        Ncount[i] = np.digitize(sample[:, i], edges[i])

    # Using digitize, values that fall on an edge are put in the right bin.
    # For the rightmost bin, we want values equal to the right
    # edge to be counted in the last bin, and not as an outlier.
    for i in np.arange(D):
        # Rounding precision
        decimal = int(-np.log10(dedges[i].min())) + 6
        # Find which points are on the rightmost edge.
        on_edge = np.where(np.around(sample[:, i], decimal)
                           == np.around(edges[i][-1], decimal))[0]
        # Shift these points one bin to the left.
        Ncount[i][on_edge] -= 1

    # Compute the sample indices in the flattened statistic matrix.
    ni = nbin.argsort()
    xy = np.zeros(N, int)
    for i in np.arange(0, D - 1):
        xy += Ncount[ni[i]] * nbin[ni[i + 1:]].prod()
    xy += Ncount[ni[-1]]

    result = np.empty(nbin.prod(), float)

    if statistic == 'mean':
        result.fill(np.nan)
        flatcount = np.bincount(xy, None)
        flatsum = np.bincount(xy, values)
        a = flatcount.nonzero()
        result[a] = flatsum[a] / flatcount[a]
    elif statistic == 'std':
        result.fill(0)
        flatcount = np.bincount(xy, None)
        flatsum = np.bincount(xy, values)
        flatsum2 = np.bincount(xy, values ** 2)
        a = flatcount.nonzero()
        result[a] = np.sqrt(flatsum2[a] / flatcount[a]
                            - (flatsum[a] / flatcount[a]) ** 2)
    elif statistic == 'count':
        result.fill(0)
        flatcount = np.bincount(xy, None)
        a = np.arange(len(flatcount))
        result[a] = flatcount
    elif statistic == 'sum':
        result.fill(0)
        flatsum = np.bincount(xy, values)
        a = np.arange(len(flatsum))
        result[a] = flatsum
    elif statistic == 'median':
        result.fill(np.nan)
        for i in np.unique(xy):
            result[i] = np.median(values[xy == i])
    elif callable(statistic):
        with warnings.catch_warnings():
            # Numpy generates a warnings for mean/std/... with empty list
            warnings.filterwarnings('ignore', category=RuntimeWarning)
            old = np.seterr(invalid='ignore')
            try:
                null = statistic([])
            except:
                null = np.nan
            np.seterr(**old)
        result.fill(null)
        for i in np.unique(xy):
            result[i] = statistic(values[xy == i])

    # Shape into a proper matrix
    result = result.reshape(np.sort(nbin))
    for i in np.arange(nbin.size):
        j = ni.argsort()[i]
        result = result.swapaxes(i, j)
        ni[i], ni[j] = ni[j], ni[i]

    # Remove outliers (indices 0 and -1 for each dimension).
    core = D * [slice(1, -1)]
    result = result[core]

    if (result.shape != nbin - 2).any():
        raise RuntimeError('Internal Shape Error')

    return result, edges, xy
