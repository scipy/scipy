from __future__ import division, print_function, absolute_import

import warnings

import numpy as np
from scipy._lib.six import callable
from collections import namedtuple

__all__ = ['binned_statistic',
           'binned_statistic_2d',
           'binned_statistic_dd']


BinnedStatisticResult = namedtuple('BinnedStatisticResult',
                                   ('statistic', 'bin_edges', 'binnumbers'))


def binned_statistic(x, values, statistic='mean',
                     bins=10, range=None):
    """
    Compute a binned statistic for one or more sets of data.

    This is a generalization of a histogram function.  A histogram divides
    the space into bins, and returns the count of the number of points in
    each bin.  This function allows the computation of the sum, mean, median,
    or other statistic of the values (or set of values) within each bin.

    Parameters
    ----------
    x : (N,) array_like
        A sequence of values to be binned.
    values : (N,) array_like or list of (N,) array_like
        The data on which the statistic will be computed.  This must be
        the same shape as `x`, or a list of sequences - each the same shape as
        `x`.  If `values` is a list of sequences, the statistic will be computed
        on each independently.
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
        If `bins` is an int, it defines the number of equal-width bins in the
        given range (10 by default).  If `bins` is a sequence, it defines the
        bin edges, including the rightmost edge, allowing for non-uniform bin
        widths.  Values in `x` that are smaller than lowest bin edge are
        assigned to bin number 0, values beyond the highest bin are assigned to
        ``bins[-1]``.  If the bin edges are specified, the number of bins will
        be, (nx = len(bins)-1).
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
    binnumbers: 1-D ndarray of ints
        Indices of the bins (corresponding to `bin_edges`) in which each value
        of `x` belongs.  Same length as `values`.

    See Also
    --------
    numpy.digitize, numpy.histogram, binned_statistic_2d, binned_statistic_dd

    Notes
    -----
    All but the last (righthand-most) bin is half-open.  In other words, if
    `bins` is ``[1, 2, 3, 4]``, then the first bin is ``[1, 2)`` (including 1,
    but excluding 2) and the second ``[2, 3)``.  The last bin, however, is
    ``[3, 4]``, which *includes* 4.

    .. versionadded:: 0.11.0

    Examples
    --------
    >>> from scipy import stats
    >>> import matplotlib.pyplot as plt

    First a basic example:

    >>> stats.binned_statistic([1, 2, 1, 2, 4], np.arange(5), statistic='mean',
    ...                        bins=3)
    (array([ 1.,  2.,  4.]), array([ 1.,  2.,  3.,  4.]), array([1, 2, 1, 2, 3]))

    As a second example, we now generate some random data of sailing boat speed
    as a function of wind speed, and then determine how fast our boat is for
    certain wind speeds:

    >>> windspeed = 8 * np.random.rand(500)
    >>> boatspeed = .3 * windspeed**.5 + .2 * np.random.rand(500)
    >>> bin_means, bin_edges, binnumbers= stats.binned_statistic(windspeed,
    ...                 boatspeed, statistic='median', bins=[1,2,3,4,5,6,7])
    >>> plt.figure()
    >>> plt.plot(windspeed, boatspeed, 'b.', label='raw data')
    >>> plt.hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='g', lw=5,
    ...            label='binned statistic of data')
    >>> plt.legend()

    Now we can use ``binnumbers`` to select all datapoints with a windspeed
    below 1:

    >>> low_boatspeed = boatspeed[binnumbers == 0]

    As a final example, we will use ``bin_edges`` and ``binnumbers`` to make a
    plot of a distribution that shows the mean and distribution around that
    mean per bin, on top of a regular histogram and the probability
    distribution function:

    >>> x = np.linspace(0, 5, num=500)
    >>> x_pdf = stats.maxwell.pdf(x)
    >>> samples = stats.maxwell.rvs(size=10000)

    >>> bin_means, bin_edges, binnumbers = stats.binned_statistic(x, x_pdf,
    ...         statistic='mean', bins=25)
    >>> bin_width = (bin_edges[1] - bin_edges[0])
    >>> bin_centers = bin_edges[1:] - bin_width/2

    >>> plt.figure()
    >>> plt.hist(samples, bins=50, normed=True, histtype='stepfilled', alpha=0.2,
    ...          label='histogram of data')
    >>> plt.plot(x, x_pdf, 'r-', label='analytical pdf')
    >>> plt.hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='g', lw=2,
    ...            label='binned statistic of data')
    >>> plt.plot((binnumbers - 0.5) * bin_width, x_pdf, 'g.', alpha=0.5)
    >>> plt.legend(fontsize=10)
    >>> plt.show()

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

    medians, edges, binnumbers = binned_statistic_dd(
        [x], values, statistic, bins, range)

    return BinnedStatisticResult(medians, edges[0], binnumbers)


BinnedStatistic2dResult = namedtuple('BinnedStatistic2dResult',
                                     ('statistic', 'x_edges', 'y_edges',
                                      'binnumbers'))


def binned_statistic_2d(x, y, values, statistic='mean',
                        bins=10, range=None, expand_binnumbers=False):
    """
    Compute a bidimensional binned statistic for one or more sets of data.

    This is a generalization of a histogram2d function.  A histogram divides
    the space into bins, and returns the count of the number of points in
    each bin.  This function allows the computation of the sum, mean, median,
    or other statistic of the values (or set of values) within each bin.

    Parameters
    ----------
    x : (N,) array_like
        A sequence of values to be binned along the first dimension.
    y : (N,) array_like
        A sequence of values to be binned along the second dimension.
    values : (N,) array_like or list of (N,) array_like
        The data on which the statistic will be computed.  This must be
        the same shape as `x`, or a list of sequences - each with the same
        shape as `x`.  If `values` is such a list, the statistic will be
        computed on each independently.
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

    bins : int or [int, int] or array_like or [array, array], optional
        The bin specification:

          * the number of bins for the two dimensions (nx = ny = bins),
          * the number of bins in each dimension (nx, ny = bins),
          * the bin edges for the two dimensions (x_edges = y_edges = bins),
          * the bin edges in each dimension (x_edges, y_edges = bins).

        If the bin edges are specified, the number of bins will be,
        (nx = len(x_edges)-1, ny = len(y_edges)-1).

    range : (2,2) array_like, optional
        The leftmost and rightmost edges of the bins along each dimension
        (if not specified explicitly in the `bins` parameters):
        [[xmin, xmax], [ymin, ymax]]. All values outside of this range will be
        considered outliers and not tallied in the histogram.
    expand_binnumbers : bool, optional
        .. versionadded:: 0.17.0
        'False' (default): the returned `binnumbers` is a shape (N,) array of
        linearized bin indices.
        'True': the returned `binnumbers` is 'unraveled' into a shape (2,N)
        ndarray, where each row gives the bin numbers in the corresponding
        dimension.
        See the `binnumbers` returned value.

    Returns
    -------
    statistic : (nx, ny) ndarray
        The values of the selected statistic in each two-dimensional bin.
    x_edges : (nx + 1) ndarray
        The bin edges along the first dimension.
    y_edges : (ny + 1) ndarray
        The bin edges along the second dimension.
    binnumbers : (N,) array of ints or (D,N) ndarray of ints
        This assigns to each element of `sample` an integer that represents the
        bin in which this observation falls.  The representation depends on the
        `expand_binnumbers` argument.
        If 'False' (default): The returned `binnumbers` is a shape (N,) array
        of linearized indices mapping each element of `sample` to its
        corresponding bin (using row-major ordering).
        If 'True': The returned `binnumbers` is a shape (2,N) ndarray where
        each row indicates where the elements of `sample` should be inserted,
        into the `x_edges` and `y_edges` arrays respectively, so as to keep
        the arrays sorted.


    See Also
    --------
    numpy.digitize, numpy.histogram2d, binned_statistic, binned_statistic_dd

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

    medians, edges, binnumbers = binned_statistic_dd(
        [x, y], values, statistic, bins, range,
        expand_binnumbers=expand_binnumbers)

    return BinnedStatistic2dResult(medians, edges[0], edges[1], binnumbers)


BinnedStatisticddResult = namedtuple('BinnedStatisticddResult',
                                     ('statistic', 'bin_edges',
                                      'binnumbers'))


def binned_statistic_dd(sample, values, statistic='mean',
                        bins=10, range=None, expand_binnumbers=False):
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
        The bin specification must be in one of the following forms:

          * A sequence of arrays describing the bin edges along each dimension.
          * The number of bins for each dimension (nx, ny, ... = bins).
          * The number of bins for all dimensions (nx = ny = ... = bins).

    range : sequence, optional
        A sequence of lower and upper bin edges to be used if the edges are
        not given explicitely in `bins`. Defaults to the minimum and maximum
        values along each dimension.
    expand_binnumbers : bool, optional
        .. versionadded:: 0.17.0
        'False' (default): the returned `binnumbers` is a shape (N,) array of
        linearized bin indices.
        'True': the returned `binnumbers` is 'unraveled' into a shape (D,N)
        ndarray, where each row gives the bin numbers in the corresponding
        dimension.
        See the `binnumbers` returned value.

    Returns
    -------
    statistic : ndarray, shape(nx1, nx2, nx3,...)
        The values of the selected statistic in each two-dimensional bin.
    bin_edges : list of ndarrays
        A list of D arrays describing the (nxi + 1) bin edges for each
        dimension.
    binnumbers : (N,) array of ints or (D,N) ndarray of ints
        This assigns to each element of `sample` an integer that represents the
        bin in which this observation falls.  The representation depends on the
        `expand_binnumbers` argument.
        If 'False' (default): The returned `binnumbers` is a shape (N,) array
        of linearized indices mapping each element of `sample` to its
        corresponding bin (using row-major ordering).
        If 'True': The returned `binnumbers` is a shape (D,N) ndarray where
        each row indicates where the elements of `sample` should be inserted
        into the corresponding `bin_edges` array so as to keep the array
        sorted.

    See Also
    --------
    numpy.digitize, numpy.histogramdd, binned_statistic, binned_statistic_2d

    Notes
    -----

    .. versionadded:: 0.11.0

    """
    known_stats = ['mean', 'median', 'count', 'sum', 'std']
    if not callable(statistic) and statistic not in known_stats:
        raise ValueError('invalid statistic %r' % (statistic,))

    # `Ndim` is the number of dimensions (e.g. `2` for `binned_statistic_2d`)
    # `Dlen` is the length of elements along each dimension.
    # This code is based on np.histogramdd
    try:
        # `sample` is an ND-array.
        Dlen, Ndim = sample.shape
    except (AttributeError, ValueError):
        # `sample` is a sequence of 1D arrays.
        sample = np.atleast_2d(sample).T
        Dlen, Ndim = sample.shape

    try:
        # `values` is an ND-array.
        Vdim, Vlen = values.shape
    except (AttributeError, ValueError):
        # Sample is a sequence of 1D arrays.
        values = np.atleast_2d(values)
        Vdim, Vlen = values.shape

    # Make sure `values` match `sample`
    if(Vlen != Dlen):
        raise AttributeError('The number of `values` elements must match the'
                             'length of each `sample` dimension.')

    nbin = np.empty(Ndim, int)    # Number of bins in each dimension
    edges = Ndim * [None]         # Bin edges for each dim (will be 2D array)
    dedges = Ndim * [None]        # Spacing between edges (will be 2D array)

    try:
        M = len(bins)
        if M != Ndim:
            raise AttributeError('The dimension of bins must be equal '
                                 'to the dimension of the sample x.')
    except TypeError:
        bins = Ndim * [bins]

    # Select range for each dimension
    # Used only if number of bins is given.
    if range is None:
        smin = np.atleast_1d(np.array(sample.min(axis=0), float))
        smax = np.atleast_1d(np.array(sample.max(axis=0), float))
    else:
        smin = np.zeros(Ndim)
        smax = np.zeros(Ndim)
        for i in np.arange(Ndim):
            smin[i], smax[i] = range[i]

    # Make sure the bins have a finite width.
    for i in np.arange(len(smin)):
        if smin[i] == smax[i]:
            smin[i] = smin[i] - .5
            smax[i] = smax[i] + .5

    # Create edge arrays
    for i in np.arange(Ndim):
        if np.isscalar(bins[i]):
            nbin[i] = bins[i] + 2  # +2 for outlier bins
            edges[i] = np.linspace(smin[i], smax[i], nbin[i] - 1)
        else:
            edges[i] = np.asarray(bins[i], float)
            nbin[i] = len(edges[i]) + 1  # +1 for outlier bins
        dedges[i] = np.diff(edges[i])

    nbin = np.asarray(nbin)

    # Compute the bin number each sample falls into, in each dimension
    sampBin = {}
    for i in np.arange(Ndim):
        sampBin[i] = np.digitize(sample[:, i], edges[i])

    # Using `digitize`, values that fall on an edge are put in the right bin.
    # For the rightmost bin, we want values equal to the right
    # edge to be counted in the last bin, and not as an outlier.
    for i in np.arange(Ndim):
        # Find the rounding precision
        decimal = int(-np.log10(dedges[i].min())) + 6
        # Find which points are on the rightmost edge.
        on_edge = np.where(np.around(sample[:, i], decimal) ==
                           np.around(edges[i][-1], decimal))[0]
        # Shift these points one bin to the left.
        sampBin[i][on_edge] -= 1

    # Compute the sample indices in the flattened statistic matrix.
    ni = nbin.argsort()
    # `binnumbers` is which bin (in linearized `Ndim` space) each sample goes
    binnumbers = np.zeros(Dlen, int)
    for i in np.arange(0, Ndim - 1):
        binnumbers += sampBin[ni[i]] * nbin[ni[i + 1:]].prod()
    binnumbers += sampBin[ni[-1]]

    result = np.empty([Vdim, nbin.prod()], float)

    if statistic == 'mean':
        result.fill(np.nan)
        flatcount = np.bincount(binnumbers, None)
        a = flatcount.nonzero()
        for vv in xrange(Vdim):
            flatsum = np.bincount(binnumbers, values[vv])
            result[vv, a] = flatsum[a] / flatcount[a]
    elif statistic == 'std':
        result.fill(0)
        flatcount = np.bincount(binnumbers, None)
        a = flatcount.nonzero()
        for vv in xrange(Vdim):
            flatsum = np.bincount(binnumbers, values[vv])
            flatsum2 = np.bincount(binnumbers, values[vv] ** 2)
            result[vv, a] = np.sqrt(flatsum2[a] / flatcount[a] -
                                    (flatsum[a] / flatcount[a]) ** 2)
    elif statistic == 'count':
        result.fill(0)
        flatcount = np.bincount(binnumbers, None)
        a = np.arange(len(flatcount))
        result[:, a] = flatcount[np.newaxis, :]
    elif statistic == 'sum':
        result.fill(0)
        for vv in xrange(Vdim):
            flatsum = np.bincount(binnumbers, values[vv])
            a = np.arange(len(flatsum))
            result[vv, a] = flatsum
    elif statistic == 'median':
        result.fill(np.nan)
        for i in np.unique(binnumbers):
            for vv in xrange(Vdim):
                result[vv, i] = np.median(values[vv, binnumbers == i])
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
        for i in np.unique(binnumbers):
            for vv in xrange(Vdim):
                result[vv, i] = statistic(values[vv, binnumbers == i])

    # Shape into a proper matrix
    result = result.reshape(np.append(Vdim, np.sort(nbin)))

    for i in np.arange(nbin.size):
        j = ni.argsort()[i]
        # Accomodate the extra `Vdim` dimension-zero with `+1`
        result = result.swapaxes(i+1, j+1)
        ni[i], ni[j] = ni[j], ni[i]

    # Remove outliers (indices 0 and -1 for each dimension).
    core = [slice(None)] + Ndim * [slice(1, -1)]
    result = result[core]

    # Unravel binnumbers into an ndarray, each row the bins for each dimension
    if(expand_binnumbers and Ndim > 1):
        binnumbers = np.asarray(np.unravel_index(binnumbers, nbin))

    if np.any(result.shape[1:] != nbin - 2):
        raise RuntimeError('Internal Shape Error')

    # Remove excess dimension-zero if ``Vdim == 1``
    result = result.squeeze()

    return BinnedStatisticddResult(result, edges, binnumbers)
