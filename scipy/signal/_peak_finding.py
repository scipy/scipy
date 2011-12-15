"""
Functions for identifying peaks in signals.
"""

import numpy as np

from scipy.signal.wavelets import cwt, ricker
from scipy.stats import scoreatpercentile


def _boolrelextrema(data, comparator,
                  axis=0, order=1, mode='clip'):
    """
    Calculate the relative extrema of `data`.

    Relative extrema are calculated by finding locations where
    comparator(data[n],data[n+1:n+order+1]) = True.

    Parameters
    ----------
    data: ndarray
    comparator: function
        function to use to compare two data points.
        Should take 2 numbers as arguments
    axis: int, optional
        axis over which to select from `data`
    order: int, optional
        How many points on each side to require
        a `comparator`(n,n+x) = True.
    mode: string, optional
        How the edges of the vector are treated.
        'wrap' (wrap around) or 'clip' (treat overflow
        as the same as the last (or first) element).
        Default 'clip'. See numpy.take

    Returns
    -------
    extrema: ndarray
        Indices of the extrema, as boolean array
        of same shape as data. True for an extrema,
        False else.

    See also
    --------
    argrelmax,argrelmin

    Examples
    --------
    >>> testdata = np.array([1,2,3,2,1])
    >>> argrelextrema(testdata, np.greater, axis=0)
    array([False, False,  True, False, False], dtype=bool)
    """

    if((int(order) != order) or (order < 1)):
        raise ValueError('Order must be an int >= 1')

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


def argrelextrema(data, comparator,
                  axis=0, order=1, mode='clip'):
    """
    Calculate the relative extrema of `data`

    Returns
    -------
    extrema: ndarray
        Indices of the extrema, as an array
        of integers (same format as argmin, argmax

    See also
    --------
    argrelmin, argrelmax

    """
    results = _boolrelextrema(data, comparator,
                              axis, order, mode)
    if ~results.any():
        return (np.array([]),) * 2
    else:
        return np.where(results)


def _identify_ridge_lines(matr, max_distances, gap_thresh):
    """
    Identify ridges in the 2D matrix. Expect that the width of
    the wavelet feature increases with increasing row number.

    Parameters
    ----------
    matr: 2-D ndarray
        Matrix in which to identify ridge lines.
    max_distances: 1-D sequence
        At each row, a ridge line is only connected
        if the relative max at row[n] is within
        `max_distances`[n] from the relative max at row[n+1].
    gap_thresh: int
        If a relative maximum is not found within `max_distances`,
        there will be a gap. A ridge line is discontinued if
        there are more than `gap_thresh` points without connecting
        a new relative maximum.

    Returns
    -------
    ridge_lines: tuple
        tuple of 2 1-D sequences. `ridge_lines`[ii][0] are the rows of the ii-th
        ridge-line, `ridge_lines`[ii][1] are the columns. Empty if none found.
        Each ridge-line will be sorted by row (increasing), but the order
        of the ridge lines is not specified

    References
    ----------
    Bioinformatics (2006) 22 (17): 2059-2065.
    doi: 10.1093/bioinformatics/btl355
    http://bioinformatics.oxfordjournals.org/content/22/17/2059.long

    Examples
    --------
    >>> data = np.random.rand(5,5)
    >>> ridge_lines = identify_ridge_lines(data, 1, 1)

    Notes:
    ------
    This function is intended to be used in conjuction with `cwt`
    as part of find_peaks_cwt.
    """

    if(len(max_distances) < matr.shape[0]):
        raise ValueError('Max_distances must have at least as many rows as matr')

    all_max_cols = _boolrelextrema(matr, np.greater, axis=1, order=1)
    #Highest row for which there are any relative maxima
    has_relmax = np.where(all_max_cols.any(axis=1))[0]
    if(len(has_relmax) == 0):
        return []
    start_row = has_relmax[-1]
    #Each ridge line is a 3-tuple:
    #rows, cols,Gap number
    ridge_lines = [[[start_row],
                   [col],
                   0] for col in np.where(all_max_cols[start_row])[0]]
    final_lines = []
    rows = np.arange(start_row - 1, -1, -1)
    cols = np.arange(0, matr.shape[1])
    for row in rows:
        this_max_cols = cols[all_max_cols[row]]

        #Increment gap number of each line,
        #set it to zero later if appropriate
        for line in ridge_lines:
            line[2] += 1

        #XXX These should always be all_max_cols[row]
        #But the order might be different. Might be an efficiency gain
        #to make sure the order is the same and avoid this iteration
        prev_ridge_cols = np.array([line[1][-1] for line in ridge_lines])
        #Look through every relative maximum found at current row
        #Attempt to connect them with existing ridge lines.
        new_lines = []
        for ind, col in enumerate(this_max_cols):
            """
            If there is a previous ridge line within
            the max_distance to connect to, do so.
            Otherwise start a new one.
            """
            line = None
            if(len(prev_ridge_cols) > 0):
                diffs = np.abs(col - prev_ridge_cols)
                closest = np.argmin(diffs)
                if diffs[closest] <= max_distances[row]:
                    line = ridge_lines[closest]
            if(line is not None):
                #Found a point close enough, extend current ridge line
                line[1].append(col)
                line[0].append(row)
                line[2] = 0
            else:
                new_line = [[row],
                            [col],
                            0]
                ridge_lines.append(new_line)

        #Remove the ridge lines with gap_number too high
        #XXX Modifying a list while iterating over it.
        #Should be safe, since we iterate backwards, but
        #still tacky.
        for ind in xrange(len(ridge_lines) - 1, -1, -1):
            line = ridge_lines[ind]
            if line[2] > gap_thresh:
                final_lines.append(line)
                del ridge_lines[ind]

    out_lines = []
    for line in (final_lines + ridge_lines):
        sortargs = np.array(np.argsort(line[0]))
        rows, cols = np.zeros_like(sortargs), np.zeros_like(sortargs)
        rows[sortargs] = line[0]
        cols[sortargs] = line[1]
        out_lines.append([rows, cols])
    return out_lines


def _filter_ridge_lines(cwt, ridge_lines, window_size=None, min_length=None,
                       min_snr=1, noise_perc=10):
    """
    Filter ridge lines according to prescribed criteria. Intended
    to be used for finding relative maxima.

    Parameters
    -------------
    cwt : 2-D ndarray
        Continuous wavelet transform from which
        the ridge_lines were defined
    ridge_lines: 1-D sequence
        Each element should contain 2 sequences, the rows and columns
        of the ridge line (respectively)
    window_size: int, optional
        Size of window to use to calculate noise floor.
        Default is `cwt`.shape[1]/20
    min_length: int, optional
        Minimum length a ridge line needs to be acceptable.
        Default is `cwt`.shape[0]/4, ie 1/4th the number of widths.
    min_snr: float, optional
        Minimum SNR ratio. Default 0. The signal is the value of
        the cwt matrix at the shortest length scale (`cwt`[0,loc]), the noise is
        the `noise_perc`th percentile of datapoints contained within
        a window of `window_size` around `cwt`[0,loc]
    noise_perc: float,optional
        When calculating the noise floor, percentile of data points
        examined below which to consider noise. Calculated using
        scipy.stats.scoreatpercentile.

    References
    ----------
    Bioinformatics (2006) 22 (17): 2059-2065. doi: 10.1093/bioinformatics/btl355
    http://bioinformatics.oxfordjournals.org/content/22/17/2059.long

    """
    num_points = cwt.shape[1]
    if min_length is None:
        min_length = np.ceil(cwt.shape[0] / 4)
    if window_size is None:
        window_size = np.ceil(num_points / 20)
    hf_window = window_size / 2

    #Filter based on SNR
    row_one = cwt[0, :]
    noises = np.zeros_like(row_one)
    for ind, val in enumerate(row_one):
        window = np.arange(max([ind - hf_window, 0]), min([ind + hf_window, num_points]))
        window = window.astype(int)
        noises[ind] = scoreatpercentile(row_one[window], per=noise_perc)

    def filt_func(line):
        if len(line[0]) < min_length:
            return False
        snr = abs(cwt[line[0][0], line[1][0]] / noises[line[1][0]])
        if snr < min_snr:
            return False
        return True

    return filter(filt_func, ridge_lines)


def find_peaks_cwt(vector, widths, wavelet=None, max_distances=None, gap_thresh=None,
               min_length=None, min_snr=1, noise_perc=10):
    """
    Attempt to find the peaks in the given 1-D array `vector`.

    The general approach is to smooth `vector` by convolving it with `wavelet(width)`
    for each width in `widths`. Relative maxima which appear at enough length scales,
    and with sufficiently high SNR, are accepted.

    Parameters
    ----------
    vector: 1-D ndarray
    widths: 1-D sequence
        Widths to use for calculating the CWT matrix. In general,
        this range should cover the expected width of peaks of interest.
    wavelet: function
        Should take a single variable and return a 1d array to convolve
        with `vector`. Should be normalized to unit area. Default
        is the ricker wavelet
    max_distances: 1-D ndarray,optional
        Default `widths`/4. See identify_ridge_lines
    gap_thresh: float, optional
        Default 2. See identify_ridge_lines
    min_length: int, optional
        Default None. See filter_ridge_lines
    min_snr: float, optional
        Default 1. See filter_ridge_lines
    noise_perc: float, optional
        Default 10. See filter_ridge_lines

    Notes
    ---------
    This approach was designed for finding sharp peaks among noisy data, however
    with proper parameter selection it should function well for different
    peak shapes.
    The algorithm is as follows:
    1. Perform a continuous wavelet transform on `vector`, for the supplied
    `widths`. This is a convolution of `vector` with `wavelet(width)` for
    each width in `widths`. See `cwt`
    2. Identify "ridge lines" in the cwt matrix. These are relative maxima
    at each row, connected across adjacent rows. See identify_ridge_lines
    3. Filter the ridge_lines using filter_ridge_lines.

    References
    ----------
    Bioinformatics (2006) 22 (17): 2059-2065. doi: 10.1093/bioinformatics/btl355
    http://bioinformatics.oxfordjournals.org/content/22/17/2059.long

    Examples
    --------
    >>> xs = np.arange(0, np.pi, 0.05)
    >>> data = np.sin(xs)
    >>> peakind = find_peaks_cwt(data, np.arange(1,10))
    >>> peakind, xs[peakind],data[peakind]
    ([32], array([ 1.6]), array([ 0.9995736]))
    """
    if gap_thresh is None:
        gap_thresh = np.ceil(widths[0])
    if max_distances is None:
        max_distances = widths / 4.0
    if wavelet is None:
        wavelet = ricker

    cwt_dat = cwt(vector, wavelet, widths)
    ridge_lines = _identify_ridge_lines(cwt_dat, max_distances, gap_thresh)
    filtered = _filter_ridge_lines(cwt_dat, ridge_lines, min_length=min_length,
                                   min_snr=min_snr, noise_perc=noise_perc)
    max_locs = map(lambda x: x[1][0], filtered)
    return sorted(max_locs)
