"""
Functions for identifying peaks in signals.
"""
from __future__ import division, print_function, absolute_import

import numpy as np

from scipy._lib.six import xrange
from scipy.signal.wavelets import cwt, ricker
from scipy.stats import scoreatpercentile


__all__ = ['argrelmin', 'argrelmax', 'argrelextrema', 'peak_prominences',
           'find_peaks_cwt']


def _boolrelextrema(data, comparator, axis=0, order=1, mode='clip'):
    """
    Calculate the relative extrema of `data`.

    Relative extrema are calculated by finding locations where
    ``comparator(data[n], data[n+1:n+order+1])`` is True.

    Parameters
    ----------
    data : ndarray
        Array in which to find the relative extrema.
    comparator : callable
        Function to use to compare two data points.
        Should take two arrays as arguments.
    axis : int, optional
        Axis over which to select from `data`.  Default is 0.
    order : int, optional
        How many points on each side to use for the comparison
        to consider ``comparator(n,n+x)`` to be True.
    mode : str, optional
        How the edges of the vector are treated.  'wrap' (wrap around) or
        'clip' (treat overflow as the same as the last (or first) element).
        Default 'clip'.  See numpy.take

    Returns
    -------
    extrema : ndarray
        Boolean array of the same shape as `data` that is True at an extrema,
        False otherwise.

    See also
    --------
    argrelmax, argrelmin

    Examples
    --------
    >>> testdata = np.array([1,2,3,2,1])
    >>> _boolrelextrema(testdata, np.greater, axis=0)
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

    Parameters
    ----------
    data : ndarray
        Array in which to find the relative minima.
    axis : int, optional
        Axis over which to select from `data`.  Default is 0.
    order : int, optional
        How many points on each side to use for the comparison
        to consider ``comparator(n, n+x)`` to be True.
    mode : str, optional
        How the edges of the vector are treated.
        Available options are 'wrap' (wrap around) or 'clip' (treat overflow
        as the same as the last (or first) element).
        Default 'clip'. See numpy.take

    Returns
    -------
    extrema : tuple of ndarrays
        Indices of the minima in arrays of integers.  ``extrema[k]`` is
        the array of indices of axis `k` of `data`.  Note that the
        return value is a tuple even when `data` is one-dimensional.

    See Also
    --------
    argrelextrema, argrelmax

    Notes
    -----
    This function uses `argrelextrema` with np.less as comparator.

    .. versionadded:: 0.11.0

    Examples
    --------
    >>> from scipy.signal import argrelmin
    >>> x = np.array([2, 1, 2, 3, 2, 0, 1, 0])
    >>> argrelmin(x)
    (array([1, 5]),)
    >>> y = np.array([[1, 2, 1, 2],
    ...               [2, 2, 0, 0],
    ...               [5, 3, 4, 4]])
    ...
    >>> argrelmin(y, axis=1)
    (array([0, 2]), array([2, 1]))

    """
    return argrelextrema(data, np.less, axis, order, mode)


def argrelmax(data, axis=0, order=1, mode='clip'):
    """
    Calculate the relative maxima of `data`.

    Parameters
    ----------
    data : ndarray
        Array in which to find the relative maxima.
    axis : int, optional
        Axis over which to select from `data`.  Default is 0.
    order : int, optional
        How many points on each side to use for the comparison
        to consider ``comparator(n, n+x)`` to be True.
    mode : str, optional
        How the edges of the vector are treated.
        Available options are 'wrap' (wrap around) or 'clip' (treat overflow
        as the same as the last (or first) element).
        Default 'clip'.  See `numpy.take`.

    Returns
    -------
    extrema : tuple of ndarrays
        Indices of the maxima in arrays of integers.  ``extrema[k]`` is
        the array of indices of axis `k` of `data`.  Note that the
        return value is a tuple even when `data` is one-dimensional.

    See Also
    --------
    argrelextrema, argrelmin

    Notes
    -----
    This function uses `argrelextrema` with np.greater as comparator.

    .. versionadded:: 0.11.0

    Examples
    --------
    >>> from scipy.signal import argrelmax
    >>> x = np.array([2, 1, 2, 3, 2, 0, 1, 0])
    >>> argrelmax(x)
    (array([3, 6]),)
    >>> y = np.array([[1, 2, 1, 2],
    ...               [2, 2, 0, 0],
    ...               [5, 3, 4, 4]])
    ...
    >>> argrelmax(y, axis=1)
    (array([0]), array([1]))
    """
    return argrelextrema(data, np.greater, axis, order, mode)


def argrelextrema(data, comparator, axis=0, order=1, mode='clip'):
    """
    Calculate the relative extrema of `data`.

    Parameters
    ----------
    data : ndarray
        Array in which to find the relative extrema.
    comparator : callable
        Function to use to compare two data points.
        Should take two arrays as arguments.
    axis : int, optional
        Axis over which to select from `data`.  Default is 0.
    order : int, optional
        How many points on each side to use for the comparison
        to consider ``comparator(n, n+x)`` to be True.
    mode : str, optional
        How the edges of the vector are treated.  'wrap' (wrap around) or
        'clip' (treat overflow as the same as the last (or first) element).
        Default is 'clip'.  See `numpy.take`.

    Returns
    -------
    extrema : tuple of ndarrays
        Indices of the maxima in arrays of integers.  ``extrema[k]`` is
        the array of indices of axis `k` of `data`.  Note that the
        return value is a tuple even when `data` is one-dimensional.

    See Also
    --------
    argrelmin, argrelmax

    Notes
    -----

    .. versionadded:: 0.11.0

    Examples
    --------
    >>> from scipy.signal import argrelextrema
    >>> x = np.array([2, 1, 2, 3, 2, 0, 1, 0])
    >>> argrelextrema(x, np.greater)
    (array([3, 6]),)
    >>> y = np.array([[1, 2, 1, 2],
    ...               [2, 2, 0, 0],
    ...               [5, 3, 4, 4]])
    ...
    >>> argrelextrema(y, np.less, axis=1)
    (array([0, 2]), array([2, 1]))

    """
    results = _boolrelextrema(data, comparator,
                              axis, order, mode)
    return np.where(results)


def peak_prominences(x, peaks, wlen=None):
    """
    Calculate the prominence of each peak in a signal.

    The prominence of a peak measures how much a peak stands out from the
    surrounding baseline of the signal and is defined as the vertical
    distance between the peak and its lowest contour line.

    Parameters
    ----------
    x : sequence
        A signal with peaks.
    peaks : sequence
        Indices of peaks in `x`.
    wlen : number, optional
        A window length in samples that limits the search for the lowest
        contour line to a symmetric interval around the evaluated peak. If not
        given the entire signal `x` is used. Use this parameter to speed up the
        calculation significantly for large vectors (see Notes).

    Returns
    -------
    prominences : ndarray
        The calculated prominences for each peak in `peaks`.
    left_bases, right_bases : ndarray
        The peaks' bases as indices in `x` to the left and right of each peak.
        The higher base of each pair is a peak's lowest contour line (see Notes
        for a more details).

    See Also
    --------
    find_peaks, peak_widths

    Notes
    -----
    Strategy to compute a peak's prominence:

    * Extend a horizontal line from the current peak to the left and right until
      the line either reaches the window end (see `wlen`) or intersects the
      signal again at the slope of a higher peak.
    * On each side find the minimal signal value within the interval defined
      above. These points are the peak's bases.
    * The higher one of the two bases marks the peak's lowest contour line. The
      prominence can then be calculated as the vertical difference between the
      peaks height itself and its lowest contour line.

    Searching for the peak's bases can be slow for large `x` because the full
    signal needs to be evaluated for each peak. This evaluation area can be
    limited with the parameter `wlen` which restricts the algorithm to a window
    around the current peak and can shorten the calculation time significantly.
    However this may stop the algorithm from finding the true global contour
    line if the peak's bases are outside this window. Instead a higher contour
    line is found within the restricted window leading to a smaller calculated
    prominence. In practice this is only relevant for the largest set of peaks
    in vector. This behavior may even be used intentionally to calculate "local"
    prominences.

    .. versionadded:: 1.1.0

    References
    ----------
    .. [1] Wikipedia Article for Topographic Prominence:
       https://en.wikipedia.org/wiki/Topographic_prominence

    Examples
    --------
    Create signal test signal with two overlayed harmonic functions

    >>> x = np.linspace(0, 6 * np.pi, 1000)
    >>> x = np.sin(x) + 0.6 * np.sin(2.6 * x)

    Find all peaks and calculate prominences

    >>> from scipy.signal import find_peaks, peak_prominences
    >>> peaks, _ = find_peaks(x)
    >>> prominences = peak_prominences(x, peaks)[0]
    >>> prominences
    array([ 1.24159486,  0.47840168,  0.28470524,  3.10716793,  0.284603  ,
            0.47822491,  2.48340261,  0.47822491])

    Calculate the height of each peak's contour line and plot the results

    >>> contour_heights = x[peaks] - prominences
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(x)
    >>> plt.plot(peaks, x[peaks], "x")
    >>> plt.vlines(peaks, ymin=contour_heights, ymax=x[peaks])
    """
    x = np.asarray(x)
    peaks = np.asarray(peaks)
    if peaks.size == 0:
        # Handle empty peaks
        return np.array([]), np.array([], dtype=int), np.array([], dtype=int)

    if x.ndim != 1:
        raise ValueError('`x` must have exactly one dimension')
    if peaks.ndim != 1:
        raise ValueError('`peaks` must have exactly one dimension')
    if x.size <= peaks.max():
        raise ValueError('an index in `peaks` exceeds the size of `x`')
    if not np.issubdtype(peaks.dtype, np.integer):
        raise ValueError('`peaks` must be an array of integers')
    if wlen is not None and wlen < 3:
            raise ValueError('`wlen` must be at least 3')

    # Prepare return arguments
    prominences = np.zeros(peaks.size)
    left_bases = np.zeros(peaks.size, dtype=int)
    right_bases = np.zeros(peaks.size, dtype=int)

    for i, peak in enumerate(peaks):
        # If wlen is twice the size of x the symmetric window always covers
        # the full x
        if wlen is not None and wlen < x.size * 2:
            wlen = int(wlen)
            # Calculate window borders around the evaluated peak
            wleft = peak - wlen // 2
            wright = peak + wlen // 2
            # Handle border cases
            wleft = 0 if wleft < 0 else wleft
            wright = x.size if x.size < wright else wright
            # Use slice for prominence calculation
            window = x[wleft:wright]
            # Correct peak position in x
            peak -= wleft
        else:
            # Use full x for prominence calculation
            window = x
            wleft = 0

        # Positions where window is larger than current peak height
        greater_peak = np.where(window > window[peak])[0]

        try:
            # Nearest position to the left of peak with
            # window[left] > window[peak]
            left = greater_peak[greater_peak < peak].max()
        except ValueError:
            left = 0
        try:
            # Nearest position to right of peak with
            # window[right] > window[peak]
            right = greater_peak[greater_peak > peak].min()
        except ValueError:
            right = None

        # Base indices to the left and right of peak in window
        left_bases[i] = window[left:peak].argmin() + left
        right_bases[i] = window[peak:right].argmin() + peak

        # Calculate lowest contour and its vertical distance to peak
        lowest_contour = max(window[left_bases[i]], window[right_bases[i]])
        prominences[i] = window[peak] - lowest_contour

        # Correct window offset
        left_bases[i] += wleft
        right_bases[i] += wleft

    return prominences, left_bases, right_bases


def _identify_ridge_lines(matr, max_distances, gap_thresh):
    """
    Identify ridges in the 2-D matrix.

    Expect that the width of the wavelet feature increases with increasing row
    number.

    Parameters
    ----------
    matr : 2-D ndarray
        Matrix in which to identify ridge lines.
    max_distances : 1-D sequence
        At each row, a ridge line is only connected
        if the relative max at row[n] is within
        `max_distances`[n] from the relative max at row[n+1].
    gap_thresh : int
        If a relative maximum is not found within `max_distances`,
        there will be a gap. A ridge line is discontinued if
        there are more than `gap_thresh` points without connecting
        a new relative maximum.

    Returns
    -------
    ridge_lines : tuple
        Tuple of 2 1-D sequences. `ridge_lines`[ii][0] are the rows of the
        ii-th ridge-line, `ridge_lines`[ii][1] are the columns. Empty if none
        found.  Each ridge-line will be sorted by row (increasing), but the
        order of the ridge lines is not specified.

    References
    ----------
    Bioinformatics (2006) 22 (17): 2059-2065.
    :doi:`10.1093/bioinformatics/btl355`
    http://bioinformatics.oxfordjournals.org/content/22/17/2059.long

    Examples
    --------
    >>> data = np.random.rand(5,5)
    >>> ridge_lines = _identify_ridge_lines(data, 1, 1)

    Notes
    -----
    This function is intended to be used in conjunction with `cwt`
    as part of `find_peaks_cwt`.

    """
    if(len(max_distances) < matr.shape[0]):
        raise ValueError('Max_distances must have at least as many rows '
                         'as matr')

    all_max_cols = _boolrelextrema(matr, np.greater, axis=1, order=1)
    # Highest row for which there are any relative maxima
    has_relmax = np.where(all_max_cols.any(axis=1))[0]
    if(len(has_relmax) == 0):
        return []
    start_row = has_relmax[-1]
    # Each ridge line is a 3-tuple:
    # rows, cols,Gap number
    ridge_lines = [[[start_row],
                   [col],
                   0] for col in np.where(all_max_cols[start_row])[0]]
    final_lines = []
    rows = np.arange(start_row - 1, -1, -1)
    cols = np.arange(0, matr.shape[1])
    for row in rows:
        this_max_cols = cols[all_max_cols[row]]

        # Increment gap number of each line,
        # set it to zero later if appropriate
        for line in ridge_lines:
            line[2] += 1

        # XXX These should always be all_max_cols[row]
        # But the order might be different. Might be an efficiency gain
        # to make sure the order is the same and avoid this iteration
        prev_ridge_cols = np.array([line[1][-1] for line in ridge_lines])
        # Look through every relative maximum found at current row
        # Attempt to connect them with existing ridge lines.
        for ind, col in enumerate(this_max_cols):
            # If there is a previous ridge line within
            # the max_distance to connect to, do so.
            # Otherwise start a new one.
            line = None
            if(len(prev_ridge_cols) > 0):
                diffs = np.abs(col - prev_ridge_cols)
                closest = np.argmin(diffs)
                if diffs[closest] <= max_distances[row]:
                    line = ridge_lines[closest]
            if(line is not None):
                # Found a point close enough, extend current ridge line
                line[1].append(col)
                line[0].append(row)
                line[2] = 0
            else:
                new_line = [[row],
                            [col],
                            0]
                ridge_lines.append(new_line)

        # Remove the ridge lines with gap_number too high
        # XXX Modifying a list while iterating over it.
        # Should be safe, since we iterate backwards, but
        # still tacky.
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
    ----------
    cwt : 2-D ndarray
        Continuous wavelet transform from which the `ridge_lines` were defined.
    ridge_lines : 1-D sequence
        Each element should contain 2 sequences, the rows and columns
        of the ridge line (respectively).
    window_size : int, optional
        Size of window to use to calculate noise floor.
        Default is ``cwt.shape[1] / 20``.
    min_length : int, optional
        Minimum length a ridge line needs to be acceptable.
        Default is ``cwt.shape[0] / 4``, ie 1/4-th the number of widths.
    min_snr : float, optional
        Minimum SNR ratio. Default 1. The signal is the value of
        the cwt matrix at the shortest length scale (``cwt[0, loc]``), the
        noise is the `noise_perc`th percentile of datapoints contained within a
        window of `window_size` around ``cwt[0, loc]``.
    noise_perc : float, optional
        When calculating the noise floor, percentile of data points
        examined below which to consider noise. Calculated using
        scipy.stats.scoreatpercentile.

    References
    ----------
    Bioinformatics (2006) 22 (17): 2059-2065. :doi:`10.1093/bioinformatics/btl355`
    http://bioinformatics.oxfordjournals.org/content/22/17/2059.long

    """
    num_points = cwt.shape[1]
    if min_length is None:
        min_length = np.ceil(cwt.shape[0] / 4)
    if window_size is None:
        window_size = np.ceil(num_points / 20)

    window_size = int(window_size)
    hf_window, odd = divmod(window_size, 2)

    # Filter based on SNR
    row_one = cwt[0, :]
    noises = np.zeros_like(row_one)
    for ind, val in enumerate(row_one):
        window_start = max(ind - hf_window, 0)
        window_end = min(ind + hf_window + odd, num_points)
        noises[ind] = scoreatpercentile(row_one[window_start:window_end],
                                        per=noise_perc)

    def filt_func(line):
        if len(line[0]) < min_length:
            return False
        snr = abs(cwt[line[0][0], line[1][0]] / noises[line[1][0]])
        if snr < min_snr:
            return False
        return True

    return list(filter(filt_func, ridge_lines))


def find_peaks_cwt(vector, widths, wavelet=None, max_distances=None,
                   gap_thresh=None, min_length=None, min_snr=1, noise_perc=10):
    """
    Find peaks in a 1-D array with wavelet transformation.

    The general approach is to smooth `vector` by convolving it with
    `wavelet(width)` for each width in `widths`. Relative maxima which
    appear at enough length scales, and with sufficiently high SNR, are
    accepted.

    Parameters
    ----------
    vector : ndarray
        1-D array in which to find the peaks.
    widths : sequence
        1-D array of widths to use for calculating the CWT matrix. In general,
        this range should cover the expected width of peaks of interest.
    wavelet : callable, optional
        Should take two parameters and return a 1-D array to convolve
        with `vector`. The first parameter determines the number of points
        of the returned wavelet array, the second parameter is the scale
        (`width`) of the wavelet. Should be normalized and symmetric.
        Default is the ricker wavelet.
    max_distances : ndarray, optional
        At each row, a ridge line is only connected if the relative max at
        row[n] is within ``max_distances[n]`` from the relative max at
        ``row[n+1]``.  Default value is ``widths/4``.
    gap_thresh : float, optional
        If a relative maximum is not found within `max_distances`,
        there will be a gap. A ridge line is discontinued if there are more
        than `gap_thresh` points without connecting a new relative maximum.
        Default is the first value of the widths array i.e. widths[0].
    min_length : int, optional
        Minimum length a ridge line needs to be acceptable.
        Default is ``cwt.shape[0] / 4``, ie 1/4-th the number of widths.
    min_snr : float, optional
        Minimum SNR ratio. Default 1. The signal is the value of
        the cwt matrix at the shortest length scale (``cwt[0, loc]``), the
        noise is the `noise_perc`th percentile of datapoints contained within a
        window of `window_size` around ``cwt[0, loc]``.
    noise_perc : float, optional
        When calculating the noise floor, percentile of data points
        examined below which to consider noise. Calculated using
        `stats.scoreatpercentile`.  Default is 10.

    Returns
    -------
    peaks_indices : ndarray
        Indices of the locations in the `vector` where peaks were found.
        The list is sorted.

    See Also
    --------
    find_peaks, cwt

    Notes
    -----
    This approach was designed for finding sharp peaks among noisy data,
    however with proper parameter selection it should function well for
    different peak shapes.

    The algorithm is as follows:
     1. Perform a continuous wavelet transform on `vector`, for the supplied
        `widths`. This is a convolution of `vector` with `wavelet(width)` for
        each width in `widths`. See `cwt`
     2. Identify "ridge lines" in the cwt matrix. These are relative maxima
        at each row, connected across adjacent rows. See identify_ridge_lines
     3. Filter the ridge_lines using filter_ridge_lines.

    .. versionadded:: 0.11.0

    References
    ----------
    .. [1] Bioinformatics (2006) 22 (17): 2059-2065.
        :doi:`10.1093/bioinformatics/btl355`
        http://bioinformatics.oxfordjournals.org/content/22/17/2059.long

    Examples
    --------
    >>> from scipy import signal
    >>> xs = np.arange(0, np.pi, 0.05)
    >>> data = np.sin(xs)
    >>> peakind = signal.find_peaks_cwt(data, np.arange(1,10))
    >>> peakind, xs[peakind], data[peakind]
    ([32], array([ 1.6]), array([ 0.9995736]))

    """
    widths = np.asarray(widths)

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
    max_locs = np.asarray([x[1][0] for x in filtered])
    max_locs.sort()

    return max_locs
