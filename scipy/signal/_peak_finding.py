"""
Functions for identifying peaks in signals.
"""
from __future__ import division, print_function, absolute_import

import numpy as np

from scipy._lib.six import xrange
from scipy.signal.wavelets import cwt, ricker
from scipy.stats import scoreatpercentile

from ._peak_finding_utils import _argmaxima1d


__all__ = ['argrelmin', 'argrelmax', 'argrelextrema', 'peak_prominences',
           'peak_widths', 'find_peaks', 'find_peaks_cwt']


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
    argrelextrema, argrelmax, find_peaks

    Notes
    -----
    This function uses `argrelextrema` with np.less as comparator. Therefore it
    requires a strict inequality on both sides of a value to consider it a
    minimum. This means flat minima (more than one sample wide) are not detected.
    In case of one-dimensional `data` `find_peaks` can be used to detect all
    local minima, including flat ones, by calling it with negated `data`.

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
    argrelextrema, argrelmin, find_peaks

    Notes
    -----
    This function uses `argrelextrema` with np.greater as comparator. Therefore
    it  requires a strict inequality on both sides of a value to consider it a
    maximum. This means flat maxima (more than one sample wide) are not detected.
    In case of one-dimensional `data` `find_peaks` can be used to detect all
    local maxima, including flat ones.

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
    find_peaks
        Find peaks inside a signal based on peak properties.
    peak_widths
        Calculate the width of peaks.

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
    >>> from scipy.signal import find_peaks, peak_prominences
    >>> import matplotlib.pyplot as plt

    Create a test signal with two overlayed harmonics

    >>> x = np.linspace(0, 6 * np.pi, 1000)
    >>> x = np.sin(x) + 0.6 * np.sin(2.6 * x)

    Find all peaks and calculate prominences

    >>> peaks, _ = find_peaks(x)
    >>> prominences = peak_prominences(x, peaks)[0]
    >>> prominences
    array([ 1.24159486,  0.47840168,  0.28470524,  3.10716793,  0.284603  ,
            0.47822491,  2.48340261,  0.47822491])

    Calculate the height of each peak's contour line and plot the results

    >>> contour_heights = x[peaks] - prominences
    >>> plt.plot(x)
    >>> plt.plot(peaks, x[peaks], "x")
    >>> plt.vlines(x=peaks, ymin=contour_heights, ymax=x[peaks])
    >>> plt.show()
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


def peak_widths(x, peaks, rel_height=0.5, prominence_data=None, wlen=None):
    """
    Calculate the width of each each peak in a signal.

    This function calculates the width of a peak in samples at a relative
    distance to the peak's height and prominence.

    Parameters
    ----------
    x : sequence
        A signal with peaks.
    peaks : sequence
        Indices of peaks in `x`.
    rel_height : float, optional
        Chooses the relative height at which the peak width is measured as a
        percentage of its prominence. 1.0 calculates the width of the peak at its
        lowest contour line while 0.5 evaluates at half the prominence height.
        Must be a number greater 0. See notes for further explanation.
    prominence_data : tuple, optional
        A tuple of three arrays matching the output of `peak_prominences` when
        called with the same arguments for `x` and `peaks`. This data is
        calculated internally if not provided (see `wlen`).
    wlen : int, optional
        A window length in samples (see `peak_prominences`). This argument is
        only used if `prominence_data` is not given in which case the missing
        data is calculated using `wlen`.

    Returns
    -------
    widths : ndarray
        The widths for each peak in samples.
    width_heights : ndarray
        The height of the contour lines at which the `widths` where evaluated.
    left_ips, right_ips : ndarray
        Interpolated positions of left and right intersection points of a
        horizontal line at the respective evaluation height.

    See Also
    --------
    find_peaks
        Find peaks inside a signal based on peak properties.
    peak_prominences
        Calculate the prominence of peaks.

    Notes
    -----
    The basic algorithm to calculate a peak's width is as follows:

    * Calculate the evaluation height :math:`h_{eval}` with the formula
      :math:`h_{eval} = h_{Peak} - P \\cdot R`, where :math:`h_{Peak}` is the
      height of the peak itself, :math:`P` is the peak's prominence and :math:`R`
      a positive ratio specified with the argument `rel_height`.
    * Draw a horizontal line at the evaluation height to both sides, starting at
      the peak's current vertical position until the lines either intersect a
      slope, the signal border or cross the vertical position of the peak's
      base (see `peak_prominences` for an definition). For the first case,
      intersection with the signal, the true intersection point is estimated with
      linear interpolation.
    * Calculate the width as the horizontal distance between the intersection
      points on both sides.

    As shown above to calculate a peaks width its prominence must be known. You
    can supply these data yourself with the arguments `prominences`, `left_bases`
    and `right_bases`. Otherwise they are internally calculated using `wlen` if
    supplied (see `peak_prominences`).

    .. versionadded:: 1.1.0

    Examples
    --------
    >>> from scipy.signal import chirp, find_peaks, peak_widths
    >>> import matplotlib.pyplot as plt

    Create a test signal with growing peak widths

    >>> x = np.linspace(0, 500, 500)
    >>> x = abs(chirp(x, 1e-4, x.max(), 1.1e-2)) + 2.0 * x / x.max()

    Find all peaks and calculate their widths at the relative height of 0.5

    >>> peaks, _ = find_peaks(x)
    >>> widths, heights, lpos, rpos = peak_widths(x, peaks, rel_height=0.5)
    >>> widths
    array([77.7462348 , 62.19574776, 45.57709222, 37.902356  , 33.33210357,
           29.81097122])

    Plot signal, peaks and contour lines at which the widths where calculated

    >>> plt.plot(x)
    >>> plt.plot(peaks, x[peaks], "x", color="C1")
    >>> plt.hlines(y=heights, xmin=lpos, xmax=rpos, color="C1")
    >>> plt.show()
    """
    x = np.asarray(x)
    peaks = np.asarray(peaks)

    if peaks.size == 0:
        # Handle empty peaks
        return tuple(np.array([]) for _ in range(4))

    if x.ndim != 1:
        raise ValueError('`x` must have exactly one dimension')
    if peaks.ndim != 1:
        raise ValueError('`peaks` must have exactly one dimension')
    if x.size <= peaks.max():
        raise ValueError('an index in `peaks` exceeds the size of `x`')
    if not np.issubdtype(peaks.dtype, np.integer):
        raise ValueError('`peaks` must be an array of integers')
    if rel_height < 0.0:
        raise ValueError('`rel_height` must be greater or equal 0.0')

    if prominence_data is None:
        # Calculate prominence if not supplied and use wlen if supplied
        prominences, left_bases, right_bases = peak_prominences(x, peaks, wlen)
    else:
        prominences, left_bases, right_bases = prominence_data

    # Calculate evaluation height for each peak
    width_heights = x[peaks] - np.asarray(prominences) * rel_height

    widths = np.zeros(peaks.size)
    left_ips = np.zeros(peaks.size)
    right_ips = np.zeros(peaks.size)
    for i, (peak, height) in enumerate(zip(peaks, width_heights)):

        # Maximal peak width is from base to base
        window = x[left_bases[i]:right_bases[i] + 1]
        peak -= left_bases[i]
        # Positions where `window` is smaller reference height
        is_smaller = np.where(window < height)[0]

        try:
            # Nearest position to the left of peak with
            # x[left] > x[peak]
            left_ip = is_smaller[is_smaller < peak].max()
        except ValueError:
            left_ip = None
        try:
            # Nearest position to right of peak with
            # x[right] > x[peak]
            right_ip = is_smaller[is_smaller > peak].min()
        except ValueError:
            right_ip = None

        # If not at window border (ip is None), interpolate sub-sample position
        # to get reasonable precision for steep slopes, do for both sides
        if left_ip is None:
            left_ip = 0
        else:
            y1, y2 = window[left_ip], window[left_ip + 1]
            left_ip += (height - y1) / (y2 - y1)
        if right_ip is None:
            right_ip = window.size - 1
        else:
            y1, y2 = window[right_ip], window[right_ip - 1]
            right_ip -= (height - y1) / (y2 - y1)

        widths[i] = right_ip - left_ip

        # Correct window offset
        left_ips[i] = left_ip + left_bases[i]
        right_ips[i] = right_ip + left_bases[i]

    return widths, width_heights, left_ips, right_ips


def _unpack_condition_args(interval, x, peaks):
    """
    Parse condition arguments for `find_peaks`.

    Parameters
    ----------
    interval : number or ndarray or sequence
        Either a number or ndarray or a 2-element sequence of the former. The
        first value is always interpreted as `imin` and the second, if supplied,
        as `imax`.
    x : ndarray
        The signal with `peaks`.
    peaks : ndarray
        An array with indices used to reduce `imin` and / or `imax` if those are
        arrays.

    Returns
    -------
    imin, imax : number or ndarray or None
        Minimal and maximal value in `argument`.

    Raises
    ------
    ValueError :
        If interval border is given as array and its size does not match the size
        of `x`.

    Notes
    -----

    .. versionadded:: 1.1.0
    """
    try:
        imin, imax = interval
    except (TypeError, ValueError):
        imin, imax = (interval, None)

    # Reduce arrays if arrays
    if isinstance(imin, np.ndarray):
        if imin.size != x.size:
            raise ValueError('array size of lower interval border must match x')
        imin = imin[peaks]
    if isinstance(imax, np.ndarray):
        if imax.size != x.size:
            raise ValueError('array size of upper interval border must match x')
        imax = imax[peaks]

    return imin, imax


def _select_by_property(peak_properties, pmin, pmax):
    """
    Evaluate where the generic property of peaks confirms to an interval.

    Parameters
    ----------
    peak_properties : ndarray
        An array with properties for each peak.
    pmin : None or number or ndarray
        Lower interval boundary for `peak_properties`. ``None`` is interpreted as
        an open border.
    pmax : None or number or ndarray
        Upper interval boundary for `peak_properties`. ``None`` is interpreted as
        an open border.

    Returns
    -------
    keep : bool
        A boolean mask evaluating to true where `peak_properties` confirms to the
        interval.

    See Also
    --------
    find_peaks

    Notes
    -----

    .. versionadded:: 1.1.0
    """
    keep = np.ones(peak_properties.size, dtype=bool)
    if pmin is not None:
        keep &= (pmin <= peak_properties)
    if pmax is not None:
        keep &= (peak_properties <= pmax)
    return keep


def _select_by_peak_threshold(x, peaks, tmin, tmax):
    """
    Evaluate which peaks fulfill the threshold condition.

    Parameters
    ----------
    x : ndarray
        A one-dimensional array which is indexable by `peaks`.
    peaks : ndarray
        Indices of peaks in `x`.
    tmin, tmax : scalar or ndarray or None
         Minimal and / or maximal required thresholds. If supplied as ndarrays
         their size must match `peaks`. ``None`` is interpreted as an open
         border.

    Returns
    -------
    keep : bool
        A boolean mask evaluating to true where `peaks` fulfill the threshold
        condition.
    left_thresholds, right_thresholds : ndarray
        Array matching `peak` containing the thresholds of each peak on
        both sides.

    Notes
    -----

    .. versionadded:: 1.1.0
    """
    # Stack thresholds on both sides to make min / max operations easier:
    # tmin is compared with the smaller, and tmax with the greater thresold to
    # each peak's side
    stacked_thresholds = np.vstack([x[peaks] - x[peaks - 1],
                                    x[peaks] - x[peaks + 1]])
    keep = np.ones(peaks.size, dtype=bool)
    if tmin is not None:
        min_thresholds = np.min(stacked_thresholds, axis=0)
        keep &= (tmin <= min_thresholds)
    if tmax is not None:
        max_thresholds = np.max(stacked_thresholds, axis=0)
        keep &= (max_thresholds <= tmax)

    return keep, stacked_thresholds[0], stacked_thresholds[1]


# Code for _select_by_peak_distance was adapted from
# https://github.com/demotu/BMC/blob/master/functions/detect_peaks.py
# by Marcos Duarte under the MIT license:
#
#     Copyright (c) 2013 Marcos Duarte
#
#     Permission is hereby granted, free of charge, to any person
#     obtaining a copy of this software and associated documentation
#     files (the "Software"), to deal in the Software without
#     restriction, including without limitation the rights to use,
#     copy, modify, merge, publish, distribute, sublicense, and/or sell
#     copies of the Software, and to permit persons to whom the
#     Software is furnished to do so, subject to the following
#     conditions:
#
#     The above copyright notice and this permission notice shall be
#     included in all copies or substantial portions of the Software.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
#     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
#     OTHER DEALINGS IN THE SOFTWARE.

def _select_by_peak_distance(peaks, priority, dmin):
    """
    Evaluate which peaks fulfill the distance condition.

    Parameters
    ----------
    peaks : ndarray
        Indices of peaks in `vector`.
    priority : ndarray
        An array with priorities matching `peaks` used to determine priority of
        peaks. A peak with a higher priority value is kept over one with a lower
        one.
    dmin : number
        Minimal distance that peaks must be spaced.

    Returns
    -------
    keep : ndarray[bool]
        A boolean mask evaluating to true where `peaks` fulfill the distance
        condition.

    Notes
    -----

    .. versionadded:: 1.1.0
    """
    # Peaks are evaluated by priority (larger first)
    eval_peaks = peaks[np.argsort(priority)][::-1]

    # Flag peaks for deletion
    del_flag = np.zeros(eval_peaks.size, dtype=bool)
    for i in range(eval_peaks.size):
        if not del_flag[i]:
            # Flag peaks in intervall +-distance around current peak
            del_flag |= (eval_peaks > eval_peaks[i] - dmin) \
                        & (eval_peaks < eval_peaks[i] + dmin)
            # Keep current peak
            del_flag[i] = False

    keep = ~del_flag[np.argsort(eval_peaks)]

    return keep


def find_peaks(x, height=None, threshold=None, distance=None,
               prominence=None, width=None, wlen=None, rel_height=0.5):
    """
    Find peaks inside a signal based on peak properties.

    This function takes a one-dimensional array and finds all local maxima by
    simple comparison of neighbouring values. Optionally, a subset of these peaks
    can be selected by specifying conditions for a peak's properties.

    Parameters
    ----------
    x : sequence
        A signal with peaks.
    height : number or ndarray or sequence, optional
        Required height of peaks. Either a number, ``None``, an array matching
        `x` or a 2-element sequence of the former. The first element is
        always interpreted as the  minimum and the second, if supplied, as the
        maximum required height.
    threshold : number or ndarray or sequence, optional
        Required threshold of peaks, the vertical distance to its neighbouring
        samples. Either a number, ``None``, an array matching `x` or a
        2-element sequence of the former. The first element is always
        interpreted as the  minimum and the second, if supplied, as the maximum
        required threshold.
    distance : number, optional
        Required minimal horizontal distance (>= 1) in samples between
        neighbouring peaks. The removal order is explained in the notes section.
    prominence : number or ndarray or sequence, optional
        Required prominence of peaks. Either a number, ``None``, an array
        matching `x` or a 2-element sequence of the former. The first
        element is always interpreted as the  minimum and the second, if
        supplied, as the maximum required prominence.
    width : number or ndarray or sequence, optional
        Required width of peaks in samples. Either a number, ``None``, an array
        matching `x` or a 2-element sequence of the former. The first
        element is always interpreted as the  minimum and the second, if
        supplied, as the maximum required prominence.
    wlen : number, optional
        Used for calculation of the peaks prominences, thus it is only used if
        one of the arguments `prominence` or `width` is given. See argument
        `wlen` in `peak_prominences` for a full description of its effects.
    rel_height : float, optional
        Used for calculation of the peaks width, thus it is only used if `width`
        is given. See argument  `rel_height` in `peak_widths` for a full
        description of its effects.

    Returns
    -------
    peaks : ndarray
        Indices of peaks in `x` that satisfy all given conditions.
    properties : dict
        A dictionary containing properties of the returned peaks which were
        calculated as intermediate results during evaluation of the specified
        conditions:

        * 'peak_heights'
              If `height` is given, the height of each peak in `x`.
        * 'left_thresholds', 'right_thresholds'
              If `threshold` is given, these keys contain a peaks vertical
              distance to its neighbouring samples.
        * 'peak_prominences', 'right_bases', 'left_bases'
              If `prominence` is given, these keys are accessible. See
              `peak_prominences` for a description of their content.
        * 'width_heights', 'left_ips', 'right_ips'
              If `width` is given, these keys are accessible. See `peak_widths`
              for a description of their content.

        To calculate and return properties without excluding peaks, provide the
        open interval ``(None, None)`` as a value to the appropriate argument
        (excluding `distance`).

    See Also
    --------
    find_peaks_cwt
        Find peaks using the wavelet transformation.
    peak_prominences
        Directly calculate the prominence of peaks.
    peak_widths
        Directly calculate the width of peaks.

    Notes
    -----
    Because this function searches for local maxima by direct sample comparison,
    the determined peak locations can be off for noisy signals if the noise
    changes the position of a local maximum. In those cases consider smoothing
    the signal before searching for peaks or using other peak finding and fitting
    methods (like `find_peaks_cwt`).

    Some additional comments on specifying conditions:

    * Almost all conditions (excluding `distance`) can be given as half-open or
      closed intervals, e.g ``1`` or ``(1, None)`` defines the half-open interval
      :math:`[1, \\infty]` while ``(None, 1)`` defines the interval
      :math:`[-\\infty, 1]`. The open interval ``(None, None)`` can be specified
      as well, which returns the matching properties without exclusion of peaks.
    * The border is always included in the interval used to select valid peaks.
    * For several conditions the interval borders can be specified with
      arrays matching `x` in shape which enables dynamic constrains based on
      the sample position.
    * The order of arguments given in the function definition above mirrors the
      actual order in which conditions are evaluated. In most cases this order is
      the fastest one because faster operations are applied first to reduce the
      number of peaks that need to be evaluated later.
    * Satisfying the distance condition is accomplished by iterating over all
      peaks in descending order based on their height and removing all lower
      peaks that are too close. This option can be quite slow if many peaks need
      to be evaluated. Try to reduce the number of peaks beforehand by
      specifying conditions that are evaluated before this one (`height` and
      `threshold`).
    * Use `wlen` to reduce the time it takes to evaluate the conditions for
      `prominence` or `width` if `x` is large or has many local maxima
      (see `peak_prominences`).

    .. versionadded:: 1.1.0

    Examples
    --------
    >>> from scipy.signal import find_peaks
    >>> import matplotlib.pyplot as plt

    Create test signal `x` using 7 harmonics:

    >>> gains = [1, -1, 0.6, 0.5, 0.4, 0.3, 0.1]
    >>> freqs = [2, 3.5, 6, 7.1, 11.1, 12, 20]
    >>> t = np.linspace(0, 6, 1000)
    >>> x = sum(g * np.sin(f * t) for g, f in zip(gains, freqs))

    Find all peaks (local maxima) in `x`

    >>> peaks, _ = find_peaks(x)
    >>> peaks
    array([ 25, 119, 205, 275, 381, 495, 586, 782, 892, 952])

    and plot the results

    >>> plt.figure()
    >>> plt.plot(x)
    >>> plt.plot(peaks, x[peaks], 'x')
    >>> plt.show()

    This time, find peaks that have a minimal prominence of 0.5 and whose peak
    base is not more than 150 samples wide (``width=(None, 150)``). For this we
    need to use the option ``rel_height=1.0`` which will evaluate the width at
    the peak's base.

    >>> peaks, prop = find_peaks(
    ...     x, prominence=0.5, width=(None, 150), rel_height=1.0)
    >>> peaks, prop['prominences'], prop['widths']
    (array([ 25, 381, 892]),
     array([1.40613266, 0.50503469, 1.23546031]),
     array([59.57018554, 57.9368999 , 86.11729375]))

    and plot the results including the calculated peak properties

    >>> plt.figure()
    >>> plt.plot(x)
    >>> plt.plot(peaks, x[peaks], 'x')
    >>> plt.vlines(x=peaks, ymin=prop['width_heights'], ymax=x[peaks])
    >>> plt.hlines(y=prop['width_heights'], xmin=prop['left_ips'],
    ...            xmax=prop['right_ips'])
    >>> plt.show()
    """
    # _argmaxima1d expects array of dtype 'float64'
    x = np.asarray(x, dtype=np.float64)
    if x.ndim != 1:
        raise ValueError('`x` must have exactly one dimension')
    if distance is not None and distance < 1:
        raise ValueError('`distance` must be greater or equal to 1')

    peaks = _argmaxima1d(x)
    properties = {}

    if height is not None:
        # Evaluate height condition
        peak_heights = x[peaks]
        hmin, hmax = _unpack_condition_args(height, x, peaks)
        keep = _select_by_property(peak_heights, hmin, hmax)
        peaks = peaks[keep]
        properties["peak_heights"] = peak_heights[keep]

    if threshold is not None:
        # Evaluate threshold condition
        tmin, tmax = _unpack_condition_args(threshold, x, peaks)
        keep, left_thresholds, right_thresholds = _select_by_peak_threshold(
            x, peaks, tmin, tmax)
        peaks = peaks[keep]
        properties["left_thresholds"] = left_thresholds
        properties["right_thresholds"] = right_thresholds
        properties = {key: array[keep] for key, array in properties.items()}

    if distance is not None:
        # Evaluate distance condition
        keep = _select_by_peak_distance(peaks, x[peaks], distance)
        peaks = peaks[keep]
        properties = {key: array[keep] for key, array in properties.items()}

    if prominence is not None or width is not None:
        # Calculate prominence (required for both conditions)
        properties.update(zip(
            ['prominences', 'left_bases', 'right_bases'],
            peak_prominences(x, peaks, wlen=wlen)
        ))

    if prominence is not None:
        # Evaluate prominence condition
        pmin, pmax = _unpack_condition_args(prominence, x, peaks)
        keep = _select_by_property(properties['prominences'], pmin, pmax)
        peaks = peaks[keep]
        properties = {key: array[keep] for key, array in properties.items()}

    if width is not None:
        # Calculate widths
        properties.update(zip(
            ['widths', 'width_heights', 'left_ips', 'right_ips'],
            peak_widths(x, peaks, rel_height, (properties['prominences'],
                                               properties['left_bases'],
                                               properties['right_bases']))
        ))
        # Evaluate width condition
        wmin, wmax = _unpack_condition_args(width, x, peaks)
        keep = _select_by_property(properties['widths'], wmin, wmax)
        peaks = peaks[keep]
        properties = {key: array[keep] for key, array in properties.items()}

    return peaks, properties


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
    cwt
        Continuous wavelet transform.
    find_peaks
        Find peaks inside a signal based on peak properties.

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
