#cython: wraparound=False
#cython: boundscheck=False
#cython: nonecheck=False

"""Utility functions for finding peaks in signals."""

import warnings

import numpy as np

cimport numpy as np
from libc.math cimport ceil


__all__ = ['_argmaxima1d', '_select_by_peak_distance', '_peak_prominences',
           '_peak_widths']


def _argmaxima1d(np.float64_t[::1] x not None):
    """
    Find indices of local maxima in a 1D array.

    This function finds all local maxima in a 1D array and returns their
    indices. For maxima who are wider than one sample the index of the center
    sample is returned (rounded down in case the number of samples is even).

    Parameters
    ----------
    x : ndarray
        The array to search for local maxima.

    Returns
    -------
    maxima : ndarray
        Indices of local maxima in `x`.

    See Also
    --------
    argrelmax

    Notes
    -----
    - Compared to `argrelmax` this function is significantly faster and can
      detect maxima that are more than one sample wide. However this comes at
      the cost of being only applicable to 1D arrays.
    - A maxima is defined as one or more samples of equal value that are
      surrounded on both sides by at least one smaller sample.

    .. versionadded:: 1.1.0
    """
    cdef:
        np.intp_t[::1] maxima
        np.intp_t m, i, i_ahead, i_max

    # Preallocate, there can't be more maxima than half the size of `x`
    maxima = np.empty(x.shape[0] // 2, dtype=np.intp)
    m = 0  # Pointer to the end of valid area in `maxima`

    with nogil:
        i = 1  # Pointer to current sample, first one can't be maxima
        i_max = x.shape[0] - 1  # Last sample can't be maxima
        while i < i_max:
            # Test if previous sample is smaller
            if x[i - 1] < x[i]:
                i_ahead = i + 1  # Index to look ahead of current sample

                # Find next sample that is unequal to x[i]
                while i_ahead < i_max and x[i_ahead] == x[i]:
                    i_ahead += 1

                # Maxima is found if next unequal sample is smaller than x[i]
                if x[i_ahead] < x[i]:
                    # Store sample in the center of flat area (round down)
                    maxima[m] = (i + i_ahead - 1) // 2
                    m += 1
                    # Skip samples that can't be maximum
                    i = i_ahead
            i += 1

    maxima.base.resize(m, refcheck=False)  # Keep only valid part of array memory.
    return maxima.base


def _select_by_peak_distance(np.intp_t[::1] peaks not None,
                             np.float64_t[::1] priority not None,
                             np.float64_t distance):
    """
    Evaluate which peaks fulfill the distance condition.

    Parameters
    ----------
    peaks : ndarray
        Indices of peaks in `vector`.
    priority : ndarray
        An array matching `peaks` used to determine priority of each peak. A
        peak with a higher priority value is kept over one with a lower one.
    distance : np.float64
        Minimal distance that peaks must be spaced.

    Returns
    -------
    keep : ndarray[bool]
        A boolean mask evaluating to true where `peaks` fulfill the distance
        condition.

    Notes
    -----
    Declaring the input arrays as C-contiguous doesn't seem to have performance
    advantages.

    .. versionadded:: 1.1.0
    """
    cdef:
        np.uint8_t[::1] keep
        np.intp_t[::1] priority_to_position
        np.intp_t i, j, k, peaks_size, distance_

    peaks_size = peaks.shape[0]
    # Round up because actual peak distance can only be natural number
    distance_ = <np.intp_t>ceil(distance)
    keep = np.ones(peaks_size, dtype=np.uint8)  # Prepare array of flags

    # Create map from `i` (index for `peaks` sorted by `priority`) to `j` (index
    # for `peaks` sorted by position). This allows to iterate `peaks` and `keep`
    # with `j` by order of `priority` while still maintaining the ability to
    # step to neighbouring peaks with (`j` + 1) or (`j` - 1).
    priority_to_position = np.argsort(priority)

    with nogil:
        # Highest priority first -> iterate in reverse order (decreasing)
        for i in range(peaks_size - 1, -1, -1):
            # "Translate" `i` to `j` which points to current peak whose
            # neighbours are to be evaluated
            j = priority_to_position[i]
            if keep[j] == 0:
                # Skip evaluation for peak already marked as "don't keep"
                continue

            k = j - 1
            # Flag "earlier" peaks for removal until minimal distance is exceeded
            while 0 <= k and peaks[j] - peaks[k] < distance_:
                keep[k] = 0
                k -= 1

            k = j + 1
            # Flag "later" peaks for removal until minimal distance is exceeded
            while k < peaks_size and peaks[k] - peaks[j] < distance_:
                keep[k] = 0
                k += 1

    return keep.base.view(dtype=np.bool)  # Return as boolean array


class PeakPropertyWarning(RuntimeWarning):
    """Calculated property of a peak has unexpected value."""
    pass


def _peak_prominences(np.float64_t[::1] x not None,
                      np.intp_t[::1] peaks not None,
                      np.intp_t wlen):
    """
    Calculate the prominence of each peak in a signal.

    Parameters
    ----------
    x : ndarray
        A signal with peaks.
    peaks : ndarray
        Indices of peaks in `x`.
    wlen : np.intp
        A window length in samples (see `peak_prominences`) which is rounded up
        to the nearest odd integer. If smaller than 2 the entire signal `x` is
        used.

    Returns
    -------
    prominences : ndarray
        The calculated prominences for each peak in `peaks`.
    left_bases, right_bases : ndarray
        The peaks' bases as indices in `x` to the left and right of each peak.

    Raises
    ------
    ValueError
        If a value in `peaks` is an invalid index for `x`.

    Warns
    -----
    PeakPropertyWarning
        If a prominence of 0 was calculated for any peak.

    Notes
    -----
    This is the inner function to `peak_prominences`.

    .. versionadded:: 1.1.0
    """
    cdef:
        np.float64_t[::1] prominences
        np.intp_t[::1] left_bases, right_bases
        np.float64_t left_min, right_min
        np.intp_t peak_nr, peak, i_min, i_max, i
        np.uint8_t show_warning

    show_warning = False
    prominences = np.empty(peaks.shape[0], dtype=np.float64)
    left_bases = np.empty(peaks.shape[0], dtype=np.intp)
    right_bases = np.empty(peaks.shape[0], dtype=np.intp)

    with nogil:
        for peak_nr in range(peaks.shape[0]):
            peak = peaks[peak_nr]
            i_min = 0
            i_max = x.shape[0] - 1
            if not i_min <= peak <= i_max:
                with gil:
                    raise ValueError("peak {} is not a valid index for `x`"
                                     .format(peak))

            if 2 <= wlen:
                # Adjust window around the evaluated peak (within bounds);
                # if wlen is even the resulting window length is is implicitly
                # rounded to next odd integer
                i_min = max(peak - wlen // 2, i_min)
                i_max = min(peak + wlen // 2, i_max)

            # Find the left base in interval [i_min, peak]
            i = left_bases[peak_nr] = peak
            left_min = x[peak]
            while i_min <= i and x[i] <= x[peak]:
                if x[i] < left_min:
                    left_min = x[i]
                    left_bases[peak_nr] = i
                i -= 1

            # Find the right base in interval [peak, i_max]
            i = right_bases[peak_nr] = peak
            right_min = x[peak]
            while i <= i_max and x[i] <= x[peak]:
                if x[i] < right_min:
                    right_min = x[i]
                    right_bases[peak_nr] = i
                i += 1

            prominences[peak_nr] = x[peak] - max(left_min, right_min)
            if prominences[peak_nr] == 0:
                show_warning = True

    if show_warning:
        warnings.warn("some peaks have a prominence of 0",
                      PeakPropertyWarning, stacklevel=2)
    # Return memoryviews as ndarrays
    return prominences.base, left_bases.base, right_bases.base


def _peak_widths(np.float64_t[::1] x not None,
                 np.intp_t[::1] peaks not None,
                 np.float64_t rel_height,
                 np.float64_t[::1] prominences not None,
                 np.intp_t[::1] left_bases not None,
                 np.intp_t[::1] right_bases not None):
    """
    Calculate the width of each each peak in a signal.

    Parameters
    ----------
    x : ndarray
        A signal with peaks.
    peaks : ndarray
        Indices of peaks in `x`.
    rel_height : np.float64
        Chooses the relative height at which the peak width is measured as a
        percentage of its prominence (see `peak_widths`).
    prominences : ndarray
        Prominences of each peak in `peaks` as returned by `peak_prominences`.
    left_bases, right_bases : ndarray
        Left and right bases of each peak in `peaks` as returned by
        `peak_prominences`.

    Returns
    -------
    widths : ndarray
        The widths for each peak in samples.
    width_heights : ndarray
        The height of the contour lines at which the `widths` where evaluated.
    left_ips, right_ips : ndarray
        Interpolated positions of left and right intersection points of a
        horizontal line at the respective evaluation height.

    Raises
    ------
    ValueError
        If the supplied prominence data doesn't satisfy the condition
        ``0 <= left_base <= peak <= right_base < x.shape[0]`` for each peak or
        if `peaks`, `left_bases` and `right_bases` don't share the same shape.
        Or if `rel_height` is not at least 0.

    Warnings
    --------
    PeakPropertyWarning
        If a width of 0 was calculated for any peak.

    Notes
    -----
    This is the inner function to `peak_widths`.

    .. versionadded:: 1.1.0
    """
    cdef:
        np.float64_t[::1] widths, width_heights, left_ips, right_ips
        np.float64_t height, left_ip, right_ip
        np.intp_t p, peak, i, i_max, i_min
        np.uint8_t show_warning

    if rel_height < 0:
        raise ValueError('`rel_height` must be greater or equal to 0.0')
    if not (peaks.shape[0] == prominences.shape[0] == left_bases.shape[0]
            == right_bases.shape[0]):
        raise ValueError("arrays in `prominence_data` must have the same shape "
                         "as `peaks`")

    show_warning = False
    widths = np.empty(peaks.shape[0], dtype=np.float64)
    width_heights = np.empty(peaks.shape[0], dtype=np.float64)
    left_ips = np.empty(peaks.shape[0], dtype=np.float64)
    right_ips = np.empty(peaks.shape[0], dtype=np.float64)

    with nogil:
        for p in range(peaks.shape[0]):
            i_min = left_bases[p]
            i_max = right_bases[p]
            peak = peaks[p]
            # Validate bounds and order
            if not 0 <= i_min <= peak <= i_max < x.shape[0]:
                with gil:
                    raise ValueError("prominence data is invalid for peak {}"
                                     .format(peak))
            height = width_heights[p] = x[peak] - prominences[p] * rel_height

            # Find intersection point on left side
            i = peak
            while i_min < i and height < x[i]:
                i -= 1
            left_ip = <np.float64_t>i
            if x[i] < height:
                # Interpolate if true intersection height is between samples
                left_ip += (height - x[i]) / (x[i + 1] - x[i])

            # Find intersection point on right side
            i = peak
            while i < i_max and height < x[i]:
                i += 1
            right_ip = <np.float64_t>i
            if  x[i] < height:
                # Interpolate if true intersection height is between samples
                right_ip -= (height - x[i]) / (x[i - 1] - x[i])

            widths[p] = right_ip - left_ip
            if widths[p] == 0:
                show_warning = True
            left_ips[p] = left_ip
            right_ips[p] = right_ip

    if show_warning:
        warnings.warn("some peaks have a width of 0",
                      PeakPropertyWarning, stacklevel=2)
    return widths.base, width_heights.base, left_ips.base, right_ips.base
