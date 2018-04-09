"""
Utility functions for finding peaks in signals.
"""

import numpy as np
import cython

cimport numpy as np
from libc.math cimport ceil


__all__ = ['_argmaxima1d', '_select_by_peak_distance', '_peak_prominences']


@cython.wraparound(False)
@cython.boundscheck(False)
def _argmaxima1d(np.float64_t[:] x not None):
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
      surrounded on both sides by atleast one smaller sample.

    .. versionadded:: 1.1.0
    """
    # Preallocate, there can't be more maxima than half the size of `x`
    cdef np.ndarray[np.intp_t, ndim=1] maxima
    maxima = np.empty(x.shape[0] // 2, dtype=np.intp)
    cdef Py_ssize_t m = 0  # Pointer to the end of valid area in `maxima`

    # Variables to loop over `x`
    cdef Py_ssize_t i = 1  # Pointer to current sample, first one can't be maxima
    cdef Py_ssize_t i_max = x.shape[0] - 1  # Last sample can't be maxima
    cdef Py_ssize_t i_ahead  # Pointer to look ahead of current sample

    while i < i_max:
        # Test if previous sample is smaller
        if x[i - 1] < x[i]:
            i_ahead = i + 1

            # Find next sample that is unequal to x[i]
            while i_ahead < i_max and x[i_ahead] == x[i]:
                i_ahead += 1

            # Maxima is found if next unequal sample is smaller than x[i]
            if x[i_ahead] < x[i]:
                # Store sample in the center of flat area (round down)
                maxima[m] = (i + i_ahead - 1) // 2
                m += 1
                # Skip samples that can't be maxima
                i = i_ahead
        i += 1

    maxima.resize(m, refcheck=False)  # Keep only valid part of array memory.
    return maxima


@cython.wraparound(False)
@cython.boundscheck(False)
def _select_by_peak_distance(np.intp_t[:] peaks not None,
                             np.float64_t[:] priority not None,
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
        np.intp_t[::1] priority_to_position
        np.int8_t[::1] keep
        np.intp_t i, j, k, peaks_size, distance_

    peaks_size = peaks.shape[0]
    # Round up because actual peak distance can only be natural number
    distance_ = <np.intp_t>ceil(distance)
    keep = np.ones(peaks_size, dtype=np.int8)  # Prepare array of flags

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


@cython.wraparound(False)
@cython.boundscheck(False)
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
        If an index in `peaks` doesn't point to a local maximum in `x`.

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
        bint raise_error

    raise_error = False
    prominences = np.empty(peaks.shape[0], dtype=np.float64)
    left_bases = np.empty(peaks.shape[0], dtype=np.intp)
    right_bases = np.empty(peaks.shape[0], dtype=np.intp)

    with nogil:
        for peak_nr in range(peaks.shape[0]):
            peak = peaks[peak_nr]
            i_min = 0
            i_max = x.shape[0] - 1

            if 2 <= wlen:
                # Adjust window around the evaluated peak (within bounds);
                # if wlen is even the resulting window length is is implicitly
                # rounded to next odd integer
                i_min = max(peak - wlen // 2, i_min)
                i_max = min(peak + wlen // 2, i_max)

            # Find the left base in interval [i_min, peak]
            i = peak
            left_min = x[peak]
            while i_min <= i and x[i] <= x[peak]:
                if x[i] < left_min:
                    left_min = x[i]
                    left_bases[peak_nr] = i
                i -= 1
            if not left_min < x[peak]:
                raise_error = True  # Raise error outside nogil statement
                break

            # Find the right base in interval [peak, i_max]
            i = peak
            right_min = x[peak]
            while i <= i_max and x[i] <= x[peak]:
                if x[i] < right_min:
                    right_min = x[i]
                    right_bases[peak_nr] = i
                i += 1
            if not right_min < x[peak]:
                raise_error = True  # Raise error outside nogil statement
                break

            prominences[peak_nr] = x[peak] - max(left_min, right_min)

    if raise_error:
        raise ValueError(str(peak) + ' is not a valid peak')

    # Return memoryviews as ndarrays
    return prominences.base, left_bases.base, right_bases.base
