"""
Utility functions for finding peaks in signals.
"""

import numpy as np
cimport numpy as np
import cython


__all__ = ['_argmaxima1d']


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
