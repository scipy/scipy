# cython: language_level=3
cimport cython
cimport numpy as cnp
import numpy as np

ctypedef fused np_numeric_t:
    cnp.int8_t
    cnp.int16_t
    cnp.int32_t
    cnp.int64_t
    cnp.uint8_t
    cnp.uint16_t
    cnp.uint32_t
    cnp.uint64_t
    cnp.float32_t
    cnp.float64_t
    cnp.complex64_t
    cnp.complex128_t

__all__ = ['get_array_bandwidth', 'issymmetric', 'ishermitian']


# %% -------------- get_array_bandwidth

@cython.embedsignature(True)
def get_array_bandwidth(a):
    """Return the lower and upper bandwidth of a numeric array.

    Parameters
    ----------
    a : ndarray
        Input array of size (N, M)

    Returns
    -------
    lu: tuple
        2-tuple of ints indicating the lower and upper bandwith. A zero
        denotes no sub- or super-diagonal on that side (triangular), and,
        say for N rows (N-1) means that side is full. Same example applies
        to the upper triangular part with (M-1).

    Notes
    -----
    This helper function simply runs over the array looking for the nonzero
    entries whether there exists a banded structure in the array or not. Hence,
    the performance depends on the density of nonzero entries and also
    memory-layout. Fortran- or C- contiguous arrays are handled best and
    otherwise suffers from extra random memory access cost.

    The strategy is to look for only untested band elements in the upper
    and lower triangular parts separately; depending on the memory layout
    we scan row-wise or column-wise. Moreover, say we are scanning rows
    and in the 6th row, 4th entry is nonzero then, on the succeeding rows
    the horizontal search is done only up to that band entries since we know
    that band is occupied. Therefore, a completely dense matrix scan cost is
    in the the order of n.

    Raises
    ------
    TypeError
        If the dtype of the array is not supported, in particular, NumPy half
        precision floats.
    ValueError
        If the input array is not a 2D NumPy array

    Examples
    --------
    >>> from scipy.linalg import get_array_bandwidth
    >>> A = np.array([[3., 0., 0., 0., 0.],
                      [0., 4., 0., 0., 0.],
                      [0., 0., 5., 1., 0.],
                      [8., 0., 0., 6., 2.],
                      [0., 9., 0., 0., 7.]])
    >>> get_array_bandwidth(A)
    (3, 1)

    """
    if a.size == 0:
        return (0, 0)

    if a.ndim != 2:
        raise ValueError('Input array must be a 2D NumPy array.')

    if a.flags['C_CONTIGUOUS']:
        l, u = get_array_bandwidth_c(a)
    elif a.flags['F_CONTIGUOUS']:
        u, l = get_array_bandwidth_c(a.T)
    else:
        l, u = get_array_bandwidth_noncontig(a)

    return l, u

@cython.initializedcheck(False)
def get_array_bandwidth_c(np_numeric_t[:, ::1]A):
    cdef int l, u
    with nogil:
        l, u = band_check_internal_c(A)
    return l, u

@cython.initializedcheck(False)
def get_array_bandwidth_noncontig(np_numeric_t[:, :]A):
    cdef int l, u
    with nogil:
        l, u = band_check_internal_noncontig(A)
    return l, u


@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline (int, int) band_check_internal_c(np_numeric_t[:, ::1]A) nogil:
    cdef Py_ssize_t n = A.shape[0], lower_band = 0, upper_band = 0, r, c
    cdef np_numeric_t zero = 0

    for r in xrange(n):
        # Only bother if outside the existing band:
        for c in xrange(r-lower_band):
            if A[r, c] != zero:
                lower_band = r - c
                break

        for c in xrange(n - 1, r + upper_band, -1):
            if A[r, c] != zero:
                upper_band = c - r
                break

    return lower_band, upper_band

@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline (int, int) band_check_internal_noncontig(np_numeric_t[:, :]A) nogil:
    cdef Py_ssize_t n = A.shape[0], lower_band = 0, upper_band = 0, r, c
    cdef np_numeric_t zero = 0

    for r in xrange(n):
        # Only bother if outside the existing band:
        for c in xrange(r-lower_band):
            if A[r, c] != zero:
                lower_band = r - c
                break

        for c in xrange(n - 1, r + upper_band, -1):
            if A[r, c] != zero:
                upper_band = c - r
                break

    return lower_band, upper_band
