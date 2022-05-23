# cython: language_level=3
cimport cython
cimport numpy as cnp
import numpy as np
from scipy.linalg._cythonized_array_utils cimport (
    lapack_t,
	np_complex_numeric_t,
	np_numeric_t
    )

__all__ = ['bandwidth', 'issymmetric', 'ishermitian']


# ====================== swap_c_and_f_layout : s, d, c, z ====================
@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.initializedcheck(False)
cdef void swap_c_and_f_layout(lapack_t *a, lapack_t *b, int r, int c, int n) nogil:
    """Recursive matrix transposition for square arrays"""
    cdef int i, j, ith_row, r2, c2
    cdef lapack_t *bb=b
    cdef lapack_t *aa=a
    if r < 16:
        for j in range(c):
            ith_row = 0
            for i in range(r):
            # Basically b[i*n+j] = a[j*n+i] without index math
                bb[ith_row] = aa[i]
                ith_row += n
            aa += n
            bb += 1
    else:
        # If tall
        if (r > c):
            r2 = r//2
            swap_c_and_f_layout(a, b, r2, c, n)
            swap_c_and_f_layout(a + r2, b+(r2)*n, r-r2, c, n)
        else:  # Nope
            c2 = c//2
            swap_c_and_f_layout(a, b, r, c2, n);
            swap_c_and_f_layout(a+(c2)*n, b+c2, r, c-c2, n)
# ============================================================================


@cython.embedsignature(True)
def bandwidth(a):
    """Return the lower and upper bandwidth of a 2D numeric array.

    Parameters
    ----------
    a : ndarray
        Input array of size (N, M)

    Returns
    -------
    lu : tuple
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
        If the dtype of the array is not supported, in particular, NumPy
        float16, float128 and complex256 dtypes.

    Examples
    --------
    >>> from scipy.linalg import bandwidth
    >>> A = np.array([[3., 0., 0., 0., 0.],
    ...               [0., 4., 0., 0., 0.],
    ...               [0., 0., 5., 1., 0.],
    ...               [8., 0., 0., 6., 2.],
    ...               [0., 9., 0., 0., 7.]])
    >>> bandwidth(A)
    (3, 1)

    """
    if a.size == 0:
        return (0, 0)

    if a.ndim != 2:
        raise ValueError('Input array must be a 2D NumPy array.')

    if a.flags['C_CONTIGUOUS']:
        l, u = bandwidth_c(a)
    elif a.flags['F_CONTIGUOUS']:
        u, l = bandwidth_c(a.T)
    else:
        l, u = bandwidth_noncontig(a)

    return l, u


@cython.initializedcheck(False)
def bandwidth_c(np_numeric_t[:, ::1]A):
    cdef int l, u
    with nogil:
        l, u = band_check_internal_c(A)
    return l, u


@cython.initializedcheck(False)
def bandwidth_noncontig(np_numeric_t[:, :]A):
    cdef int l, u
    with nogil:
        l, u = band_check_internal_noncontig(A)
    return l, u


@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline (int, int) band_check_internal_c(np_numeric_t[:, ::1]A) nogil:
    cdef Py_ssize_t n = A.shape[0], m = A.shape[1]
    cdef Py_ssize_t lower_band = 0, upper_band = 0, r, c
    cdef np_numeric_t zero = 0

    # lower triangular part
    for r in range(n-1, 0, -1):
        # Only bother if outside the existing band:
        for c in range(min(r-lower_band, m - 1)):
            if A[r, c] != zero:
                lower_band = r - c
                break
        # If the first element is full don't go higher we are done
        if lower_band == r:
            break

    # upper triangular part
    for r in range(n-1):
        for c in range(m - 1, r + upper_band, -1):
            if A[r, c] != zero:
                upper_band = c - r
                break
        if upper_band == c:
            break

    return lower_band, upper_band


@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline (int, int) band_check_internal_noncontig(np_numeric_t[:, :]A) nogil:
    cdef Py_ssize_t n = A.shape[0], m = A.shape[1]
    cdef Py_ssize_t lower_band = 0, upper_band = 0, r, c
    cdef np_numeric_t zero = 0

    # lower triangular part
    for r in range(n-1, 0, -1):
        # Only bother if outside the existing band:
        for c in range(min(r-lower_band, m - 1)):
            if A[r, c] != zero:
                lower_band = r - c
                break
        # If the first element is full don't go higher; we are done
        if lower_band == r:
            break

    # upper triangular part
    for r in range(n-1):
        for c in range(m - 1, r + upper_band, -1):
            if A[r, c] != zero:
                upper_band = c - r
                break
        if upper_band == c:
            break

    return lower_band, upper_band


@cython.embedsignature(True)
def issymmetric(a, atol=None, rtol=None):
    """Check if a square 2D array is symmetric.

    Parameters
    ----------
    a : ndarray
        Input array of size (N, N).

    atol : float, optional
        Absolute error bound

    rtol : float, optional
        Relative error bound

    Returns
    -------
    sym : bool
        Returns True if the array symmetric.

    Notes
    -----
    For square empty arrays the result is returned True by convention. Complex
    valued arrays are tested for symmetricity and not for being Hermitian (see
    examples)

    The diagonal of the array is not scanned. Thus if there are infs, NaNs or
    similar problematic entries on the diagonal, they will be ignored. However,
    `numpy.inf` will be treated as a number, that is to say ``[[1, inf],
    [inf, 2]]`` will return ``True``. On the other hand `numpy.NaN` is never
    symmetric, say, ``[[1, nan], [nan, 2]]`` will return ``False``.

    When ``atol`` and/or ``rtol`` are set to , then the comparison is performed
    by `numpy.allclose` and the tolerance values are passed to it. Otherwise an
    exact comparison against zero is performed by internal functions. Hence
    performance can improve or degrade depending on the size and dtype of the
    array. If one of ``atol`` or ``rtol`` given the other one is automatically
    set to zero.

    See Also
    --------
    ishermitian : Check if a square 2D array is Hermitian

    Raises
    ------
    TypeError
        If the dtype of the array is not supported, in particular, NumPy
        float16, float128 and complex256 dtypes for exact comparisons.

    Examples
    --------
    >>> from scipy.linalg import issymmetric
    >>> A = np.arange(9).reshape(3, 3)
    >>> A = A + A.T
    >>> issymmetric(A)
    True
    >>> Ac = np.array([[1. + 1.j, 3.j], [3.j, 2.]])
    >>> issymmetric(Ac)  # not Hermitian but symmetric
    True

    """
    if a.ndim != 2:
        raise ValueError('Input array must be a 2D NumPy array.')
    if not np.equal(*a.shape):
        raise ValueError('Input array must be square.')
    if a.size == 0:
        return True

    # It's not worth going element-by-element basis if comparison is inexact
    # Also integers do not have tolerances
    if (atol or rtol) and not np.issubdtype(a.dtype, np.integer):
        # We ended up here because one of atol/rtol is given
        # Don't send None to allclose if not provided; replace with 0.
        return np.allclose(a, a.T,
                           atol=atol if atol else 0.,
                           rtol=rtol if rtol else 0.)

    if a.flags['C_CONTIGUOUS']:
        s = is_sym_her_real_c(a)
    elif a.flags['F_CONTIGUOUS']:
        s = is_sym_her_real_c(a.T)
    else:
        s = is_sym_her_real_noncontig(a)

    return s


@cython.initializedcheck(False)
def is_sym_her_real_c(np_numeric_t[:, ::1]A):
    cdef bint s
    with nogil:
        s = is_sym_her_real_c_internal(A)
    return s


@cython.initializedcheck(False)
def is_sym_her_real_noncontig(np_numeric_t[:, :]A):
    cdef bint s
    with nogil:
        s = is_sym_her_real_noncontig_internal(A)
    return s


@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline bint is_sym_her_real_c_internal(np_numeric_t[:, ::1]A) nogil:
    cdef Py_ssize_t n = A.shape[0], r, c

    for r in xrange(n):
        for c in xrange(r):
            if A[r, c] != A[c, r]:
                return False
    return True


@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline bint is_sym_her_real_noncontig_internal(np_numeric_t[:, :]A) nogil:
    cdef Py_ssize_t n = A.shape[0], r, c

    for r in xrange(n):
        for c in xrange(r):
            if A[r, c] != A[c, r]:
                return False
    return True


@cython.embedsignature(True)
def ishermitian(a, atol=None, rtol=None):
    """Check if a square 2D array is Hermitian.

    Parameters
    ----------
    a : ndarray
        Input array of size (N, N)

    atol : float, optional
        Absolute error bound

    rtol : float, optional
        Relative error bound

    Returns
    -------
    her : bool
        Returns True if the array Hermitian.

    Notes
    -----
    For square empty arrays the result is returned True by convention.

    `numpy.inf` will be treated as a number, that is to say ``[[1, inf],
    [inf, 2]]`` will return ``True``. On the other hand `numpy.NaN` is never
    symmetric, say, ``[[1, nan], [nan, 2]]`` will return ``False``.

    When ``atol`` and/or ``rtol`` are set to , then the comparison is performed
    by `numpy.allclose` and the tolerance values are passed to it. Otherwise an
    exact comparison against zero is performed by internal functions. Hence
    performance can improve or degrade depending on the size and dtype of the
    array. If one of ``atol`` or ``rtol`` given the other one is automatically
    set to zero.

    Raises
    ------
    TypeError
        If the dtype of the array is not supported, in particular, NumPy
        float16, float128 and complex256 dtypes.

    See Also
    --------
    issymmetric : Check if a square 2D array is symmetric

    Examples
    --------
    >>> from scipy.linalg import ishermitian
    >>> A = np.arange(9).reshape(3, 3)
    >>> A = A + A.T
    >>> ishermitian(A)
    True
    >>> A = np.array([[1., 2. + 3.j], [2. - 3.j, 4.]])
    >>> ishermitian(A)
    True
    >>> Ac = np.array([[1. + 1.j, 3.j], [3.j, 2.]])
    >>> ishermitian(Ac)  # not Hermitian but symmetric
    False
    >>> Af = np.array([[0, 1 + 1j], [1 - (1+1e-12)*1j, 0]])
    >>> ishermitian(Af)
    False
    >>> ishermitian(Af, atol=5e-11) # almost hermitian with atol
    True

    """
    if a.ndim != 2:
        raise ValueError('Input array must be a 2D NumPy array.')
    if not np.equal(*a.shape):
        raise ValueError('Input array must be square.')
    if a.size == 0:
        return True

    # It's not worth going element-by-element basis if comparison is inexact
    # Also integers do not have tolerances
    if (atol or rtol) and not np.issubdtype(a.dtype, np.integer):
        # We ended up here because one of atol/rtol is given
        # Don't send None to allclose if not provided; replace with 0.
        return np.allclose(a, a.conj().T,
                           atol=atol if atol else 0.,
                           rtol=rtol if rtol else 0.)


    if np.iscomplexobj(a):
        # complex entries on the diagonal
        if a.flags['C_CONTIGUOUS']:
            s = is_sym_her_complex_c(a)
        elif a.flags['F_CONTIGUOUS']:
            s = is_sym_her_complex_c(a.T)
        else:
            s = is_sym_her_complex_noncontig(a)

    else:  # real branch; delegate to issymmetric
        if a.flags['C_CONTIGUOUS']:
            s = is_sym_her_real_c(a)
        elif a.flags['F_CONTIGUOUS']:
            s = is_sym_her_real_c(a.T)
        else:
            s = is_sym_her_real_noncontig(a)

    return s


@cython.initializedcheck(False)
def is_sym_her_complex_c(np_complex_numeric_t[:, ::1]A):
    cdef bint s
    with nogil:
        s = is_sym_her_complex_c_internal(A)
    return s

@cython.initializedcheck(False)
def is_sym_her_complex_noncontig(np_complex_numeric_t[:, :]A):
    cdef bint s
    with nogil:
        s = is_sym_her_complex_noncontig_internal(A)
    return s

@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline bint is_sym_her_complex_c_internal(np_complex_numeric_t[:, ::1]A) nogil:
    cdef Py_ssize_t n = A.shape[0], r, c

    for r in xrange(n):
        for c in xrange(r+1):
            if A[r, c] != A[c, r].conjugate():
                return False
    return True

@cython.initializedcheck(False)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline bint is_sym_her_complex_noncontig_internal(np_complex_numeric_t[:, :]A) nogil:
    cdef Py_ssize_t n = A.shape[0], r, c

    for r in xrange(n):
        for c in xrange(r+1):
            if A[r, c] != A[c, r].conjugate():
                return False
    return True
