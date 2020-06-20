from __future__ import division, absolute_import

cimport cython
cimport numpy as cnp

import os
import math
import numpy as np
from scipy.stats import norm

cdef int MAXDIM = 21201  # max number of dimensions
cdef int MAXDEG = 18  # max polynomial degree
cdef int MAXBIT = 30  # max number of bits

cdef int poly[21201]
cdef int vinit[21201][18]

cdef int LARGEST_NUMBER = 2 ** MAXBIT  # largest possible integer
cdef float RECIPD = 1.0 / LARGEST_NUMBER  # normalization constant

cdef bint is_initialized = False


# Load direction numbers (taken from https://web.maths.unsw.edu.au/~fkuo/sobol/)
def initialize_direction_numbers():
    global is_initialized, poly, vinit
    if not is_initialized:
        dns = np.load(os.path.join(os.path.dirname(__file__), "_sobol_direction_numbers.npz"))
        poly = dns["poly"]
        vinit = dns["vinit"]
        is_initialized = True


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int bit_length(const int n):
    cdef int bits = 0
    cdef int nloc = n
    while nloc != 0:
        nloc >>= 1
        bits += 1
    return bits


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int low_0_bit(const int x) nogil:
    """Get the position of the right-most 0 bit for an integer.

    Examples:
        >>> low_0_bit(0)
        1
        >>> low_0_bit(1)
        2
        >>> low_0_bit(2)
        1
        >>> low_0_bit(5)
        2
        >>> low_0_bit(7)
        4

    Parameters
    ----------
        x: int
            An integer.

    Returns
    -------
        position: int
            Position of the right-most 0 bit.

    """
    cdef int z = x
    cdef int i = 0
    while True:
        i += 1
        if z % 2 == 0:
            break
        z = z // 2
    return i


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int ibits(const int x, const int pos, const int length) nogil:
    """Extract a sequence of bits from the bit representation of an integer.

    Extract the sequence from position `pos` (inclusive) to `pos + length`
    (not inclusive), leftwise.

    Examples:
        >>> ibits(1, 0, 1)
        1
        >>> ibits(1, 1, 1)
        0
        >>> ibits(2, 0, 1)
        0
        >>> ibits(2, 0, 2)
        2
        >>> ibits(25, 1, 5)
        12


    Parameters
    ----------
        x: int
            Integer to convert to bit representation.
        pos: int
            Starting position of sequence in bit representation of integer.
        length: int
            Length of sequence (number of bits).

    Returns
    -------
        ibits: int
            Integer value corresponding to bit sequence.

    """
    return (x >> pos) & ((1 << length) - 1)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void initialize_v(cnp.int_t[:, :] v, const int dim):
    cdef int d, i, j, k, m, p, newv, pow2

    if dim == 0:
        return

    # first row of v is all 1s
    for i in range(MAXBIT):
        v[0, i] = 1

    # Remaining rows of v (row 2 through dim, indexed by [1:dim])
    for d in range(1, dim):
        p = poly[d]
        m = bit_length(p) - 1

        # First m elements of row d comes from vinit
        for j in range(m):
            v[d, j] = vinit[d][j]

        # Fill in remaining elements of V per Bratley and Fox, Section 2
        # @TODO: UPDATE
        for j in range(m, MAXBIT):
            newv = v[d, j - m]
            pow2 = 1
            for k in range(m):
                pow2 = pow2 << 1
                if (p >> (m - 1 - k)) & 1:
                    newv = newv ^ (pow2 * v[d, j - k - 1])
            v[d, j] = newv

    # Multiply each column of v by power of 2:
    # v * [2^(maxbit-1), 2^(maxbit-2),..., 2, 1]
    pow2 = 1
    for d in range(MAXBIT):
        for i in range(dim):
            v[i, MAXBIT - 1 - d] *= pow2
        pow2 = pow2 << 1


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void _draw(
    const int n,
    const int num_gen,
    const int dim,
    cnp.int_t[:, :] sv,
    cnp.int_t[:] quasi,
    cnp.float_t[:, :] result,
) nogil:
    cdef int i, j, l, qtmp
    cdef int num_gen_loc = num_gen
    for i in range(n):
        l = low_0_bit(num_gen_loc)
        for j in range(dim):
            qtmp = quasi[j] ^ sv[j, l - 1]
            quasi[j] = qtmp
            result[i, j] = qtmp * RECIPD
        num_gen_loc += 1


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void _fast_forward(
    const int n,
    const int num_gen,
    const int dim,
    cnp.int_t[:, :] sv,
    cnp.int_t[:] quasi,
) nogil:
    cdef int i, j, l
    cdef int num_gen_loc = num_gen
    for i in range(n):
        l = low_0_bit(num_gen_loc)
        for j in range(dim):
            quasi[j] = quasi[j] ^ sv[j, l - 1]
        num_gen_loc += 1


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int cdot_pow2(cnp.int_t[:] a) nogil:
    cdef int i
    cdef int size = a.shape[0]
    cdef int z = 0
    cdef int pow2 = 1
    for i in range(size):
        z += a[size - 1 - i] * pow2
        pow2 *= 2
    return z


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void _cscramble(
    const int dim,
    cnp.int_t[:, :, :] ltm,
    cnp.int_t[:, :] sv,
) nogil:
    cdef int d, i, j, k, l, lsm, lsmdp, p, t1, t2, vdj

    # Set diagonals of maxbit x maxbit arrays to 1
    for d in range(dim):
        for i in range(MAXBIT):
            ltm[d, i, i] = 1

    for d in range(dim):
        for j in range(MAXBIT):
            vdj = sv[d, j]
            l = 1
            t2 = 0
            for p in range(MAXBIT - 1, -1, -1):
                lsmdp = cdot_pow2(ltm[d, p, :])
                t1 = 0
                for k in range(MAXBIT):
                    t1 += ibits(lsmdp, k, 1) * ibits(vdj, k, 1)
                t1 = t1 % 2
                t2 = t2 + t1 * l
                l = 2 * l
            sv[d, j] = t2


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void _fill_p_cumulative(
    cnp.float_t[:] p,
    cnp.float_t[:] p_cumulative,
) nogil:
    cdef int i
    cdef int len_p = p.shape[0]
    cdef float tot = 0
    cdef float t
    for i in range(len_p):
        t = tot + p[i]
        p_cumulative[i] = t
        tot = t


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void _categorize(
    cnp.float_t[:] draws,
    cnp.float_t[:] p_cumulative,
    cnp.int_t[:] result,
) nogil:
    cdef int i
    cdef int n_p = p_cumulative.shape[0]
    for i in range(draws.shape[0]):
        j = _find_index(p_cumulative, n_p, draws[i])
        result[j] = result[j] + 1


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int _find_index(
    cnp.float_t[:] p_cumulative,
    const int size,
    const float value,
) nogil:
    cdef int l = 0
    cdef int r = size - 1
    cdef int m
    while r > l:
        m = (l + r) // 2
        if value > p_cumulative[m]:
            l = m + 1
        else:
            r = m
    return r


def _test_find_index(p_cumulative, size, value):
    # type: (np.ndarray, int, float) -> int
    """Wrapper for testing in python"""
    return _find_index(p_cumulative, size, value)
