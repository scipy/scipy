from __future__ import division, absolute_import

cimport cython
cimport numpy as cnp

import os
import numpy as np

# Parameters are linked to the direction numbers list.
# See `initialize_direction_numbers` for more details.
# Declared using DEF to be known at compilation time for ``poly`` et ``vinit``
DEF MAXDIM = 21201  # max number of dimensions
DEF MAXDEG = 18  # max polynomial degree
DEF MAXBIT = 30  # max number of bits


# Needed to be accessed with python
cdef extern from *:
    """
    int MAXDIM_DEFINE = 21201;
    int MAXDEG_DEFINE = 18;
    int MAXBIT_DEFINE = 30;
    """
    int MAXDIM_DEFINE  # max number of dimensions
    int MAXDEG_DEFINE  # max polynomial degree
    int MAXBIT_DEFINE  # max number of bits

_MAXDIM = MAXDIM_DEFINE
_MAXDEG = MAXDEG_DEFINE
_MAXBIT = MAXBIT_DEFINE

cdef int poly[MAXDIM]
cdef int vinit[MAXDIM][MAXDEG]

cdef int LARGEST_NUMBER = 2 ** MAXBIT  # largest possible integer
cdef float RECIPD = 1.0 / LARGEST_NUMBER  # normalization constant

cdef bint is_initialized = False


def initialize_direction_numbers():
    """Load direction numbers.

    Direction numbers obtained using the search criterion D(6)
    up to the dimension 21201. This is the recommended choice by the authors.

    Original data can be found at https://web.maths.unsw.edu.au/~fkuo/sobol/.
    For additional details on the quantities involved, see [1].

    [1] S. Joe and F. Y. Kuo. Remark on algorithm 659: Implementing sobol's
        quasirandom sequence generator. ACM Trans. Math. Softw., 29(1):49-57,
        Mar. 2003.

    The C-code generated from putting the numbers in as literals is obscenely
    large/inefficient. The data file was thus packaged and save as an .npz data
    file for fast loading using the following code (this assumes that the file
    https://web.maths.unsw.edu.au/~fkuo/sobol/new-joe-kuo-6.21201 is present in
    the working directory):

        import pandas as pd
        import numpy as np

        # read in file content
        with open("./new-joe-kuo-6.21201", "r") as f:
            lines = f.readlines()

        rows = []

        # parse data from file line by line
        for l in lines[1:]:
            nums = [int(n) for n in l.replace(" \n", "").split()]
            d, s, a = nums[:3]
            vs = {f"v{i}": int(v) for i,v in enumerate(nums[3:])}
            rows.append({"d": d, "s": s, "a": a, **vs})


        # read in as dataframe, explicitly use zero values
        df = pd.DataFrame(rows).fillna(0).astype(int)

        # peform conversion
        df["poly"] = 2 * df["a"] + 2 ** df["s"] + 1

        # ensure columns are properly ordered
        vs = df[[f"v{i}" for i in range(18)]].values

        # add the degenerate d=1 column (not included in the data file)
        vs = np.vstack([vs[0][np.newaxis, :], vs])
        poly = np.concatenate([[1], df["poly"].values])

        # save as compressed .npz file to minimize size of distribution
        np.savez_compressed("./_sobol_direction_numbers", vinit=vs, poly=poly)

    """
    cdef int[:] dns_poly
    cdef int[:, :] dns_vinit

    global is_initialized
    if not is_initialized:
        dns = np.load(os.path.join(os.path.dirname(__file__), "_sobol_direction_numbers.npz"))
        dns_poly = dns["poly"].astype(np.intc)
        dns_vinit = dns["vinit"].astype(np.intc)
        for i in range(MAXDIM):
            poly[i] = dns_poly[i]
        for i in range(MAXDIM):
            for j in range(MAXDEG):
                vinit[i][j] = dns_vinit[i, j]
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
    cdef int i = 0
    while x & (1 << i) != 0:
        i += 1
    return i + 1


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int ibits(const int x, const int pos, const int length) nogil:
    """Extract a sequence of bits from the bit representation of an integer.

    Extract the sequence from position `pos` (inclusive) to ``pos + length``
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

        # Fill in remaining elements of v as in Section 2 (top of pg. 90) of:
        #
        # P. Bratley and B. L. Fox. Algorithm 659: Implementing sobol's
        # quasirandom sequence generator. ACM Trans.
        # Math. Softw., 14(1):88-100, Mar. 1988.
        #
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
cpdef void _draw(const int n,
                 const int num_gen,
                 const int dim,
                 cnp.int_t[:, :] sv,
                 cnp.int_t[:] quasi,
                 cnp.float_t[:, :] result) nogil:
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
cpdef void _fast_forward(const int n,
                         const int num_gen,
                         const int dim,
                         cnp.int_t[:, :] sv,
                         cnp.int_t[:] quasi) nogil:
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
cpdef void _cscramble(const int dim,
                      cnp.int_t[:, :, :] ltm,
                      cnp.int_t[:, :] sv) nogil:
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
cpdef void _fill_p_cumulative(cnp.float_t[:] p,
                              cnp.float_t[:] p_cumulative) nogil:
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
cpdef void _categorize(cnp.float_t[:] draws,
                       cnp.float_t[:] p_cumulative,
                       cnp.int_t[:] result) nogil:
    cdef int i
    cdef int n_p = p_cumulative.shape[0]
    for i in range(draws.shape[0]):
        j = _find_index(p_cumulative, n_p, draws[i])
        result[j] = result[j] + 1


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int _find_index(cnp.float_t[:] p_cumulative,
                     const int size,
                     const float value) nogil:
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
