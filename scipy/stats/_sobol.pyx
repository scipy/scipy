# cython: language_level=3
# cython: cdivision=True
from __future__ import division, absolute_import

cimport cython
cimport numpy as cnp

import os
import numpy as np

cnp.import_array()

# Parameters are linked to the direction numbers list.
# See `_initialize_direction_numbers` for more details.
# Declared using DEF to be known at compilation time for ``poly`` et ``vinit``
DEF MAXDIM = 21201  # max number of dimensions
DEF MAXDEG = 18  # max polynomial degree


ctypedef fused uint_32_64:
    cnp.uint32_t
    cnp.uint64_t


# Needed to be accessed with python
cdef extern from *:
    """
    int MAXDIM_DEFINE = 21201;
    int MAXDEG_DEFINE = 18;
    """
    int MAXDIM_DEFINE  # max number of dimensions
    int MAXDEG_DEFINE  # max polynomial degree

_MAXDIM = MAXDIM_DEFINE
_MAXDEG = MAXDEG_DEFINE

_poly_dict = {}
_vinit_dict = {}


def get_poly_vinit(kind, dtype):
    """Initialize and cache the direction numbers.

    Uses a dictionary to store the arrays. `kind` allows to select which
    dictionary to pull. The key of each dictionary corresponds to the `dtype`.
    If the key is not present in any of the dictionary, both dictionaries are
    initialized with `_initialize_direction_numbers`, for the given `dtype`.

    This is only used during the initialization step in `_initialize_v`.

    Parameters
    ----------
    kind : {'poly', 'vinit'}
        Select which dictionary to pull.
    dtype : {np.uint32, np.uint64}
        Which dtype to use.

    Returns
    -------
    poly_vinit : np.ndarray
        Either ``poly`` or ``vinit`` matrix.

    """
    if kind == 'poly':
        poly_vinit = _poly_dict.get(dtype)
    else:
        poly_vinit = _vinit_dict.get(dtype)

    if poly_vinit is None:
        _poly_dict[dtype] = np.empty((MAXDIM,), dtype=dtype)
        _vinit_dict[dtype] = np.empty((MAXDIM, MAXDEG), dtype=dtype)

        _initialize_direction_numbers(_poly_dict[dtype], _vinit_dict[dtype], dtype)

        if kind == 'poly':
            poly_vinit = _poly_dict.get(dtype)
        else:
            poly_vinit = _vinit_dict.get(dtype)

    return poly_vinit


def _initialize_direction_numbers(poly, vinit, dtype):
    """Load direction numbers into two arrays.

    Parameters
    ----------
    poly, vinit : np.ndarray
        Direction numbers arrays to fill.
    dtype : {np.uint32, np.uint64}
        Which dtype to use.

    Notes
    -----
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

        # perform conversion
        df["poly"] = 2 * df["a"] + 2 ** df["s"] + 1

        # ensure columns are properly ordered
        vs = df[[f"v{i}" for i in range(18)]].values

        # add the degenerate d=1 column (not included in the data file)
        vs = np.vstack([vs[0][np.newaxis, :], vs])
        poly = np.concatenate([[1], df["poly"].values])

        # save as compressed .npz file to minimize size of distribution
        np.savez_compressed("./_sobol_direction_numbers", vinit=vs, poly=poly)

    """
    dns = np.load(os.path.join(os.path.dirname(__file__),
                  "_sobol_direction_numbers.npz"))
    dns_poly = dns["poly"].astype(dtype)
    dns_vinit = dns["vinit"].astype(dtype)
    poly[...] = dns_poly
    vinit[...] = dns_vinit


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int bit_length(uint_32_64 n):
    cdef int bits = 0
    cdef uint_32_64 nloc = n
    while nloc != 0:
        nloc >>= 1
        bits += 1
    return bits


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int low_0_bit(uint_32_64 x) nogil:
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
    x : int
        An integer.

    Returns
    -------
    position : int
        Position of the right-most 0 bit.

    """
    cdef int i = 0
    while x & (1 << i) != 0:
        i += 1
    return i + 1


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int ibits(uint_32_64 x, const int pos, const int length) nogil:
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
    x : int
        Integer to convert to bit representation.
    pos : int
        Starting position of sequence in bit representation of integer.
    length : int
        Length of sequence (number of bits).

    Returns
    -------
    ibits : int
        Integer value corresponding to bit sequence.

    """
    return (x >> pos) & ((1 << length) - 1)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void _initialize_v(
    uint_32_64[:, ::1] v, const int dim, const int bits
):
    """Initialize matrix of size ``dim * bits`` with direction numbers."""
    cdef int d, i, j, k, m
    cdef uint_32_64 p, newv, pow2
    cdef uint_32_64[:] poly = get_poly_vinit(
        'poly',
        np.uint32 if uint_32_64 is cnp.uint32_t else np.uint64
    )
    cdef uint_32_64[:, ::1] vinit = get_poly_vinit(
        'vinit',
        np.uint32 if uint_32_64 is cnp.uint32_t else np.uint64
    )

    if dim == 0:
        return

    # first row of v is all 1s
    for i in range(bits):
        v[0, i] = 1

    # Remaining rows of v (row 2 through dim, indexed by [1:dim])
    for d in range(1, dim):
        p = poly[d]
        m = bit_length(p) - 1

        # First m elements of row d comes from vinit
        v[d, :m] = vinit[d, :m]

        # Fill in remaining elements of v as in Section 2 (top of pg. 90) of:
        #
        # P. Bratley and B. L. Fox. Algorithm 659: Implementing sobol's
        # quasirandom sequence generator. ACM Trans.
        # Math. Softw., 14(1):88-100, Mar. 1988.
        #
        for j in range(m, bits):
            newv = v[d, j - m]
            pow2 = 1
            for k in range(m):
                pow2 = pow2 << 1
                if (p >> (m - 1 - k)) & 1:
                    newv = newv ^ (pow2 * v[d, j - k - 1])
            v[d, j] = newv

    # Multiply each column of v by power of 2:
    # v * [2^(bits-1), 2^(bits-2),..., 2, 1]
    pow2 = 1
    for d in range(bits):
        for i in range(dim):
            v[i, bits - 1 - d] *= pow2
        pow2 = pow2 << 1


def _draw(
    n,
    num_gen,
    const int dim,
    const cnp.float64_t scale,
    uint_32_64[:, ::1] sv,
    uint_32_64[::1] quasi,
    cnp.float64_t[:, ::1] sample
):
    # necessary wrapper to guide Cython for n, num_gen and scale
    cdef uint_32_64 n_ = n
    cdef uint_32_64 num_gen_ = num_gen
    draw(n_, num_gen_, dim, scale, sv, quasi, sample)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void draw(
    const uint_32_64 n,
    const uint_32_64 num_gen,
    const int dim,
    const cnp.float64_t scale,
    uint_32_64[:, ::1] sv,
    uint_32_64[::1] quasi,
    cnp.float64_t[:, ::1] sample
) nogil:
    cdef int j, l
    cdef uint_32_64 num_gen_loc = num_gen
    cdef uint_32_64 i, qtmp

    for i in range(n):
        l = low_0_bit(num_gen_loc)
        for j in range(dim):
            qtmp = quasi[j] ^ sv[j, l - 1]
            quasi[j] = qtmp
            sample[i, j] = qtmp * scale
        num_gen_loc += 1


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void _fast_forward(const uint_32_64 n,
                         const uint_32_64 num_gen,
                         const int dim,
                         uint_32_64[:, ::1] sv,
                         uint_32_64[::1] quasi) nogil:
    cdef int j, l
    cdef uint_32_64 num_gen_loc = num_gen
    cdef uint_32_64 i
    for i in range(n):
        l = low_0_bit(num_gen_loc)
        for j in range(dim):
            quasi[j] = quasi[j] ^ sv[j, l - 1]
        num_gen_loc += 1


@cython.boundscheck(False)
@cython.wraparound(False)
cdef uint_32_64 cdot_pow2(uint_32_64[::1] a) nogil:
    cdef int i
    cdef int size = a.shape[0]
    cdef uint_32_64 z = 0
    cdef uint_32_64 pow2 = 1
    for i in range(size):
        z += a[size - 1 - i] * pow2
        pow2 *= 2
    return z


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void _cscramble(const int dim,
                      const int bits,
                      uint_32_64[:, :, ::1] ltm,
                      uint_32_64[:, ::1] sv) nogil:
    """Scrambling using (left) linear matrix scramble (LMS)."""
    cdef int d, i, j, k, p
    cdef uint_32_64 l, lsmdp, t1, t2, vdj

    # Set diagonals of bits x bits arrays to 1
    for d in range(dim):
        for i in range(bits):
            ltm[d, i, i] = 1

    for d in range(dim):
        for j in range(bits):
            vdj = sv[d, j]
            l = 1
            t2 = 0
            for p in range(bits - 1, -1, -1):
                lsmdp = cdot_pow2(ltm[d, p, :])
                t1 = 0
                for k in range(bits):
                    t1 += ibits(lsmdp, k, 1) * ibits(vdj, k, 1)
                t1 = t1 % 2
                t2 = t2 + t1 * l
                l = 2 * l
            sv[d, j] = t2


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef void _fill_p_cumulative(cnp.float_t[::1] p,
                              cnp.float_t[::1] p_cumulative) nogil:
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
cpdef void _categorize(cnp.float_t[::1] draws,
                       cnp.float_t[::1] p_cumulative,
                       cnp.int_t[::1] result) nogil:
    cdef int i
    cdef int n_p = p_cumulative.shape[0]
    for i in range(draws.shape[0]):
        j = _find_index(p_cumulative, n_p, draws[i])
        result[j] = result[j] + 1


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int _find_index(cnp.float_t[::1] p_cumulative,
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
