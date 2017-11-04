from __future__ import absolute_import

cdef extern from "limits.h":
    unsigned long ULONG_MAX


def _comb_int(N, k):
    # Fast path with machine integers
    try:
        r = _comb_int_long(N, k)
        if r != 0:
            return r
    except (OverflowError, TypeError):
        pass

    # Fallback
    N = int(N)
    k = int(k)

    if k > N or N < 0 or k < 0:
        return 0

    M = N + 1
    nterms = min(k, N - k)

    numerator = 1
    denominator = 1
    for j in xrange(1, nterms + 1):
        numerator *= M - j
        denominator *= j

    return numerator // denominator


cdef unsigned long _comb_int_long(unsigned long N, unsigned long k):
    """
    Compute binom(N, k) for integers.
    Returns 0 if error/overflow encountered.
    """
    cdef unsigned long val, j, M, nterms

    if k > N or N == ULONG_MAX:
        return 0

    M = N + 1
    nterms = min(k, N - k)

    val = 1

    for j in xrange(1, nterms + 1):
        # Overflow check
        if val > ULONG_MAX // (M - j):
            return 0

        val *= M - j
        val //= j

    return val
