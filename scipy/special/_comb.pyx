cdef extern from "limits.h":
    long LONG_MAX

# Compute N choose k in integer arithmetic.
def _comb(N, k):
    cdef:
        object numerator, denominator
        object j, M, nterms
    
    if (k > N) or (N < 0) or (k < 0):
        return 0

    M = N + 1
    nterms = min(k, N - k)
    
    if nterms == 0:
        return 1

    if M < pow(LONG_MAX, 1./nterms):
        # delegate to a fully typed version unless overflow
        # The condition is rather loose: numerator < M**nterms
        return _comb_int_unsafe(N, k)
    elif M < LONG_MAX:
        # both N and k fit into a C long, (but
        # intermediaries would overflow)
        return _comb_interm_obj(N, k)

    numerator = 1
    denominator = 1
    for j in xrange(1, nterms + 1):
        numerator *= M - j
        denominator *= j

    return numerator // denominator


# N and k fit into a C long, use object intermediaries
cdef object _comb_interm_obj(long N, long k):
    cdef:
        object numerator, denominator
        long j, M, nterms
    
    if (k > N) or (N < 0) or (k < 0):
        return 0

    M = N + 1
    nterms = min(k, N - k)
    
    numerator = 1
    denominator = 1
    for j in xrange(1, nterms + 1):
        numerator *= M - j
        denominator *= j
        
    return numerator // denominator


# not only N and k but also intermediaries fit into a C long
cdef long _comb_int_unsafe(long N, long k):
    cdef:
        long numerator, denominator
        long j, M, nterms
    
    if (k > N) or (N < 0) or (k < 0):
        return 0

    M = N + 1
    nterms = min(k, N - k)  
    
    numerator = 1
    denominator = 1
    for j in xrange(1, nterms + 1):
        numerator *= M - j
        denominator *= j
        
    return numerator // denominator

