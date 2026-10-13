# Exact Wasserstein-1 integration for validated binary64 inputs.
# Included by _stats.pyx; no additional extension or external dependency.
from libc.math cimport frexp, ldexp, fabs, isfinite
from libc.stdint cimport int64_t, uint32_t, uint64_t
cimport numpy as cnp
from math import gcd as _w1_gcd


cdef double _w1_round_ratio(object numerator, object denominator, int exponent) except *:
    """Exact nearest-even rounding of a nonnegative rational times 2**exponent."""
    cdef int log, quantum, shift
    cdef object quotient, remainder, divisor
    if not numerator:
        return 0.0
    log = numerator.bit_length() - denominator.bit_length()
    if log >= 0:
        log -= numerator < (denominator << log)
    else:
        log -= (numerator << -log) < denominator
    quantum = max(log + exponent - 52, -1074)
    shift = exponent - quantum
    if shift >= 0:
        quotient, remainder = divmod(numerator << shift, denominator)
        divisor = denominator
    else:
        divisor = denominator << -shift
        quotient, remainder = divmod(numerator, divisor)
    if 2 * remainder > divisor or (2 * remainder == divisor and quotient % 2):
        quotient += 1
    # quotient has at most 54 bits and is exactly representable: if it has
    # 54 bits it is exactly 2**53. ldexp applies an exact power-of-two scale,
    # including subnormal outputs; overflow yields positive infinity.
    return ldexp(float(quotient), quantum)


cdef int _w1_minimum_exponent(const double[:] values):
    cdef int minimum = 1024, exponent
    cdef Py_ssize_t i
    for i in range(values.shape[0]):
        if values[i] != 0:
            frexp(values[i], &exponent)
            if exponent - 53 < minimum:
                minimum = exponent - 53
    return minimum if minimum != 1024 else 0


cdef object _w1_big_integer(double value, int minimum):
    cdef int exponent
    cdef int64_t mantissa
    if value == 0:
        return 0
    mantissa = <int64_t>ldexp(frexp(value, &exponent), 53)
    return (<object>mantissa) << (exponent - 53 - minimum)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef object _w1_exact_distance(const double[:] u, const double[:] v, uw=None, vw=None):
    cdef Py_ssize_t n = u.shape[0], m = v.shape[0], i, j
    cdef const double[:] wu, wv
    cdef const double[:] values = np.concatenate((u, v))
    cdef const cnp.intp_t[:] order = np.argsort(values, kind='stable')
    cdef int eu, ev, ex = _w1_minimum_exponent(values)
    cdef list weights = []
    cdef object totalu = 0, totalv = 0, weight, factoru, factorv, common
    cdef object imbalance = 0, cost = 0, previous = 0, position, change, updated
    if uw is None:
        weights = [1] * n
        totalu = n
    else:
        wu = uw
        eu = _w1_minimum_exponent(wu)
        for i in range(n):
            weight = _w1_big_integer(wu[i], eu)
            totalu += weight
            weights.append(weight)
    if vw is None:
        weights.extend([1] * m)
        totalv = m
    else:
        wv = vw
        ev = _w1_minimum_exponent(wv)
        for i in range(m):
            weight = _w1_big_integer(wv[i], ev)
            totalv += weight
            weights.append(weight)
    common = _w1_gcd(totalu, totalv)
    factoru, factorv = totalv // common, totalu // common
    for i in range(n + m):
        j = order[i]
        change = weights[j] * (factoru if j < n else -factorv)
        position = _w1_big_integer(values[j], ex)
        if imbalance:
            cost += abs(imbalance) * (position - previous)
        imbalance += change
        previous = position
    return _w1_round_ratio(cost, totalu * factoru, ex)


cdef struct _w1_UInt:
    uint32_t w[8]
    int used


cdef void _w1_zero(_w1_UInt* a) noexcept nogil:
    cdef int i
    for i in range(8):
        a.w[i] = 0
    a.used = 0


cdef void _w1_trim(_w1_UInt* a) noexcept nogil:
    while a.used and not a.w[a.used - 1]:
        a.used -= 1


cdef void _w1_add(_w1_UInt* a, const _w1_UInt* b) noexcept nogil:
    cdef int i, n = max(a.used, b.used)
    cdef uint64_t t, carry = 0
    for i in range(n):
        t = <uint64_t>a.w[i] + b.w[i] + carry
        a.w[i] = <uint32_t>t
        carry = t >> 32
    a.used = n
    if carry and n < 8:
        a.w[n] = <uint32_t>carry
        a.used += 1


cdef int _w1_compare(const _w1_UInt* a, const _w1_UInt* b) noexcept nogil:
    cdef int i
    if a.used != b.used:
        return 1 if a.used > b.used else -1
    for i in range(a.used - 1, -1, -1):
        if a.w[i] != b.w[i]:
            return 1 if a.w[i] > b.w[i] else -1
    return 0


cdef void _w1_subtract(_w1_UInt* a, const _w1_UInt* b) noexcept nogil:
    # Requires a >= b.
    cdef int i
    cdef uint64_t t, borrow = 0
    for i in range(a.used):
        t = <uint64_t>b.w[i] + borrow
        borrow = <uint64_t>a.w[i] < t
        a.w[i] = <uint32_t>(<uint64_t>a.w[i] - t)
    _w1_trim(a)


cdef void _w1_multiply(const _w1_UInt* a, const _w1_UInt* b, _w1_UInt* out) noexcept nogil:
    cdef int i, j
    cdef uint64_t t, carry
    _w1_zero(out)
    if not a.used or not b.used:
        return
    for i in range(a.used):
        carry = 0
        for j in range(b.used):
            t = <uint64_t>a.w[i] * b.w[j] + out.w[i + j] + carry
            out.w[i + j] = <uint32_t>t
            carry = t >> 32
        if i + b.used < 8:
            out.w[i + b.used] = <uint32_t>carry
    out.used = min(8, a.used + b.used)
    _w1_trim(out)


cdef int _w1_trailing_zeros(uint64_t v) noexcept nogil:
    cdef int count = 0
    if not (v & <uint64_t>0xffffffff):
        count += 32
        v >>= 32
    if not (v & 0xffff):
        count += 16
        v >>= 16
    if not (v & 0xff):
        count += 8
        v >>= 8
    if not (v & 0xf):
        count += 4
        v >>= 4
    if not (v & 3):
        count += 2
        v >>= 2
    if not (v & 1):
        count += 1
    return count


@cython.boundscheck(False)
@cython.wraparound(False)
cdef bint _w1_range_info(const double[:] values, int* quantum, int* bits) noexcept nogil:
    cdef Py_ssize_t i
    cdef int e, low, smallest = 1024, largest = -1074
    cdef uint64_t m
    for i in range(values.shape[0]):
        if not isfinite(values[i]):
            return False
        if values[i] != 0:
            m = <uint64_t>ldexp(frexp(fabs(values[i]), &e), 53)
            low = e - 53 + _w1_trailing_zeros(m)
            smallest = min(smallest, low)
            largest = max(largest, e)
    if smallest == 1024:
        quantum[0] = 0
        bits[0] = 0
    else:
        quantum[0] = smallest
        bits[0] = largest - smallest
    return True


cdef void _w1_integer_value(double value, int quantum, _w1_UInt* out) noexcept nogil:
    cdef int e, shift, offset, residual, i
    cdef uint64_t mantissa, part
    _w1_zero(out)
    if value == 0:
        return
    mantissa = <uint64_t>ldexp(frexp(fabs(value), &e), 53)
    shift = e - 53 - quantum
    if shift < 0:
        mantissa >>= -shift
        shift = 0
    offset, residual = shift // 32, shift % 32
    for i in range(2):
        part = (mantissa & <uint64_t>0xffffffff) << residual
        out.w[offset + i] |= <uint32_t>part
        if part >> 32:
            out.w[offset + i + 1] |= <uint32_t>(part >> 32)
        mantissa >>= 32
        if not mantissa:
            break
    out.used = min(8, offset + i + 2)
    _w1_trim(out)


cdef object _w1_python_integer(const _w1_UInt* a):
    cdef object value = 0
    cdef int i
    for i in range(a.used - 1, -1, -1):
        value = (value << 32) + a.w[i]
    return value


@cython.boundscheck(False)
@cython.wraparound(False)
def _wasserstein_distance_finite(const double[:] u, const double[:] v, uw=None, vw=None):
    cdef Py_ssize_t n = u.shape[0], m = v.shape[0], i, j
    cdef const double[:] values = np.concatenate((u, v))
    cdef const double[:] wu, wv, sortedu, sortedv, active, active_weights
    cdef const cnp.intp_t[:] order
    cdef int eu = 0, ev = 0, ex, bu, bv, bx, sign = 0, direction
    cdef double x, previous_x = 0, target
    cdef bint same, uniform
    cdef int active_quantum
    cdef _w1_UInt tu, tv, weight, change, imbalance, position, previous, delta, term, cost, denom
    cdef object numerator, denominator
    if not _w1_range_info(values, &ex, &bx):
        return None
    if n == m:
        same = True
        for i in range(n):
            if u[i] != v[i]:
                same = False
                break
        if same:
            if uw is None and vw is None:
                return 0.0
            if uw is not None and vw is not None:
                wu, wv = uw, vw
                for i in range(n):
                    if wu[i] != wv[i]:
                        same = False
                        break
                if same:
                    return 0.0
    if uw is None:
        bu = int(n).bit_length()
    else:
        wu = uw
        _w1_range_info(wu, &eu, &bu)
        bu += int(n - 1).bit_length()
    if vw is None:
        bv = int(m).bit_length()
    else:
        wv = vw
        _w1_range_info(wv, &ev, &bv)
        bv += int(m - 1).bit_length()
    # Total mass product bounds every absolute prefix imbalance. The sum of
    # widths is at most twice the largest absolute coordinate. Thus this
    # bound also protects every multiplication and the complete integral.
    if bu + bv + bx + 1 > 256:
        return _w1_exact_distance(u, v, uw, vw)
    if n == 1 or m == 1:
        if m == 1:
            active, target, active_quantum = u, v[0], eu
            uniform = uw is None
            if not uniform:
                active_weights = wu
        else:
            active, target, active_quantum = v, u[0], ev
            uniform = vw is None
            if not uniform:
                active_weights = wv
        _w1_zero(&cost)
        _w1_zero(&tu)
        _w1_integer_value(target, ex, &previous)
        with nogil:
            for i in range(active.shape[0]):
                _w1_integer_value(1.0 if uniform else active_weights[i],
                              active_quantum, &weight)
                _w1_add(&tu, &weight)
                _w1_integer_value(active[i], ex, &position)
                if (active[i] < 0) != (target < 0):
                    _w1_add(&position, &previous)
                elif _w1_compare(&position, &previous) >= 0:
                    _w1_subtract(&position, &previous)
                else:
                    delta = previous
                    _w1_subtract(&delta, &position)
                    position = delta
                _w1_multiply(&weight, &position, &term)
                _w1_add(&cost, &term)
        numerator, denominator = _w1_python_integer(&cost), _w1_python_integer(&tu)
        return _w1_round_ratio(numerator, denominator, ex)
    if uw is None and vw is None and n == m:
        sortedu = np.sort(np.asarray(u))
        sortedv = np.sort(np.asarray(v))
        _w1_zero(&cost)
        with nogil:
            for i in range(n):
                _w1_integer_value(sortedu[i], ex, &position)
                _w1_integer_value(sortedv[i], ex, &previous)
                if (sortedu[i] < 0) != (sortedv[i] < 0):
                    _w1_add(&position, &previous)
                elif _w1_compare(&position, &previous) >= 0:
                    _w1_subtract(&position, &previous)
                else:
                    _w1_subtract(&previous, &position)
                    position = previous
                _w1_add(&cost, &position)
        numerator, denominator = _w1_python_integer(&cost), int(n)
        return _w1_round_ratio(numerator, denominator, ex)
    order = np.argsort(values)
    _w1_zero(&tu)
    _w1_zero(&tv)
    _w1_zero(&imbalance)
    _w1_zero(&previous)
    _w1_zero(&cost)
    with nogil:
        for i in range(n):
            _w1_integer_value(1.0 if uw is None else wu[i], eu, &weight)
            _w1_add(&tu, &weight)
        for i in range(m):
            _w1_integer_value(1.0 if vw is None else wv[i], ev, &weight)
            _w1_add(&tv, &weight)
        for i in range(n + m):
            j = order[i]
            x = values[j]
            _w1_integer_value(x, ex, &position)
            if imbalance.used:
                if previous_x >= 0:
                    delta = position
                    _w1_subtract(&delta, &previous)
                elif x <= 0:
                    delta = previous
                    _w1_subtract(&delta, &position)
                else:
                    delta = previous
                    _w1_add(&delta, &position)
                _w1_multiply(&imbalance, &delta, &term)
                _w1_add(&cost, &term)
            if j < n:
                _w1_integer_value(1.0 if uw is None else wu[j], eu, &weight)
                _w1_multiply(&weight, &tv, &change)
                direction = 1
            else:
                _w1_integer_value(1.0 if vw is None else wv[j - n], ev, &weight)
                _w1_multiply(&weight, &tu, &change)
                direction = -1
            if sign == direction or not imbalance.used:
                _w1_add(&imbalance, &change)
                sign = direction
            elif _w1_compare(&imbalance, &change) >= 0:
                _w1_subtract(&imbalance, &change)
            else:
                _w1_subtract(&change, &imbalance)
                imbalance = change
                sign = direction
            previous = position
            previous_x = x
        _w1_multiply(&tu, &tv, &denom)
    numerator, denominator = _w1_python_integer(&cost), _w1_python_integer(&denom)
    return _w1_round_ratio(numerator, denominator, ex)
