# Implement sin(pi*z) and cos(pi*z) for complex z. Since the periods
# of these functions are integral (and thus more representable in
# floating point), it's possible to compute them with greater accuracy
# than sin(z), cos(z).

from libc.math cimport ceil, floor, M_PI
from _complexstuff cimport number_t, zreal, zsin, zcos, zabs

DEF tol = 2.220446049250313e-16


cdef inline number_t sinpi(number_t z) nogil:
    """
    Compute sin(pi*z) by finding zn so that Re(zn) is in [0, 0.5]
    and sin(pi*z) = sin(pi*zn) or sin(pi*z) = -sin(pi*zn).

    """
    cdef:
        double p = ceil(zreal(z))
        double hp = p/2
        number_t zn

    if z == p:
        return 0
    elif zabs(z) < 0.5:
        # If z is small shifting costs us accuracy
        return zsin(M_PI*z)

    # Make p the even integer to the left of z
    if hp != ceil(hp):
        p -= 1
    # zn.real is in [-1, 1)
    zn = z - p
    # Reflect zn.real in (0.5, 1) to (0, 0.5).
    if zreal(zn) > 0.5:
        zn = 1 - zn
    # Reflect zn.real in [-1, -0.5) to (-0.5, 0]
    if zreal(zn) < -0.5:
        zn = -1 - zn
    return zsin(M_PI*zn)


cdef inline number_t cospi(number_t z) nogil:
    """
    Compute cos(pi*z) by finding zn so that Re(zn) is in [0, 0.5] and
    cos(pi*z) = cos(pi*zn) or cos(pi*z) = -cos(pi*zn).

    """
    cdef:
        int sgn = 1
        double p = floor(zreal(z))
        double hp = p/2
        number_t zn

    # Make p the even integer to the left of z
    if hp != floor(hp):
        p -= 1
    # zn.real is in [0, 2).
    zn = z - p
    # Shift z.real to [0, 1], negate.
    if zreal(zn) > 1:
        zn -= 1
        sgn *= -1
    if zabs(zn - 0.5) < 0.1:
        return sgn*cospi_taylor(zn)
    else:
        return sgn*zcos(M_PI*zn)


cdef inline number_t cospi_taylor(number_t z) nogil:
    """
    Taylor series for cos(pi*z) around z = 0.5. Since the root is
    exactly representable in double precision we get gains over
    just using cos(z) here.

    """
    cdef:
        int n
        number_t zz, term, res

    z = M_PI*(z - 0.5)
    zz = z*z
    term = -z
    res = term
    for n in range(1, 16):
        term *= -zz/((2*n + 1)*(2*n))
        res += term
        if zabs(term) <= tol*zabs(res):
            break
    return res
