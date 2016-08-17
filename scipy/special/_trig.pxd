# Implement sin(pi*z) and cos(pi*z) for complex z. Since the periods
# of these functions are integral (and thus better representable in
# floating point), it's possible to compute them with greater accuracy
# than sin(z), cos(z).

from libc.math cimport ceil, floor, M_PI
from _complexstuff cimport number_t, zreal, zsin, zcos, zabs

DEF tol = 2.220446049250313e-16


cdef inline number_t sinpi(number_t z) nogil:
    """Compute sin(pi*z) by shifting z to (-0.5, 0.5]."""
    cdef:
        double p = ceil(zreal(z))
        double hp = p/2

    # Make p the even integer closest to z
    if hp != ceil(hp):
        p -= 1
    # zn.real is in (-1, 1]
    z -= p
    # Reflect zn.real in (0.5, 1] to [0, 0.5).
    if zreal(z) > 0.5:
        z = 1 - z
    # Reflect zn.real in (-1, -0.5) to (-0.5, 0)
    if zreal(z) < -0.5:
        z = -1 - z
    return zsin(M_PI*z)


cdef inline number_t cospi(number_t z) nogil:
    """Compute cos(pi*z) by shifting z to (-1, 1]."""
    cdef:
        double p = ceil(zreal(z))
        double hp = p/2

    # Make p the even integer closest to z
    if hp != ceil(hp):
        p -= 1
    # zn.real is in (-1, 1].
    z -= p
    if zabs(z - 0.5) < 0.2:
        return cospi_taylor(z)
    elif zabs(z + 0.5) < 0.2:
        return cospi_taylor(-z)
    else:
        return zcos(M_PI*z)


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
    for n in range(1, 20):
        term *= -zz/((2*n + 1)*(2*n))
        res += term
        if zabs(term) <= tol*zabs(res):
            break
    return res
