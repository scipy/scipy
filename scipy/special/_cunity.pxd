cimport numpy as np
from libc.math cimport fabs, sin, cos, exp, atan2

from ._complexstuff cimport (
    zisfinite, zabs, zpack, npy_cdouble_from_double_complex,
    double_complex_from_npy_cdouble)

cdef extern from "_complexstuff.h":
    np.npy_cdouble npy_clog(np.npy_cdouble x) nogil
    np.npy_cdouble npy_cexp(np.npy_cdouble x) nogil


cdef extern from "dd_real_wrappers.h":
    ctypedef struct double2:
        double hi
        double lo

    double2 dd_create_d(double x) nogil
    double2 dd_add(const double2* a, const double2* b) nogil
    double2 dd_mul(const double2* a, const double2* b) nogil
    double dd_to_double(const double2* a) nogil

cdef extern from "special_wrappers.h" nogil:
    double cephes_cosm1_wrap(double x)
    double cephes_expm1_wrap(double x)
    double cephes_log1p_wrap(double x)

# log(z + 1) = log(x + 1 + 1j*y)
#             = log(sqrt((x+1)**2 + y**2)) + 1j*atan2(y, x+1)
#
# Using atan2(y, x+1) for the imaginary part is always okay.  The real part
# needs to be calculated more carefully.  For |z| large, the naive formula
# log(z + 1) can be used.  When |z| is small, rewrite as
#
# log(sqrt((x+1)**2 + y**2)) = 0.5*log(x**2 + 2*x +1 + y**2)
#       = 0.5 * log1p(x**2 + y**2 + 2*x)
#       = 0.5 * log1p(hypot(x,y) * (hypot(x, y) + 2*x/hypot(x,y)))
#
# This expression suffers from cancellation when x < 0 and
# y = +/-sqrt(2*fabs(x)). To get around this cancellation problem, we use
# double-double precision when necessary.
cdef inline double complex clog1p(double complex z) noexcept nogil:
    cdef double zr, zi, x, y, az, azi
    cdef np.npy_cdouble ret

    if not zisfinite(z):
        z = z + 1
        ret = npy_clog(npy_cdouble_from_double_complex(z))
        return double_complex_from_npy_cdouble(ret)

    zr = z.real
    zi = z.imag

    if zi == 0.0 and zr >= -1.0:
        return zpack(cephes_log1p_wrap(zr), 0.0)

    az = zabs(z)
    if az < 0.707:
        azi = fabs(zi)
        if zr < 0 and fabs(-zr - azi*azi/2)/(-zr) < 0.5:
            return clog1p_ddouble(zr, zi)
        else:
            x = 0.5 * cephes_log1p_wrap(az*(az + 2*zr/az))
            y = atan2(zi, zr + 1.0)
            return zpack(x, y)

    z = z + 1.0
    ret = npy_clog(npy_cdouble_from_double_complex(z))
    return double_complex_from_npy_cdouble(ret)

cdef inline double complex clog1p_ddouble(double zr, double zi) noexcept nogil:
    cdef double x, y
    cdef double2 r, i, two, rsqr, isqr, rtwo, absm1

    r = dd_create_d(zr)
    i = dd_create_d(zi)
    two = dd_create_d(2.0)

    rsqr = dd_mul(&r,& r)
    isqr = dd_mul(&i, &i)
    rtwo = dd_mul(&two, &r)
    absm1 = dd_add(&rsqr, &isqr)
    absm1 = dd_add(&absm1, &rtwo)

    x = 0.5 * cephes_log1p_wrap(dd_to_double(&absm1))
    y = atan2(zi, zr+1.0)
    return zpack(x, y)

# cexpm1(z) = cexp(z) - 1
#
# The imaginary part of this is easily computed via exp(z.real)*sin(z.imag)
# The real part is difficult to compute when there is cancellation e.g. when
# z.real = -log(cos(z.imag)).  There isn't a way around this problem  that
# doesn't involve computing exp(z.real) and/or cos(z.imag) to higher
# precision.
cdef inline double complex cexpm1(double complex z) noexcept nogil:
    cdef double zr, zi, ezr, x, y
    cdef np.npy_cdouble ret

    if not zisfinite(z):
        ret = npy_cexp(npy_cdouble_from_double_complex(z))
        return double_complex_from_npy_cdouble(ret) - 1.0

    zr = z.real
    zi = z.imag

    if zr <= -40:
        x = -1.0
    else:
        ezr = cephes_expm1_wrap(zr)
        x = ezr*cos(zi) + cephes_cosm1_wrap(zi)
    # don't compute exp(zr) too, unless necessary
    if zr > -1.0:
        y = (ezr + 1.0)*sin(zi)
    else:
        y = exp(zr)*sin(zi)

    return zpack(x, y)
