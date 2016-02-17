cimport numpy as np
from libc.math cimport fabs, sin, cos, exp
from _complexstuff cimport zisfinite, zabs, zpack

cdef extern from "_complexstuff.h":
    double npy_atan2(double y, double x) nogil
    np.npy_cdouble npy_clog(np.npy_cdouble x) nogil
    np.npy_cdouble npy_cexp(np.npy_cdouble x) nogil
    

cdef extern from "c_misc/double2.h":
    ctypedef struct double2_t:
        double x[2]

    void double2_init(double2_t* a, double y) nogil
    void double2_add(double2_t* a, double2_t* b, double2_t* c) nogil
    void double2_mul(double2_t* a, double2_t* b, double2_t* c) nogil
    double double2_double(double2_t* a) nogil

cdef extern from "cephes.h":
    double log1p(double x) nogil
    double expm1(double x) nogil
    double cosm1(double x) nogil

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
# This expression suffers from cancelation when x < 0 and
# y = +/-sqrt(2*fabs(x)). To get around this cancelation problem, we use
# double-double precision when necessary.
cdef inline double complex clog1p(double complex z) nogil:
    cdef double zr, zi, x, y, az, azi
    cdef np.npy_cdouble ret

    if not zisfinite(z):
        z = z + 1 
        ret = npy_clog((<np.npy_cdouble*>&z)[0])
        return (<double complex*>&ret)[0]

    zr = z.real
    zi = z.imag

    if zi == 0.0 and zr >= -1.0:
        return zpack(log1p(zr), 0.0) 

    az = zabs(z)
    if az < 0.707:
        azi = fabs(zi)
        if zr < 0 and fabs(-zr - azi*azi/2)/(-zr) < 0.5:
            return clog1p_ddouble(zr, zi)
        else:
            x = 0.5 * log1p(az*(az + 2*zr/az))
            y = npy_atan2(zi, zr + 1.0)
            return zpack(x, y)

    z = z + 1.0
    ret = npy_clog((<np.npy_cdouble*>&z)[0])
    return (<double complex*>&ret)[0]

cdef inline double complex clog1p_ddouble(double zr, double zi) nogil:
    cdef double x, y
    cdef double2_t r, i, two, rsqr, isqr, rtwo, absm1

    double2_init(&r, zr)
    double2_init(&i, zi)
    double2_init(&two, 2.0)
    
    double2_mul(&r, &r, &rsqr)
    double2_mul(&i, &i, &isqr)
    double2_mul(&two, &r, &rtwo)
    double2_add(&rsqr, &isqr, &absm1)
    double2_add(&absm1, &rtwo, &absm1)

    x = 0.5 * log1p(double2_double(&absm1))
    y = npy_atan2(zi, zr+1.0)
    return zpack(x, y) 
 
# cexpm1(z) = cexp(z) - 1
# 
# The imaginary part of this is easily computed via exp(z.real)*sin(z.imag)
# The real part is difficult to compute when there is cancelation e.g. when
# z.real = -log(cos(z.imag)).  There isn't a way around this problem  that
# doesn't involve computing exp(z.real) and/or cos(z.imag) to higher
# precision.
cdef inline double complex cexpm1(double complex z) nogil:
    cdef double zr, zi, ezr, x, y
    cdef np.npy_cdouble ret

    if not zisfinite(z):
        ret = npy_cexp((<np.npy_cdouble*>&z)[0])
        return (<double complex*>&ret)[0] - 1.0

    zr = z.real
    zi = z.imag

    if zr <= -40:
        x = -1.0
    else:
        ezr = expm1(zr)
        x = ezr*cos(zi) + cosm1(zi)
    # don't compute exp(zr) too, unless necessary
    if zr > -1.0:
        y = (ezr + 1.0)*sin(zi)
    else:
        y = exp(zr)*sin(zi)

    return zpack(x, y)

