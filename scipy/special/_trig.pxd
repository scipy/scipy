# Implement sin(pi*z) and cos(pi*z) for complex z. Since the periods
# of these functions are integral (and thus better representable in
# floating point), it's possible to compute them with greater accuracy
# than sin(z), cos(z).
#
from libc.math cimport sin, cos, sinh, cosh, exp, fabs, fmod, M_PI

from ._cephes cimport sinpi as dsinpi, cospi as dcospi
from ._complexstuff cimport number_t, double_complex, zpack

cdef extern from "numpy/npy_math.h":
    double npy_copysign(double x, double y) nogil
    double NPY_INFINITY


cdef inline double complex csinpi(double complex z) nogil:
    """Compute sin(pi*z) for complex arguments."""
    cdef:
        double x = z.real
        double piy = M_PI*z.imag
        double abspiy = fabs(piy)
        double sinpix = sinpi(x)
        double cospix = cospi(x)
        double exphpiy, coshfac, sinhfac

    if abspiy < 700:
        return zpack(sinpix*cosh(piy), cospix*sinh(piy))

    # Have to be careful--sinh/cosh could overflow while cos/sin are
    # small. At this large of values
    #
    # cosh(y) ~ exp(y)/2
    # sinh(y) ~ sgn(y)*exp(y)/2
    #
    # so we can compute exp(y/2), scale by the right factor of sin/cos
    # and then multiply by exp(y/2) to avoid overflow.
    exphpiy = exp(abspiy/2)
    if exphpiy == NPY_INFINITY:
        if sinpix == 0:
            # Preserve the sign of zero
            coshfac = npy_copysign(0.0, sinpix)
        else:
            coshfac = npy_copysign(NPY_INFINITY, sinpix)
        if cospix == 0:
            sinhfac = npy_copysign(0.0, cospix)
        else:
            sinhfac = npy_copysign(NPY_INFINITY, cospix)
        return zpack(coshfac, sinhfac)

    coshfac = 0.5*sinpix*exphpiy
    sinhfac = 0.5*cospix*exphpiy
    return zpack(coshfac*exphpiy, sinhfac*exphpiy)


cdef inline double complex ccospi(double complex z) nogil:
    """Compute cos(pi*z) for complex arguments."""
    cdef:
        double x = z.real
        double piy = M_PI*z.imag
        double abspiy = fabs(piy)
        double sinpix = sinpi(x)
        double cospix = cospi(x)
        double exphpiy, coshfac, sinhfac

    if abspiy < 700:
        return zpack(cospix*cosh(piy), -sinpix*sinh(piy))

    # See csinpi(z) for an idea of what's going on here
    exphpiy = exp(abspiy/2)
    if exphpiy == NPY_INFINITY:
        if sinpix == 0:
            coshfac = npy_copysign(0.0, cospix)
        else:
            coshfac = npy_copysign(NPY_INFINITY, cospix)
        if cospix == 0:
            sinhfac = npy_copysign(0.0, sinpix)
        else:
            sinhfac = npy_copysign(NPY_INFINITY, sinpix)
        return zpack(coshfac, sinhfac)

    coshfac = 0.5*cospix*exphpiy
    sinhfac = 0.5*sinpix*exphpiy
    return zpack(coshfac*exphpiy, sinhfac*exphpiy)


cdef inline number_t sinpi(number_t z) nogil:
    if number_t is double:
        return dsinpi(z)
    else:
        return csinpi(z)


cdef inline number_t cospi(number_t z) nogil:
    if number_t is double:
        return dcospi(z)
    else:
        return ccospi(z)
