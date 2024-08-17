/* Translated from Cython into C++ by SciPy developers in 2023.
 *
 * Original author: Josh Wilson, 2016.
 * Original author: Takuma Yoshimura, 2024.
 */

/* Implement sin(pi*z), cos(pi*z), tan(pi*z) and cot(pi*z) for complex z.
 * Since the periods of these functions are integral (and thus better
 * representable in floating point), it's possible to compute them with
 * greater accuracy than sin(z), cos(z), tan(z), cot(z).
 */

#pragma once

#include "cephes/trig.h"
#include "config.h"
#include "evalpoly.h"

namespace xsf {

template <typename T>
XSF_HOST_DEVICE T sinpi(T x) {
    return cephes::sinpi(x);
}

template <typename T>
XSF_HOST_DEVICE std::complex<T> sinpi(std::complex<T> z) {
    T x = z.real();
    T piy = M_PI * z.imag();
    T abspiy = std::abs(piy);
    T sinpix = cephes::sinpi(x);
    T cospix = cephes::cospi(x);

    if (abspiy < 700) {
        return {sinpix * std::cosh(piy), cospix * std::sinh(piy)};
    }

    /* Have to be careful--sinh/cosh could overflow while cos/sin are small.
     * At this large of values
     *
     * cosh(y) ~ exp(y)/2
     * sinh(y) ~ sgn(y)*exp(y)/2
     *
     * so we can compute exp(y/2), scale by the right factor of sin/cos
     * and then multiply by exp(y/2) to avoid overflow. */
    T exphpiy = std::exp(abspiy / 2);
    T coshfac;
    T sinhfac;
    if (exphpiy == std::numeric_limits<T>::infinity()) {
        if (sinpix == 0.0) {
            // Preserve the sign of zero.
            coshfac = std::copysign(0.0, sinpix);
        } else {
            coshfac = std::copysign(std::numeric_limits<T>::infinity(), sinpix);
        }
        if (cospix == 0.0) {
            // Preserve the sign of zero.
            sinhfac = std::copysign(0.0, cospix);
        } else {
            sinhfac = std::copysign(std::numeric_limits<T>::infinity(), cospix);
        }
        return {coshfac, sinhfac};
    }

    coshfac = 0.5 * sinpix * exphpiy;
    sinhfac = 0.5 * cospix * exphpiy;
    return {coshfac * exphpiy, sinhfac * exphpiy};
}

template <typename T>
XSF_HOST_DEVICE T cospi(T x) {
    return cephes::cospi(x);
}

template <typename T>
XSF_HOST_DEVICE std::complex<T> cospi(std::complex<T> z) {
    T x = z.real();
    T piy = M_PI * z.imag();
    T abspiy = std::abs(piy);
    T sinpix = cephes::sinpi(x);
    T cospix = cephes::cospi(x);

    if (abspiy < 700) {
        return {cospix * std::cosh(piy), -sinpix * std::sinh(piy)};
    }

    // See csinpi(z) for an idea of what's going on here.
    T exphpiy = std::exp(abspiy / 2);
    T coshfac;
    T sinhfac;
    if (exphpiy == std::numeric_limits<T>::infinity()) {
        if (sinpix == 0.0) {
            // Preserve the sign of zero.
            coshfac = std::copysign(0.0, cospix);
        } else {
            coshfac = std::copysign(std::numeric_limits<T>::infinity(), cospix);
        }
        if (cospix == 0.0) {
            // Preserve the sign of zero.
            sinhfac = std::copysign(0.0, sinpix);
        } else {
            sinhfac = std::copysign(std::numeric_limits<T>::infinity(), sinpix);
        }
        return {coshfac, sinhfac};
    }

    coshfac = 0.5 * cospix * exphpiy;
    sinhfac = 0.5 * sinpix * exphpiy;
    return {coshfac * exphpiy, sinhfac * exphpiy};
}

template <typename T>
XSF_HOST_DEVICE T tanpi(T x) {
    return cephes::tanpi(x);
}

template <typename T>
XSF_HOST_DEVICE std::complex<T> tanpi(std::complex<T> z) {
    T u = std::exp(-std::abs(2.0 * z.imag() * M_PI));

    if (u == 1.0) {
        return std::complex<T>(tanpi(z.real()), 0.0);
    }

    T n = 1.0 + u * (2.0 * cospi(2.0 * z.real()) + u);

    T r = 2.0 * u * sinpi(2.0 * z.real()) / n;
    T i = (u + 1.0) * (u - 1.0) / n;

    std::complex<T> y = (z.imag() > 0.0)
        ? std::complex<T>(r, -i) 
        : std::complex<T>(r, i);

    return y;
}

template <typename T>
XSF_HOST_DEVICE T cotpi(T x) {
    return cephes::cotpi(x);
}

template <typename T>
XSF_HOST_DEVICE std::complex<T> cotpi(std::complex<T> z) {
    T u = std::exp(-std::abs(2.0 * z.imag() * M_PI));
    
    if (u == 1.0) {
        return std::complex<T>(cotpi(z.real()), 0.0);
    }

    T n = 1.0 + u * (-2.0 * cospi(2.0 * z.real()) + u);

    T r = 2.0 * u * sinpi(2.0 * z.real()) / n;
    T i = (u + 1.0) * (u - 1.0) / n;

    std::complex<T> y = (z.imag() > 0.0)
        ? std::complex<T>(r, i) 
        : std::complex<T>(r, -i);

    return y;
}

} // namespace xsf
