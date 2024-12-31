#pragma once

#include "xsf/cephes/exp10.h"
#include "xsf/cephes/exp2.h"

namespace xsf {

inline double expm1(double x) { return cephes::expm1(x); }

inline float expm1(float x) { return expm1(static_cast<double>(x)); }

// cexpm1(z) = cexp(z) - 1
//
// The imaginary part of this is easily computed via exp(z.real)*sin(z.imag)
// The real part is difficult to compute when there is cancellation e.g. when
// z.real = -log(cos(z.imag)).  There isn't a way around this problem  that
// doesn't involve computing exp(z.real) and/or cos(z.imag) to higher
// precision.
inline std::complex<double> expm1(std::complex<double> z) {
    if (!std::isfinite(std::real(z)) || !std::isfinite(std::imag(z))) {
        return std::exp(z) - 1.0;
    }

    double x;
    double ezr = 0;
    if (std::real(z) <= -40) {
        x = -1.0;
    } else {
        ezr = expm1(std::real(z));
        x = ezr * std::cos(std::imag(z)) + cosm1(std::imag(z));
    }

    // don't compute exp(zr) too, unless necessary
    double y;
    if (std::real(z) > -1.0) {
        y = (ezr + 1.0) * sin(std::imag(z));
    } else {
        y = exp(std::real(z)) * sin(std::imag(z));
    }

    return std::complex<double>{x, y};
}

inline std::complex<float> expm1(std::complex<float> z) {
    return static_cast<std::complex<float>>(expm1(static_cast<std::complex<double>>(z)));
}

double exp2(double x) { return cephes::exp2(x); }

float exp2(float x) { return exp2(static_cast<double>(x)); }

double exp10(double x) { return cephes::exp10(x); }

float exp10(float x) { return exp10(static_cast<double>(x)); }

} // namespace xsf
