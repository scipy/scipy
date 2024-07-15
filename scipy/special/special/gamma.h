#pragma once

#include "cephes/gamma.h"
#include "loggamma.h"

namespace special {

template <typename T>
SPECFUN_HOST_DEVICE T gamma(T x) {
    return cephes::Gamma(x);
}

template <typename T>
SPECFUN_HOST_DEVICE T gammaln(T x) {
    return cephes::lgam(x);
}

SPECFUN_HOST_DEVICE inline std::complex<double> gamma(std::complex<double> z) {
    // Compute Gamma(z) using loggamma.
    if (z.real() <= 0 && z == std::floor(z.real())) {
        // Poles
        set_error("gamma", SF_ERROR_SINGULAR, NULL);
        return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }
    return std::exp(loggamma(z));
}

SPECFUN_HOST_DEVICE inline std::complex<float> gamma(std::complex<float> z) {
    return static_cast<std::complex<float>>(gamma(static_cast<std::complex<double>>(z)));
}

} // namespace special
