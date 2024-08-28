#pragma once

#include "cephes/gamma.h"
#include "loggamma.h"

namespace xsf {

template <typename T>
XSF_HOST_DEVICE T gamma(T x) {
    return cephes::Gamma(x);
}

XSF_HOST_DEVICE inline double gammaln(double x) { return cephes::lgam(x); }

XSF_HOST_DEVICE inline float gammaln(float x) { return gammaln(static_cast<double>(x)); }

XSF_HOST_DEVICE inline double gammasgn(double x) { return cephes::gammasgn(x); }

XSF_HOST_DEVICE inline float gammasgn(float x) { return gammasgn(static_cast<double>(x)); }

XSF_HOST_DEVICE inline std::complex<double> gamma(std::complex<double> z) {
    // Compute Gamma(z) using loggamma.
    if (z.real() <= 0 && z == std::floor(z.real())) {
        // Poles
        set_error("gamma", SF_ERROR_SINGULAR, NULL);
        return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }
    return std::exp(loggamma(z));
}

XSF_HOST_DEVICE inline std::complex<float> gamma(std::complex<float> z) {
    return static_cast<std::complex<float>>(gamma(static_cast<std::complex<double>>(z)));
}

template <typename T>
T gamma_ratio(T a, T b) {
    return std::tgamma(a) / std::tgamma(b);
}

} // namespace xsf
