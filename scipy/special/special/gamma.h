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

template <typename T>
SPECFUN_HOST_DEVICE inline std::complex<T> gamma(std::complex<T> z) {
    return gamma(z);
}

template <typename T>
T gamma_ratio(T a, T b) {
    return std::tgamma(a) / std::tgamma(b);
}

} // namespace special
