#pragma once

#include "cephes/gamma.h"
#include "specfun/specfun.h"

namespace special {

template <typename T>
SPECFUN_HOST_DEVICE T gamma(T x) {
    return cephes::Gamma(x);
}

template <typename T>
std::complex<T> gamma(std::complex<T> z) {
    return specfun::cgama(z, 1);
}

template <typename T>
SPECFUN_HOST_DEVICE T gammaln(T x) {
    return cephes::lgam(x);
}

template <typename T>
std::complex<T> gammaln(std::complex<T> z) {
    return specfun::cgama(z, 0);
}

} // namespace special
