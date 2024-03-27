#pragma once

#include "cephes/gamma.h"

namespace special {

template <typename T>
SPECFUN_HOST_DEVICE T gamma(T x) {
    return cephes::Gamma(x);
}

template <>
SPECFUN_HOST_DEVICE inline float gamma(float xf) {
    double x = xf;

    return gamma(x);
}

template <typename T>
SPECFUN_HOST_DEVICE T gammaln(T x);

template <>
SPECFUN_HOST_DEVICE inline double gammaln(double x) {
    return cephes::lgam(x);
}

template <>
SPECFUN_HOST_DEVICE inline float gammaln(float xf) {
    double x = xf;

    return gammaln(x);
}

} // namespace special
