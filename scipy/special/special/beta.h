#pragma once

#include "cephes/beta.h"

namespace special {

template <typename T>
SPECFUN_HOST_DEVICE T beta(T a, T b) {
    return cephes::beta(a, b);
}

template <>
SPECFUN_HOST_DEVICE inline float beta(float af, float bf) {
    double a = af;
    double b = bf;

    return beta(a, b);
}

template <typename T>
SPECFUN_HOST_DEVICE T betaln(T a, T b) {
    return cephes::lbeta(a, b);
}

template <>
SPECFUN_HOST_DEVICE inline float betaln(float af, float bf) {
    double a = af;
    double b = bf;

    return betaln(a, b);
}

} // namespace special