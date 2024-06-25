#pragma once

#include "cephes/zeta.h"

namespace special {

template <typename T>
SPECFUN_HOST_DEVICE T zeta(T x, T q) {
    return cephes::zeta(x, q);
}

template <>
SPECFUN_HOST_DEVICE inline float zeta(float xf, float qf) {
    double x = xf;
    double q = qf;

    return zeta(x, q);
}

} // namespace special