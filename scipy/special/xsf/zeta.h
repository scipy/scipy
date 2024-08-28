#pragma once

#include "cephes/zeta.h"
#include "cephes/zetac.h"

namespace xsf {

template <typename T>
XSF_HOST_DEVICE T zeta(T x, T q) {
    return cephes::zeta(x, q);
}

template <>
XSF_HOST_DEVICE inline float zeta(float x, float q) {
    return zeta(static_cast<double>(x), static_cast<double>(q));
}

template <typename T>
XSF_HOST_DEVICE T zetac(T x) {
    return cephes::zetac(x);
}

template <>
XSF_HOST_DEVICE inline float zetac(float x) {
    return zetac(static_cast<double>(x));
}

} // namespace xsf
