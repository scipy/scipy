#pragma once

#include "specfun.h"

namespace special {

template <typename T>
T itstruve0(T x) {
    if (x < 0) {
        x = -x;
    }

    T out = specfun::itsh0(x);
    SPECFUN_CONVINF("itstruve0", out);
    return out;
}

template <typename T>
T it2struve0(T x) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }

    T out = specfun::itth0(x);
    SPECFUN_CONVINF("it2struve0", out);
    if (flag) {
        out = M_PI - out;
    }
    return out;
}

template <typename T>
T itmodstruve0(T x) {
    if (x < 0) {
        x = -x;
    }

    T out = specfun::itsl0(x);
    SPECFUN_CONVINF("itmodstruve0", out);
    return out;
}

} // namespace special
