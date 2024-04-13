#pragma once

#include <cmath>

#include "config.h"

namespace special {

template <typename T>
T expit(T x) {
    return 1 / (1 + std::exp(-x));
};

template <typename T>
T exprel(T x) {
    if (std::abs(x) < std::numeric_limits<T>::epsilon()) {
        return 1;
    }

    if (x >= M_LN10 * std::numeric_limits<T>::max_exponent10) {
        return std::numeric_limits<T>::infinity(); // std::log(10) * std::numeric_limits<T>::max_exponent10 is the
                                                   // largest x such that std::exp(x) is representable
    }

    return std::expm1(x) / x;
}

template <typename T>
T logit(T x) {
    return std::log(x / (1 - x));
};

//
// The logistic sigmoid function 'expit' is
//
//     S(x) = 1/(1 + exp(-x))     = exp(x)/(exp(x) + 1)
//
// so
//
// log S(x) = -log(1 + exp(-x))   = x - log(exp(x) + 1)
//          = -log1p(exp(-x))     = x - log1p(exp(x))
//
// By using -log1p(exp(-x)) for x >= 0 and x - log1p(exp(x))
// for x < 0, we extend the range of x values for which we
// obtain accurate results (compared to the naive implementation
// log(expit(x))).
//
template <typename T>
T log_expit(T x) {
    if (x < 0) {
        return x - std::log1p(std::exp(x));
    }

    return -std::log1p(std::exp(-x));
};

} // namespace special
