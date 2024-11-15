#pragma once

#include <cmath>

#include "config.h"

namespace xsf {

template <typename T>
T expit(T x) {
    return 1 / (1 + std::exp(-x));
};

inline double exprel(double x) {
    if (std::abs(x) < std::numeric_limits<double>::epsilon()) {
        return 1;
    }

    if (x > 717) { // near log(DBL_MAX)
        return std::numeric_limits<double>::infinity();
    }

    return std::expm1(x) / x;
}

inline float exprel(float x) { return exprel(static_cast<double>(x)); }

template <typename T>
T logit(T x) {
    // The standard formula is log(x/(1 - x)), but this expression
    // loses precision near x=0.5, as does log(x) - log1p(-x).
    // We use the standard formula away from p=0.5, and use
    // log1p(2*(x - 0.5)) - log1p(-2*(x - 0.5)) around p=0.5, which
    // provides very good precision in this interval.
    if (x < 0.3 || x > 0.65) {
        return std::log(x/(1 - x));
    }
    else {
        T s = 2*(x - 0.5);
        return std::log1p(s) - std::log1p(-s);
    }
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

} // namespace xsf
