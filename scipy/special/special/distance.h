// Functions related to distance metrics.

#pragma once

#include "config.h"

namespace special {

SPECFUN_HOST_DEVICE inline double js_div(double a, double b) {

    if (std::isnan(a)) {
        return a;
    }
    if (std::isnan(b)) {
        return b;
    }
    if (!(a >= 0 && std::isfinite(a) && b >= 0 && std::isfinite(b))) {
        return std::numeric_limits<double>::infinity();
    }

    // Handle the case where at least one input is zero.
    const double HALF_LN2 = 0.34657359027997264; // 0.5*log(2)
    if (a == 0 && b == 0) {
        return 0.0;
    }
    if (a == 0) {
        return b * HALF_LN2;
    }
    if (b == 0) {
        return a * HALF_LN2;
    }

    // Now both inputs are positive.
    const double c = (a+b)/2;
    const double t = (b-a)/(b+a);
    if (std::abs(t) <= 0.5) {
        return c*( t*std::atanh(t) + 0.5*std::log1p(-t*t) );  // fma?
    } else {
        return 0.5*( a*std::log(a/c) + b*std::log(b/c) );
    }
}

SPECFUN_HOST_DEVICE inline float js_div(float a, float b) {
    return js_div(static_cast<double>(a), static_cast<double>(b));
}

} // namespace special
