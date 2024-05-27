// Functions related to distance metrics.

#pragma once

#include "config.h"

namespace special {

SPECFUN_HOST_DEVICE inline double js_div(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) {
        return std::numeric_limits<double>::quiet_NaN();
    } else if (!(a >= 0 && a < INFINITY && b >= 0 && b < INFINITY)) {
        return INFINITY;
    } else if (a == 0 || b == 0) {
        return (a + b) * (0.5 * std::log(2.0)) + 0.0; // avoid -0.0
    } else {
        const double c = (a + b == INFINITY) ? 0.5*a+0.5*b : 0.5*(a+b);
        const double t = (a + b == INFINITY) ? 0.5*(a-b)/c : (a-b)/(a+b);
        if (std::abs(t) <= 0.5) {
            return c * (t * std::atanh(t) + 0.5 * std::log1p(-t * t));
        } else {
            return 0.5 * (a * std::log(a / c) + b * std::log(b / c));
        }
    }
}

SPECFUN_HOST_DEVICE inline float js_div(float a, float b) {
    return js_div(static_cast<double>(a), static_cast<double>(b));
}

} // namespace special
