#pragma once

#include "cephes/igam.h"
#include "config.h"
#include "error.h"
#include "tools.h"

namespace xsf {

XSF_HOST_DEVICE inline double gdtrib(double a, double p, double x) {
    if (std::isnan(p) || std::isnan(a) || std::isnan(x)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (!((0 <= p) && (p <= 1))) {
        set_error("gdtrib", SF_ERROR_DOMAIN, "Input parameter p is out of range");
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (!(a > 0)) {
        set_error("gdtrib", SF_ERROR_DOMAIN, "Input parameter a is out of range");
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (!(x >= 0)) {
        set_error("gdtrib", SF_ERROR_DOMAIN, "Input parameter x is out of range");
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (x == 0.0) {
        if (p == 0.0) {
            /* CDF evaluated at x = 0 will be zero except for degenerate b = 0
             * case. Following our convention, we take the smallest normalized value
             * for b. */
            return std::numeric_limits<double>::min();
        } else if (p == 1.0) {
            /* The CDF evaluated at x = 0 will be one for the degenerate b = 0
             * case, where we get a point distribution at zero. */
            return 0.0;
        }
        // When x = 0, there can be no b for which 0 < p < 1.
        return std::numeric_limits<double>::quiet_NaN();
    }

    double q = 1.0 - p;
    std::function<double(double)> func;
    if (p <= q) {
        func = [p, a, x](double b) { return cephes::igam(b, a * x) - p; };
    } else {
        func = [q, a, x](double b) { return q - cephes::igamc(b, a * x); };
    }
    double lower_bound = std::numeric_limits<double>::min();
    double upper_bound = std::numeric_limits<double>::max();
    auto [x_left, x_right, bracket_status] = detail::bracket_root(func, 1.0, 5.0, lower_bound, upper_bound, 8.0, false);
    if (bracket_status == 1) {
        set_error("gdtrib", SF_ERROR_UNDERFLOW, NULL);
        return 0.0;
    }
    if (bracket_status == 2) {
        set_error("gdtrib", SF_ERROR_OVERFLOW, NULL);
        return std::numeric_limits<double>::infinity();
    }
    if (bracket_status == 3) {
        set_error("gdtrib", SF_ERROR_OTHER, "Computational Error");
        ;
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (func(x_right) == 0.0) {
        return detail::find_closest_root(func, x_left, x_right);
    }
    auto [result, root_status] = detail::find_root_bus_dekker_r(func, x_left, x_right);
    if (root_status) {
        set_error("gdtrib", SF_ERROR_OTHER, "Computational Error");
        ;
        return std::numeric_limits<double>::quiet_NaN();
    }
    return result;
}

} // namespace xsf
