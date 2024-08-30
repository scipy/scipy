
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
    if (!(a > 0) || std::isinf(a)) {
        set_error("gdtrib", SF_ERROR_DOMAIN, "Input parameter a is out of range");
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (!(x >= 0) || std::isinf(x)) {
        set_error("gdtrib", SF_ERROR_DOMAIN, "Input parameter x is out of range");
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (x == 0.0) {
        if (p == 0.0) {
            set_error("gdtrib", SF_ERROR_DOMAIN, "Indeterminate result for (x, p) == (0, 0).");
            return std::numeric_limits<double>::quiet_NaN();
        }
        /* gdtrib(a, p, x) tends to 0 as x -> 0 when p > 0 */
        return 0.0;
    }
    if (p == 0.0) {
        /* gdtrib(a, p, x) tends to infinity as p -> 0 from the right when x > 0. */
        set_error("gdtrib", SF_ERROR_SINGULAR, NULL);
        return std::numeric_limits<double>::infinity();
    }
    if (p == 1.0) {
        /* gdtrib(a, p, x) tends to 0 as p -> 1.0 from the left when x > 0. */
        return 0.0;
    }

    double q = 1.0 - p;
    double (*func)(double, std::tuple<double, double, double>);
    std::tuple<double, double, double> args;
    if (p <= q) {
        args = {a, p, x};
        func = +[](double b, std::tuple<double, double, double> args) -> double {
            auto [a, p, x] = args;
            return cephes::igam(b, a * x) - p;
        };
    } else {
        args = {a, q, x};
        func = +[](double b, std::tuple<double, double, double> args) -> double {
            auto [a, q, x] = args;
            return q - cephes::igamc(b, a * x);
        };
    }
    double lower_bound = std::numeric_limits<double>::min();
    double upper_bound = std::numeric_limits<double>::max();
    auto [xl, xr, f_xl, f_xr, bracket_status] = detail::bracket_root_for_cdf_inversion(
        func, 1.0, lower_bound, upper_bound, -0.875, 7.0, 0.125, 8, false, 342, args
    );
    if (bracket_status == 1) {
        set_error("gdtrib", SF_ERROR_UNDERFLOW, NULL);
        return 0.0;
    }
    if (bracket_status == 2) {
        set_error("gdtrib", SF_ERROR_OVERFLOW, NULL);
        return std::numeric_limits<double>::infinity();
    }
    if (bracket_status >= 3) {
        set_error("gdtrib", SF_ERROR_OTHER, "Computational Error");
        ;
        return std::numeric_limits<double>::quiet_NaN();
    }
    auto [result, root_status] = detail::find_root_chandrupatla(
        func, xl, xr, f_xl, f_xr, std::numeric_limits<double>::epsilon(), 1e-100, 100, args
    );
    if (root_status) {
        /* The root finding return should only fail if there's a bug in our code. */
        set_error("gdtrib", SF_ERROR_OTHER, "Computational Error, (%.17g, %.17g, %.17g)", a, p, x);
        ;
        return std::numeric_limits<double>::quiet_NaN();
    }
    return result;
}

} // namespace xsf
