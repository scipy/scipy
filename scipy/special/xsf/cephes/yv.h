/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*
 * Cephes Math Library Release 2.8: June, 2000
 * Copyright 1984, 1987, 2000 by Stephen L. Moshier
 */
#pragma once

#include "../config.h"
#include "../error.h"

#include "const.h"
#include "jv.h"
#include "yn.h"

namespace xsf {
namespace cephes {

    /*
     * Bessel function of noninteger order
     */
    XSF_HOST_DEVICE inline double yv(double v, double x) {
        double y, t;
        int n;

        n = v;
        if (n == v) {
            y = yn(n, x);
            return (y);
        } else if (v == std::floor(v)) {
            /* Zero in denominator. */
            set_error("yv", SF_ERROR_DOMAIN, NULL);
            return std::numeric_limits<double>::quiet_NaN();
        }

        t = M_PI * v;
        y = (std::cos(t) * jv(v, x) - jv(-v, x)) / std::sin(t);

        if (std::isinf(y)) {
            if (v > 0) {
                set_error("yv", SF_ERROR_OVERFLOW, NULL);
                return -std::numeric_limits<double>::infinity();
            } else if (v < -1e10) {
                /* Whether it's +inf or -inf is numerically ill-defined. */
                set_error("yv", SF_ERROR_DOMAIN, NULL);
                return std::numeric_limits<double>::quiet_NaN();
            }
        }

        return (y);
    }
} // namespace cephes
} // namespace xsf
