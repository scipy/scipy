/* Translated into C++ by SciPy developers in 2024. */
#pragma once

#include "../config.h"
#include "../error.h"

#include "ndtri.h"

namespace special {
namespace cephes {

    /*
     * Inverse of the error function.
     *
     * Computes the inverse of the error function on the restricted domain
     * -1 < y < 1. This restriction ensures the existence of a unique result
     * such that erf(erfinv(y)) = y.
     */
    SPECFUN_HOST_DEVICE inline double erfinv(double y) {
        constexpr double domain_lb = -1;
        constexpr double domain_ub = 1;

        constexpr double thresh = 1e-7;

        /*
         * For small arguments, use the Taylor expansion
         * erf(y) = 2/\sqrt{\pi} (y - y^3 / 3 + O(y^5)),    y\to 0
         * where we only retain the linear term.
         * Otherwise, y + 1 loses precision for |y| << 1.
         */
        if ((-thresh < y) && (y < thresh)) {
            return y / M_2_SQRTPI;
        }
        if ((domain_lb < y) && (y < domain_ub)) {
            return ndtri(0.5 * (y + 1)) * M_SQRT1_2;
        } else if (y == domain_lb) {
            return -std::numeric_limits<double>::infinity();
        } else if (y == domain_ub) {
            return std::numeric_limits<double>::infinity();
        } else if (std::isnan(y)) {
            set_error("erfinv", SF_ERROR_DOMAIN, NULL);
            return y;
        } else {
            set_error("erfinv", SF_ERROR_DOMAIN, NULL);
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    /*
     * Inverse of the complementary error function.
     *
     * Computes the inverse of the complimentary error function on the restricted
     * domain 0 < y < 2. This restriction ensures the existence of a unique result
     * such that erfc(erfcinv(y)) = y.
     */
    SPECFUN_HOST_DEVICE inline double erfcinv(double y) {
        constexpr double domain_lb = 0;
        constexpr double domain_ub = 2;

        if ((domain_lb < y) && (y < domain_ub)) {
            return -ndtri(0.5 * y) * M_SQRT1_2;
        } else if (y == domain_lb) {
            return std::numeric_limits<double>::infinity();
        } else if (y == domain_ub) {
            return -std::numeric_limits<double>::infinity();
        } else if (std::isnan(y)) {
            set_error("erfcinv", SF_ERROR_DOMAIN, NULL);
            return y;
        } else {
            set_error("erfcinv", SF_ERROR_DOMAIN, NULL);
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

} // namespace cephes
} // namespace special
