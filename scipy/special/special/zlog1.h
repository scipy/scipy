#pragma once

/* Author: Josh Wilson 2016
 * Translated from Cython to C++ by SciPy developers 2023.
 */

#include "config.h"

namespace special {
namespace detail {

    SPECFUN_HOST_DEVICE inline std::complex<double> zlog1(std::complex<double> z) {
        /* Compute log, paying special attention to accuracy around 1. We
         * implement this ourselves because some systems (most notably the
         * Travis CI machines) are weak in this regime. */
        const double tol = 2.220446092504131e-16; // machine epsilon
        std::complex<double> coeff = -1.0;
        std::complex<double> res = 0.0;

        if (std::abs(z - 1.0) > 0.1) {
            return std::log(z);
        }

        z -= 1.0;
        for (int n = 1; n < 17; n++) {
            coeff *= -z;
            res += coeff / static_cast<double>(n);
            if (std::abs(res / coeff) < tol) {
                break;
            }
        }
        return res;
    }
} // namespace detail
} // namespace special
