// Binomial coefficient
#pragma once

#include "config.h"

#include "cephes/beta.h"
#include "cephes/gamma.h"

namespace special {

SPECFUN_HOST_DEVICE inline double binom(double n, double k) {
    double kx, nx, num, den, dk, sgn;

    if (n < 0) {
        nx = STD::floor(n);
        if (n == nx) {
            // Undefined
            return STD::numeric_limits<double>::quiet_NaN();
        }
    }

    kx = STD::floor(k);
    if (k == kx && (STD::abs(n) > 1E-8 || n == 0)) {
        /* Integer case: use multiplication formula for less rounding
         * error for cases where the result is an integer.
         *
         * This cannot be used for small nonzero n due to loss of
         * precision. */
        nx = STD::floor(n);
        if (nx == n && kx > nx / 2 && nx > 0) {
            // Reduce kx by symmetry
            kx = nx - kx;
        }

        if (kx >= 0 && kx < 20) {
            num = 1.0;
            den = 1.0;
            for (int i = 1; i < 1 + static_cast<int>(kx); i++) {
                num *= i + n - kx;
                den *= i;
                if (STD::abs(num) > 1E50) {
                    num /= den;
                    den = 1.0;
                }
            }
            return num / den;
        }
    }

    // general case
    if (n >= 1E10 * k and k > 0) {
        // avoid under/overflows intermediate results
        return STD::exp(-cephes::lbeta(1 + n - k, 1 + k) - STD::log(n + 1));
    }
    if (k > 1E8 * STD::abs(n)) {
        // avoid loss of precision
        num = cephes::Gamma(1 + n) / STD::abs(k) + cephes::Gamma(1 + n) * n / (2 * k * k); // + ...
        num /= M_PI * STD::pow(STD::abs(k), n);
        if (k > 0) {
            kx = STD::floor(k);
            if (static_cast<int>(kx) == kx) {
                dk = k - kx;
                sgn = (static_cast<int>(kx) % 2 == 0) ? 1 : -1;
            } else {
                dk = k;
                sgn = 1;
            }
            return num * STD::sin((dk - n) * M_PI) * sgn;
        }
        kx = STD::floor(k);
        if (static_cast<int>(kx) == kx) {
            return 0;
        }
        return num * STD::sin(k * M_PI);
    }
    return 1 / (n + 1) / cephes::beta(1 + n - k, 1 + k);
}

} // namespace special
