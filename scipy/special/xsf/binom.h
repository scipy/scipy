/* Translated from Cython into C++ by SciPy developers in 2024.
 *
 * Original authors: Pauli Virtanen, Eric Moore
 */

// Binomial coefficient

#pragma once

#include "config.h"
#include "digamma.h"

#include "cephes/beta.h"
#include "cephes/gamma.h"
#include "cephes/unity.h"

namespace xsf {

XSF_HOST_DEVICE inline double binom(double n, double k) {
    double kx, nx, num, den, dk, sgn;

    if (n < 0) {
        nx = std::floor(n);
        if (n == nx) {
            // Undefined
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    kx = std::floor(k);
    if (k == kx && (std::abs(n) > 1E-8 || n == 0)) {
        /* Integer case: use multiplication formula for less rounding
         * error for cases where the result is an integer.
         *
         * This cannot be used for small nonzero n due to loss of
         * precision. */
        nx = std::floor(n);
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
                if (std::abs(num) > 1E50) {
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
        return std::exp(-cephes::lbeta(1 + n - k, 1 + k) - std::log(n + 1));
    }
    if (k > 1E8 * std::abs(n)) {
        // avoid loss of precision
        num = cephes::Gamma(1 + n) / std::abs(k) + cephes::Gamma(1 + n) * n / (2 * k * k); // + ...
        num /= M_PI * std::pow(std::abs(k), n);
        if (k > 0) {
            kx = std::floor(k);
            if (static_cast<int>(kx) == kx) {
                dk = k - kx;
                sgn = (static_cast<int>(kx) % 2 == 0) ? 1 : -1;
            } else {
                dk = k;
                sgn = 1;
            }
            return num * std::sin((dk - n) * M_PI) * sgn;
        }
        kx = std::floor(k);
        if (static_cast<int>(kx) == kx) {
            return 0;
        }
        return num * std::sin(k * M_PI);
    }
    return 1 / (n + 1) / cephes::beta(1 + n - k, 1 + k);
}

XSF_HOST_DEVICE inline float binom(float n, float k) {
    return binom(static_cast<double>(n), static_cast<double>(k));
}

XSF_HOST_DEVICE inline double binomln(double n, double k) {
    /* Natural logarithm of absolute value of binomial coefficient */
    double a = n - k + 1;
    bool denominator_pole = (k <= - 1 && k == std::floor(k)) || (a <= 0 && a == std::floor(a));
    if (n <= -1 && n == std::floor(n)) {
	return denominator_pole ? std::numeric_limits<double>::quiet_NaN() : std::numeric_limits<double>::infinity();
    }
    if (denominator_pole) {
	return -std::numeric_limits<double>::infinity();
    }
    if (n == k) {
	return 0;
    }
    if (std::abs(k) > 1e-7) {
	double result = n >= -1 ? -std::log1p(n) : -std::log(-n - 1);
        result -= cephes::lbeta(1 + n - k, 1 + k);
	return result;
    }
    if (std::abs(n) > 1e-7) {
	/* binomln(n, k) = gammaln(n + 1) - gammaln(n + 1 - k) - gammaln(1 + k)
	 * gammaln(1 + k) can be computed accurately for small k using lgam1p.
	 * We can then use the Taylor series expansion
	 * gammaln(n + 1 - k) = gammaln(n + 1) - k*digamma(n + 1) +
	 * k^2 * polygamma(3, n + 1) / 2 + ... + (-k)^j * polygamma(j, n + 1) / j! + ...
	 */
	return -cephes::lgam1p(k) + k*digamma(n + 1);
    }
    /* n and k are both small. */
    return cephes::lgam1p(n) - cephes::lgam1p(k) - cephes::lgam1p(n - k);
}

} // namespace xsf
