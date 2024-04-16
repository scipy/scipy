// Numerically stable computation of iv(v+1, x) / iv(v, x)

#pragma once

#include "special/tools.h"
#include "special/error.h"
#include <cstdint>    // for std::uint64_t
#include <cmath>      // for std::frexp, std::ldexp, std::fmax
#include <utility>    // for std::pair

/* Generates the Perron continued fraction for `iv(v, x) / iv(v-1, x)` for
 * v >= 1 and x >= 0.
 *
 * The formulas are given in [1].  We additionally perform an equivalent
 * transform of the c.f. to avoid overflow.
 *
 *    iv(v, x)      xc    -(2vc+c)(xc) -(2vc+3c)(xc) -(2vc+5c)(xc)
 *   --------- = -------- ------------ ------------- ------------- ...
 *   iv(v-1,x)   2vc+xc + 2(vc+xc)+c + 2(vc+xc)+2c + 2(vc+xc)+3c +
 *
 * [1] Gautschi, W. and Slavik, J. (1978). "On the computation of modified
 * Bessel function ratios." Mathematics of Computation, 32(143):865-875.
 */
class IvRatioCFGenerator {
    using frac_type = std::pair<double, double>;

public:

    // It is assumed that v >= 1, x >= 0, and both are finite.
    IvRatioCFGenerator(double v, double x) noexcept {
        int e;
        std::frexp(std::fmax(v, x), &e);
        double c = std::ldexp(1, 1-e); // rescaling multiplier

        double vc = v * c;
        double xc = x * c;
        frac_ = {xc, 2*vc+xc};
        a0_ = -(2*vc-c)*xc;
        as_ = -2*c*xc;
        b0_ = 2*(vc+xc);
        bs_ = c;
        k_ = 0;
    }

    frac_type operator()() {
        frac_type frac = frac_;
        double k = ++k_;
        frac_ = {std::fma(k, as_, a0_), std::fma(k, bs_, b0_)};
        return frac;
    }

private:
    frac_type frac_;  // current fraction
    double a0_, as_;  // a[k] == a0 + as*k, k >= 1
    double b0_, bs_;  // b[k] == b0 + bs*k, k >= 1
    std::uint64_t k_; // current index (0-based)
};

inline double iv_ratio(double v, double x) {

    if (std::isnan(v) || std::isnan(x)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (v < 0 || x < 0) {
        special::set_error("iv_ratio", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<double>::signaling_NaN();
    }
    if (std::isinf(v) && std::isinf(x)) {
        // There is not a unique limit as both v and x tends to infinity.
        special::set_error("iv_ratio", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<double>::signaling_NaN();
    }
    if (std::isinf(v)) {
        return 0.0;
    }
    if (std::isinf(x)) {
        return 1.0;
    }

    // Now v >= 0 and x >= 0 and both are finite.
    IvRatioCFGenerator cf(v+1, x);

    auto [result, terms] = special::detail::series_eval_kahan(
        special::detail::continued_fraction_series(cf),
        std::numeric_limits<double>::epsilon() * 0.5,
        1000);

    if (terms == 0) { // failed to converge; should not happen
        special::set_error("iv_ratio", SF_ERROR_NO_RESULT, NULL);
        return std::numeric_limits<double>::signaling_NaN();
    }

    return result;
}
