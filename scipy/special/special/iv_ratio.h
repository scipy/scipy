// Numerically stable computation of iv(v+1, x) / iv(v, x)

#pragma once

#include "config.h"
#include "tools.h"
#include "error.h"

namespace special {

/* Generates the "tail" of Perron's continued fraction for `iv(v,x)/iv(v-1,x)`
 * for v >= 1 and x >= 0.
 *
 * The Perron continued fraction is studied in [1].  It is given by
 *
 *         iv(v, x)      x    -(2v+1)x   -(2v+3)x   -(2v+5)x
 *   R := --------- = ------ ---------- ---------- ---------- ...
 *        iv(v-1,x)   2v+x + 2(v+x)+1 + 2(v+x)+2 + 2(v+x)+3 +
 *
 * Rearrange the expression by making an equivalent transform to prevent
 * floating point overflow and extracting the first fraction to simplify
 * the recurrence relation.  This leads to
 *
 *        xc                -(2vc+c)(xc) -(2vc+3c)(xc) -(2vc+5c)(xc)
 *   R = -----,  fc = 2vc + ------------ ------------- ------------- ...
 *       xc+fc              2(vc+xc)+c + 2(vc+xc)+2c + 2(vc+xc)+3c +
 *
 * This class generates the fractions of fc after 2vc.
 *
 * [1] Gautschi, W. and Slavik, J. (1978). "On the computation of modified
 *     Bessel function ratios." Mathematics of Computation, 32(143):865-875.
 */
struct IvRatioCFTailGenerator {

    // It is assumed that v >= 1, x >= 0, c > 0, and all are finite.
    IvRatioCFTailGenerator(double vc, double xc, double c) noexcept {
        a0_ = -(2*vc-c)*xc;
        as_ = -2*c*xc;
        b0_ = 2*(vc+xc);
        bs_ = c;
        k_ = 0;
    }

    std::pair<double, double> operator()() {
        ++k_;
        return {std::fma(k_, as_, a0_), std::fma(k_, bs_, b0_)};
    }

private:
    double a0_, as_;  // a[k] == a0 + as*k, k >= 1
    double b0_, bs_;  // b[k] == b0 + bs*k, k >= 1
    std::uint64_t k_; // current index
};

SPECFUN_HOST_DEVICE inline double iv_ratio(double v, double x) {

    if (std::isnan(v) || std::isnan(x)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (v < 1 || x < 0) {
        set_error("iv_ratio", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (std::isinf(v) && std::isinf(x)) {
        // There is not a unique limit as both v and x tends to infinity.
        set_error("iv_ratio", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (x == 0.0) {
        // If x is +/-0.0, return +/-0.0 to agree with the limiting behavior.
        return x;
    }
    if (std::isinf(v)) {
        return 0.0;
    }
    if (std::isinf(x)) {
        return 1.0;
    }

    // Now v >= 1 and x >= 0 and both are finite.
    int e;
    std::frexp(std::fmax(v, x), &e);
    double c = std::ldexp(1, 2-e); // rescaling multiplier
    double vc = v * c;
    double xc = x * c;

    IvRatioCFTailGenerator cf(vc, xc, c);
    auto [fc, terms] = detail::series_eval_kahan(
        detail::continued_fraction_series(cf),
        std::numeric_limits<double>::epsilon() * 0.5,
        1000,
        2*vc);

    if (terms == 0) { // failed to converge; should not happen
        set_error("iv_ratio", SF_ERROR_NO_RESULT, NULL);
        return std::numeric_limits<double>::quiet_NaN();
    }

    return xc / (xc + fc);
}

SPECFUN_HOST_DEVICE inline float iv_ratio(float v, float x) {
    return iv_ratio(static_cast<double>(v), static_cast<double>(x));
}

} // namespace special
