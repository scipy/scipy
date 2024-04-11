// Numerically stable computation of iv(v+1, x) / iv(v, x)

#pragma once

#include "special/tools.h"
#include <algorithm>
#include <iterator>
#include <cstddef>
#include <cmath>
#include <utility>

class IvRatioCFIterator {
  /* Generates the Perron continued fraction for iv(v, x) / iv(v-1, x), for
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

public:
    using iterator_category = std::input_iterator_tag;
    using value_type = std::pair<double, double>;
    using difference_type = std::ptrdiff_t;
    using pointer = const value_type *;
    using reference = const value_type &;

    // It is assumed that v >= 1, x >= 0, and both are finite.
    IvRatioCFIterator(double v, double x) noexcept {
        int e;
        std::frexp(std::max(v, x), &e);
        double c = std::ldexp(1, 1-e); // rescaling multiplier

        double vc = v * c;
        double xc = x * c;
        _frac = {xc, 2*vc+xc};
        _a0 = -(2*vc-c)*xc;
        _as = -2*c*xc;
        _b0 = 2*(vc+xc);
        _bs = c;
        _k = 0;
    }

    reference operator*() const { return _frac; }

    IvRatioCFIterator & operator++() /* prefix ++ */ {
        ++_k;
        _frac = {std::fma(_k, _as, _a0), std::fma(_k, _bs, _b0)};
        return *this;
    }

    IvRatioCFIterator operator++(int) /* postfix ++ */ {
        IvRatioCFIterator self = *this;
        ++_k;
        _frac = {std::fma(_k, _as, _a0), std::fma(_k, _bs, _b0)};
        return self;
    }

private:
    std::pair<double, double> _frac;  // current fraction
    double _a0, _as;  // a[k] == _a0 + _as*k, k >= 1
    double _b0, _bs;  // b[k] == _b0 + _bs*k, k >= 1
    double _k;        // current index (0-based)
};

inline double iv_ratio(double v, double x) {
//
//    std::pair<double, double> cf[] = {
//        {4., 1.}, {1., 3.}, {4., 5.}, {9., 7.}, {16., 9.},
//        {25., 11.}, {36, 13}, {49, 15}, {64, 17}, {81, 19},
//        {100, 21}, {0, 23},
//    };
//    return special::detail::continued_fraction_eval_series(
////    return special::detail::continued_fraction_eval_lentz(
//        0., &cf[0], 0., 5000).first;

    if (std::isnan(v) || std::isnan(x)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (v < 0 || x < 0) {
        // TODO: raise domain error
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (std::isinf(v) && std::isinf(x)) {
        // TBC
        // There is not a unique limit as both v and x tends to infinity.
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (std::isinf(v)) {
        return 0.0;
    }
    if (std::isinf(x)) {
        return 1.0;
    }

    // Now v >= 0 and x >= 0 and both are finite.  Evaluate the c.f. using
    // series method.  The termination condition is
    IvRatioCFIterator cf(v+1, x);
    double tol = std::numeric_limits<double>::epsilon() * (v+0.75);
    tol = tol / std::max(tol, x); // prevent overflow
    auto result = special::detail::continued_fraction_eval_series(
        special::detail::KahanSummer<double>(), cf, tol, 5000);
//        0.0, cf, tol*0, 5000);
    if (result.second == 0) { // failed to converge; not expected
        return std::numeric_limits<double>::quiet_NaN();
    }
    return static_cast<double>(result.first);
}
