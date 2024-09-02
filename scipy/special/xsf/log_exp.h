#pragma once

#include <cmath>
#include <cstdio>
#include <cfloat>

//#ifndef FLT_EVAL_METHOD
//#error FLT_EVAL_METHOD not defined; too old compiler?
//#endif
//
//#if FLT_EVAL_METHOD != 0
//#error The pow1p function requires FLT_EVAL_METHOD to be zero
//#endif

#include "config.h"

namespace xsf {

template <typename T>
T expit(T x) {
    return 1 / (1 + std::exp(-x));
};

inline double exprel(double x) {
    if (std::abs(x) < std::numeric_limits<double>::epsilon()) {
        return 1;
    }

    if (x > 717) { // near log(DBL_MAX)
        return std::numeric_limits<double>::infinity();
    }

    return std::expm1(x) / x;
}

inline float exprel(float x) { return exprel(static_cast<double>(x)); }

template <typename T>
T logit(T x) {
    // The standard formula is log(x/(1 - x)), but this expression
    // loses precision near x=0.5, as does log(x) - log1p(-x).
    // We use the standard formula away from p=0.5, and use
    // log1p(2*(x - 0.5)) - log1p(-2*(x - 0.5)) around p=0.5, which
    // provides very good precision in this interval.
    if (x < 0.3 || x > 0.65) {
        return std::log(x/(1 - x));
    }
    else {
        T s = 2*(x - 0.5);
        return std::log1p(s) - std::log1p(-s);
    }
};

//
// The logistic sigmoid function 'expit' is
//
//     S(x) = 1/(1 + exp(-x))     = exp(x)/(exp(x) + 1)
//
// so
//
// log S(x) = -log(1 + exp(-x))   = x - log(exp(x) + 1)
//          = -log1p(exp(-x))     = x - log1p(exp(x))
//
// By using -log1p(exp(-x)) for x >= 0 and x - log1p(exp(x))
// for x < 0, we extend the range of x values for which we
// obtain accurate results (compared to the naive implementation
// log(expit(x))).
//
template <typename T>
T log_expit(T x) {
    if (x < 0) {
        return x - std::log1p(std::exp(x));
    }

    return -std::log1p(std::exp(-x));
};

#ifndef FLT_EVAL_METHOD
#error Missing definition of FLT_EVAL_METHOD
#endif

//#if FLT_EVAL_METHOD == 0
// Floating point arithmetic in target precision
template <typename T>
void fast_two_sum(const T &a, const T &b, T &s, T &t) {
    s = a + b;
    t = s - a;
    t = b - t;
}
//#else
//template <typename T>
//void fast_two_sum(const T &a, const T &b, T &s, T &t) {
//    using std::fma;
//    s = fma(T(1), a, b);
//    t = fma(T(-1), a, s);
//    t = fma(T(-1), t, b);
//}
//#endif

template <typename T>
void two_sum(const T &a, const T &b, T &s, T &t) {
    using std::abs;
    if (abs(a) >= abs(b)) {
        fast_two_sum(a, b, s, t);
    } else {
        fast_two_sum(b, a, s, t);
    }
}

// Compute pow(1+x,y) accurately for real x and y.
template <typename T>
T pow1p_impl(T x, T y) {
    auto out = stderr;
    fprintf(out, "[pow1p] ========\n");
    fprintf(out, "[pow1p] x=%.16e y=%.16e\n", x, y);

    // The special values follow the spec of the pow() function defined in
    // IEEE-754, Section 9.2.1.  The order of the `if` statements matters.
    if (y == T(0)) { // pow(x, +/-0)
        return T(1);
    }
    if (x == T(-1)) { // pow(+/-0, y)
        if (std::isfinite(y) && y < 0) {
            // TODO: signal divideByZero exception
        }
        return std::pow(T(0), y);
    }
    if (x == T(-2) && std::isinf(y)) { // pow(-1, +/-inf)
        return T(1);
    }
    if (x == T(0)) { // pow(+1, y)
        return T(1);
    }
    if (y == std::numeric_limits<T>::infinity()) { // pow(x, +inf)
        return (x < 0 && x > -2) ? T(0) : std::numeric_limits<T>::infinity();
    }
    if (y == -std::numeric_limits<T>::infinity()) { // pow(x, -inf)
        return (x < 0 && x > -2) ? std::numeric_limits<T>::infinity() : T(0);
    }
    if (std::isinf(x)) { // pow(+/-inf, y)
        return std::pow(x, y);
    }

    // Up to this point, (1+x) = {+/-0, +/-1, +/-inf} have been handled, and
    // and y = {+/-0, +/-inf} have been handled.  Next, we handle `nan` and
    // y = {+/-1}.
    if (std::isnan(x) || std::isnan(y)) {
        return std::numeric_limits<T>::quiet_NaN();
    }
    if (y == T(1)) {
        return T(1) + x;
    }
    if (y == T(-1)) {
        return T(1) / (T(1) + x); // guaranteed no overflow
    }

    // Handle (1+x) < 0
    T sign = T(1);
    if (x < -1) {
        if (std::fmod(y, T(1)) != T(0)) { // y is not an integer
            // TODO: signal invalid operation exception
            return std::numeric_limits<T>::quiet_NaN();
        }
        constexpr T two_over_eps = T(2) / std::numeric_limits<T>::epsilon();
        // For 1 < (-x) <= 2^p = 2/eps, (1+x) is exact;
        // For (-x) > 2^(2p) = (2/eps)^2, (1+x)^y ~= x^y.
        // In both cases, we may delegate to pow() without loss of accuracy.
        if (x >= -two_over_eps || x < -two_over_eps*two_over_eps) {
            return std::pow(1+x, y);
        }
        // Otherwise, we take the "slow path" as in the case (1+x) > 0.
        sign = (std::fmod(y, T(2)) == T(0)) ? T(1) : T(-1);
    }

    // Now x, y are finite and not equal to 0 or +/-1.
    // To compute (1+x)^y, write |1+x| == (s+t) where s is equal to |1+x|
    // rounded toward 1, and t is the (exact) rounding error.
    T s, t;
    two_sum(T(1), x, s, t);
    if (t != T(0) && std::signbit(t) != std::signbit(x)) {
        T delta = s - std::nextafter(s, T(1)); // exact
        s -= delta; // exact
        t += delta; // exact
    }
    if (x < -1) {
        s = -s;
        t = -t;
    }
    fprintf(out, "[pow1p] s=%.16e t=%.16e\n", s, t);

    // Because x > -1 and s is rounded toward 1, s is guaranteed to be > 0.

    // Write (1+x)^y == (s+t)^y == (s^y)*((1+t/s)^y) == term1*term2.
    // It can be shown that either both terms <= 1 or both >= 1; so
    // if the first term over/underflows, then the result over/underflows.
    // And of course, term2 == 1 if t == 0.
    T term1 = sign * std::pow(s, y);
    fprintf(out, "[pow1p] term1=%.16e\n", term1);

    if (t == T(0) || term1 == T(0) || std::isinf(term1)) {
        return term1;
    }

    // (1+t/s)^y == exp(y*log(1+t/s)).  The relative error of the result is
    // equal to the absolute error of the exponent (plus the relative error
    // of 'exp').  Therefore, when the exponent is small, it is accurately
    // evaluated to machine epsilon using T arithmetic.  In addition, since
    // t/s <= epsilon, log(1+t/s) is well approximated by t/s to first order.
    T u = t / s;
    T w = y * u;
    T term2 = std::exp(w);
    fprintf(out, "[pow1p] term2=%.16e\n", term2);
    if (std::abs(w) <= 0.5) {
        return term1 * term2;
    }

    // Now y*log(1+t/s) is large, and its relative error is "magnified" by
    // the exponent.  To reduce the error, we use double-T arithmetic, and
    // expand log(1+t/s) to second order.

    // (u + uu) ~= t/s.
    T r1 = std::fma(-u, s, t);
    T uu = r1 / s;

    // (u + vv) ~= log(1+(u+uu)) ~= log(1+t/s).
    T vv = std::fma(-0.5*u, u, uu);

    // (w + ww) ~= y*(u+vv) ~= y*log(1+t/s).
    T r2 = std::fma(y, u, -w);
    T ww = std::fma(y, vv, r2);

    // TODO: maybe ww is small enough such that exp(ww) ~= 1+ww.
    T term3 = std::exp(ww);
    fprintf(out, "[pow1p] term3=%.16e\n", term3);
    return term1 * term2 * term3;
}

inline double pow1p(double x, double y) {
    return pow1p_impl(x, y);
}

inline float pow1p(float x, float y) {
    return pow1p_impl(static_cast<double>(x), static_cast<double>(y));
}

} // namespace xsf
