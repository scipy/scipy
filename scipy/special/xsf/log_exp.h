#pragma once

#include <cmath>

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

// Compute pow(1+x,y) accurately for real x and y.
template <typename T>
T pow1p_impl(T x, T y) {
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
    if (x < -1) {
        if (std::fmod(y, T(1)) != T(0)) { // y is not an integer
            // TODO: signal invalid operation exception
            return std::numeric_limits<T>::quiet_NaN();
        }
        // TODO: maybe use (1+x)^y == [(1+x)^2]^(y/2)
        return std::pow(1+x, y);
    }

    // Now x, y are finite and not equal to 0 or +/-1, and x > -1.
    // To compute (1+x)^y, write (1+x) == (s+t) where s is equal to (1+x)
    // rounded toward 1, and t is the (exact) rounding error.
    T s, t;
    if (x < 0) {
        s = T(1) + x;
        t = x - (s - T(1));
        if (t > 0) {
            s = std::nextafter(s, T(1));
            t = x - (s - T(1));
        }
    } else if (x < 1) {
        s = T(1) + x;
        t = x - (s - T(1));
        if (t < 0) {
            s = std::nextafter(s, T(0));
            t = x - (s - T(1));
        }
    } else {
        s = x + T(1);
        t = T(1) - (s - x);
        if (t < 0) {
            s = std::nextafter(s, T(0));
            t = T(1) - (s - x);
        }
    }

    // Because x > -1 and s is rounded toward 1, s is guaranteed to be > 0.
    // Then (1+x)^y == (s+t)^y == (s^y)*((1+u)^y), where u := t / s.
    // It can be shown that either both terms <= 1 or both >= 1, so
    // if the first term over/underflows, then the result over/underflows.
    T u = t / s;
    T term1 = std::pow(s, y);
    if (term1 == T(0) || std::isinf(term1)) {
        return term1;
    }

    // (1+u)^y == exp(y*log(1+u)).  Since u is close to machine epsilon,
    // log(1+u) ~= u.  Let y*u == z+w, where z is the rounded result and
    // w is the rounding error.  This improves accuracy when y is large.
    // Then exp(y*u) == exp(z)*exp(w).
    T z = y * u;
    T w = std::fma(y, u, -z);
    T term2 = std::exp(z) * std::exp(w);
    return term1 * term2;
}

inline double pow1p(double x, double y) {
    return pow1p_impl(x, y);
}

inline float pow1p(float x, float y) {
    return pow1p_impl(static_cast<double>(x), static_cast<double>(y));
}

} // namespace xsf
