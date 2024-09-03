#pragma once

#include <cmath>
#include <cstdio>
#include <cfloat>
#include <iostream>

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
#error "Missing definition of FLT_EVAL_METHOD; C++11 or above required"
#endif

// Compute the exact sum of floating point numbers a and b and store the
// result in s and t such that (a + b) == (s + t) exactly and |s| >> |t|.
// Requires |a| >= |b|.
//
// The following implementation assumes "+" and "-" are correctly rounded
// to nearest with ties resolved in any way.
template <typename T>
void fast_two_sum(const T &a, const T &b, T &s, T &t) {
    s = a + b;
    t = s - a;
    t = b - t;
}

// fast_two_sum that handles excess precision by force rounding.
template <typename T>
void fast_two_sum_workaround(const T &a, const T &b, T &s, T &t) {
    using std::fma;
    s = fma(T(1), a, b);
    t = fma(T(-1), a, s);
    t = fma(T(-1), t, b);
}

#if FLT_EVAL_METHOD != 0
inline void fast_two_sum(const float &a, const float &b, float &s, float &t) {
    fast_two_sum_workaround(a, b, s, t);
}
#endif

#if FLT_EVAL_METHOD != 0 && FLT_EVAL_METHOD != 1
inline void fast_two_sum(const double &a, const double &b, double &s, double &t) {
    fast_two_sum_workaround(a, b, s, t);
}
#endif

template <typename T>
void two_sum(const T &a, const T &b, T &s, T &t) {
    using std::abs;
    if (abs(a) >= abs(b)) {
        fast_two_sum(a, b, s, t);
    } else {
        fast_two_sum(b, a, s, t);
    }
}

// Handles special values of pow(1+x,y) according to the spec of the special
// values of pow() defined in IEEE-754, Section 9.2.1, plus a few others.
//
// If the input is handled, stores the (IEEE-754 conformant) return value in
// 'result', and returns 'true'.
//
// If the input is not handled, returns 'false'.  In this case, it is
// guaranteed that (1) x and y are finite, (2) x != 0, x != -1, y != 0, and
// (3) y is an integer if x < -1.
template <typename T>
bool pow1p_special(const T &x, const T &y, T &result) {
    result = 0;

    // Handle y == 0 or (1+x) == 1.  IEEE-754 requires the result to be 1
    // even if the other operand is nan.  Therefore this test must precede
    // the test for nan.
    if (y == 0 || x == 0) {
        result = 1;
        return true;
    }

    // Handle nan.
    if (std::isnan(x)) {
        result = x;
        return true;
    }
    if (std::isnan(y)) {
        result = y;
        return true;
    }

    // Handle y == +/-inf.  This test must precede the test for x == +/-inf.
    if (y == std::numeric_limits<T>::infinity()) {
        result = (x == -2) ? T(1) :
                 (x < 0 && x > -2) ? T(0) :
                 std::numeric_limits<T>::infinity();
        return true;
    }
    if (y == -std::numeric_limits<T>::infinity()) {
        result = (x == -2) ? T(1) :
                 (x < 0 && x > -2) ? std::numeric_limits<T>::infinity() :
                 T(0);
        return true;
    }

    // Handle 1+x == +/- inf.
    if (x == std::numeric_limits<T>::infinity()) {
        result = (y < 0) ? T(0) : std::numeric_limits<T>::infinity();
        return true;
    }
    if (x == -std::numeric_limits<T>::infinity()) {
        bool y_is_odd_int = (std::abs(std::fmod(y, T(2))) == 1);
        if (y_is_odd_int) {
            result = (y < 0)? -T(0) : -std::numeric_limits<T>::infinity();
        } else {
            result = (y < 0) ? T(0) : std::numeric_limits<T>::infinity();
        }
        return true;
    }

    // Handle 1+x == +/-0 (can only be +0 in fact).
    if (x == -1) {
        if (y < 0) {
            // TODO: signal divideByZero exception
            result = std::numeric_limits<T>::infinity();
        } else {
            // y == 0 is already handled above
            result = 0;
        }
        return true;
    }

    // Handle 1+x < 0 and y not an integer
    if ((x < -1) && std::fmod(y, T(1)) != 0) {
        // TODO: signal invalid operation exception
        result = std::numeric_limits<T>::quiet_NaN();
        return true;
    }

    // Not a special input.
    return false;
}

// Compute pow(x1+x2,y) by decomposing x1+x2 == s+t exactly such that |t| << s
// and sign(t)*sign(s-1) >= 0.  The caller must ensure that the decomposition
// is actually possible.
template <typename T>
T pow1p_decomp(const T &x1, const T &x2, const T &y) {
    auto out = stderr;
    fprintf(out, "[pow1p] ========\n");
    fprintf(out, "[pow1p] x1=%.15e x2=%.15e y=%.15e\n", x1, x2, y);

    // Write (x+dx) == (s+t) where s is equal to (x+dx) rounded toward 1
    // and t is the (exact) rounding error.
    T s, t;
    two_sum(x1, x2, s, t);
    if ((t > 0 && s < 1) || (t < 0 && s > 1)) {
        T delta = s - std::nextafter(s, T(1)); // exact
        s -= delta; // exact
        t += delta; // exact under precondition
    }
    fprintf(out, "[pow1p] s=%.15e t=%.15e\n", s, t);

    // Write (s+t)^y == (s^y)*((1+t/s)^y) == term1*term2.  Since s is rounded
    // toward one, both terms <= 1 or both >= 1.  So if term1 over/underflows,
    // then term1*term2 necessarily over/underflows.
    // And of course, term2 == 1 if t == 0.
    T term1 = std::pow(s, y);
    fprintf(out, "[pow1p] term1=%.15e\n", term1);

    if (t == 0 || term1 == 0 || std::isinf(term1)) {
        // TODO: handle floating point exceptions
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
    fprintf(out, "[pow1p] term2=%.15e\n", term2);
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
    T vv = std::fma(-0.5 * u, u, uu);

    // (w + ww) ~= y*(u+vv) ~= y*log(1+t/s).
    T r2 = std::fma(y, u, -w);
    T ww = std::fma(y, vv, r2);

    // TODO: maybe ww is small enough such that exp(ww) ~= 1+ww.
    T term3 = std::exp(ww);
    fprintf(out, "[pow1p] term3=%.15e\n", term3);
    return term1 * term2 * term3;
}

// Compute pow(1+x,y) to full precision of T.
template <typename T>
T pow1p_precise(T x, T y) {
    auto out = stderr;
    fprintf(out, "[pow1p] ========\n");
    fprintf(out, "[pow1p] x=%.16e y=%.16e\n", x, y);

    // Handle special values.
    T answer;
    if (pow1p_special(x, y, answer)) {
        return answer;
    }

    // Now (1) x and y are finite, (2) x != 0, x != -1, y != 0, and
    // (3) y is an integer if x < -1.

    // Handle (1+x) < 0 (and integer y).
    if (x < -1) {
        // If 1 < (-x) <= 2^p, ((-x)-1) is exact; if 2^p < (-x) <= 2^(2p),
        // ((-x)-1) can be written as (s+t) exactly where t << s and t > 0.
        // For both cases we delegate to pow1p_decomp().
        //
        // If (-x) > 2^(2p), ((-x)-1) ~= (-x) and the total error is small.
        // We call pow() directly in this case.
        constexpr T exp2p = T(2) / std::numeric_limits<T>::epsilon();
        if (-x <= exp2p*exp2p) {
            T sign = (std::fmod(y, T(2)) == 0) ? 1 : -1;
            return sign * pow1p_decomp(-x, T(-1), y);
        } else {
            // TODO: handle floating point exceptions
            return std::pow(T(1) + x, y);
        }
    }

    // Now x, y are finite and (1+x) > 0.
    return pow1p_decomp(x, T(1), y);
}

// Compute pow(1+x,y) using exp(y*log1p(x)) or similar.
template <typename T>
T pow1p_approx(T x, T y) {
    auto &out = std::cerr;
    out.precision(std::numeric_limits<T>::digits10 - 1);
    out.setf(std::ios::scientific, std::ios::floatfield);

    // Note: different T has different denorm and overflow range.
    out << "[pow1p] ========" << std::endl;
    out << "[pow1p] x=" << x << " y=" << y << std::endl;

    // Handle special values.
    T answer;
    if (pow1p_special(x, y, answer)) {
        return answer;
    }

    if (x > -1) {
        T z = std::log1p(x);
        out << "[pow1] z=" << z << std::endl;
        T yz = y * z;
        out << "[pow1] y*z=" << yz << std::endl;
        return std::exp(yz);
    } else {
        // x < -1
        return std::pow(x+1, y);
    }
}

// Return the extra number of bits needed for using the quick method.
template <class T>
constexpr int extra_bits_needed() {
    constexpr int p = std::numeric_limits<T>::digits;
    constexpr int emin = std::numeric_limits<T>::min_exponent;
    constexpr int emax = std::numeric_limits<T>::max_exponent;
    constexpr int m = (std::max)(p - emin, emax + 1);
    return 1 + std::ilogb(1 + 1.05*m); // only constexpr since C++23
}

inline double pow1p(double x, double y) {
    /*constexpr*/ int q = extra_bits_needed<double>();
    if /*constexpr*/(std::numeric_limits<long double>::digits -
                 std::numeric_limits<double>::digits >= q) {
        return pow1p_approx<long double>(x, y);
    } else {
        return pow1p_precise(x, y);
    }
}

inline float pow1p(float x, float y) {
    /*constexpr*/ int q = extra_bits_needed<float>();
    if /*constexpr*/(std::numeric_limits<double>::digits -
                 std::numeric_limits<float>::digits >= q) {
        return pow1p_approx<double>(x, y);
    } else if /*constexpr*/(std::numeric_limits<long double>::digits -
                        std::numeric_limits<float>::digits >= q) {
        return pow1p_approx<long double>(x, y);
    } else {
        return pow1p_precise(x, y);
    }
}

} // namespace xsf
