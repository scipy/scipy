#pragma once

#include "error.h"
#include "specfun.h"

namespace special {

template <typename T>
T exp1(T x) {
    T out = specfun::e1xb(x);
    SPECFUN_CONVINF("exp1", out);
    return out;
}

template <typename T>
std::complex<T> exp1(std::complex<T> z) {
    std::complex<T> outz = specfun::e1z(z);
    SPECFUN_ZCONVINF("exp1", outz);
    return outz;
}

template <typename T>
T expi(T x) {
    T out = specfun::eix(x);
    SPECFUN_CONVINF("expi", out);
    return out;
}

template <typename T>
std::complex<T> expi(std::complex<T> z) {
    std::complex<T> outz = specfun::eixz(z);
    SPECFUN_ZCONVINF("cexpi", outz);
    return outz;
}

namespace detail {

    //
    // Compute a factor of the exponential integral E1.
    // This is used in scaled_exp1(x) for moderate values of x.
    //
    // The function uses the continued fraction expansion given in equation 5.1.22
    // of Abramowitz & Stegun, "Handbook of Mathematical Functions".
    // For n=1, this is
    //
    //    E1(x) = exp(-x)*C(x)
    //
    // where C(x), expressed in the notation used in A&S, is the continued fraction
    //
    //            1    1    1    2    2    3    3
    //    C(x) = ---  ---  ---  ---  ---  ---  ---  ...
    //           x +  1 +  x +  1 +  x +  1 +  x +
    //
    // Here, we pull a factor of 1/z out of C(x), so
    //
    //    E1(x) = (exp(-x)/x)*F(x)
    //
    // and a bit of algebra gives the continued fraction expansion of F(x) to be
    //
    //            1    1    1    2    2    3    3
    //    F(x) = ---  ---  ---  ---  ---  ---  ---  ...
    //           1 +  x +  1 +  x +  1 +  x +  1 +
    //
    inline double expint1_factor_cont_frac(double x) {
        // The number of terms to use in the truncated continued fraction
        // depends on x.  Larger values of x require fewer terms.
        int m = 20 + (int) (80.0 / x);
        double t0 = 0.0;
        for (int k = m; k > 0; --k) {
            t0 = k / (x + k / (1 + t0));
        }
        return 1 / (1 + t0);
    }

} // namespace detail

//
// Scaled version  of the exponential integral E_1(x).
//
// Factor E_1(x) as
//
//    E_1(x) = exp(-x)/x * F(x)
//
// This function computes F(x).
//
// F(x) has the properties:
//  * F(0) = 0
//  * F is increasing on [0, inf)
//  * lim_{x->inf} F(x) = 1.
//
inline double scaled_exp1(double x) {
    if (x < 0) {
        return NAN;
    }

    if (x == 0) {
        return 0.0;
    }

    if (x <= 1) {
        // For small x, the naive implementation is sufficiently accurate.
        return x * std::exp(x) * exp1(x);
    }

    if (x <= 1250) {
        // For moderate x, use the continued fraction expansion.
        return detail::expint1_factor_cont_frac(x);
    }

    // For large x, use the asymptotic expansion.  This is equation 5.1.51
    // from Abramowitz & Stegun, "Handbook of Mathematical Functions".
    return 1 + (-1 + (2 + (-6 + (24 - 120 / x) / x) / x) / x) / x;
}

inline float scaled_exp1(float x) { return scaled_exp1(static_cast<double>(x)); }

} // namespace special
