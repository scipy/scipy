#include <cmath>
#include <cstdint>
#include <limits>
#include <xsf/zeta.h>  // for xsf::cephes::zeta(a, n)
#include "sf_error.h"


//
// If zeta(a, n + 1) / zeta(a, 1) exceeds zeta_ratio_threshold, the two
// values are close enough that the loss of precision in the subtraction
// zeta(a, 1) - zeta(a, n + 1) should be avoided by using the direct
// sum of powers instead.  The most conservative value would 0.5, but
// experimentation shows that 0.9 maintains a relative error of less
// than 5e-15.
//
static const double zeta_ratio_threshold = 0.9;

//
// *** NOTE ***
// In the following templated functions, T will be either long long int
// or double.  When T is double, the *values* of the variables with that
// type are still expected to be positive integers less than 2**53.
// This is necessary so that, for example, incrementing such a variable with
// the ++ operator will, in fact, change the value by 1.
// The bounds and integrality of the inputs are not checked in the functions;
// it is expected that the Python code that ends up calling the ufuncs will
// ensure that these conditions are met.
//

//
// Compute sum_{i=m}^{n} i**-a.
//
// This function assumes 1 <= m <= n and `a` is finite.
//
// See the NOTE above regarding T.
//
template <typename T>
static inline double
sum_powers(T m, T n, double a)
{
    double sum = 0.0;
    if (a >= 0) {
        for (T i = n; i >= m; --i) {
            sum += std::pow(i, -a);
        }
    }
    else {
        for (T i = m; i <= n; ++i) {
            sum += std::pow(i, -a);
        }
    }
    return sum;
}

//
// Compute
//
//   sum_{i=1}^{n} i**-a
//
// This is the generalized harmonic number H_{n, a} [1].
//
// * If n < 1, NAN is returned with SF_ERROR_DOMAIN set.
// * If a > 1 and not too close to 1, the sum is computed using
//   the Hurwitz zeta function.  Otherwise it is computed as a simple
//   sum of the powers.
//
// See the NOTE above regarding T.
//
// [1] https://en.wikipedia.org/wiki/Harmonic_number#Generalized_harmonic_numbers
//
template <typename T>
static inline double
gen_harmonic(T n, double a)
{
    if constexpr(std::is_same_v<T, double>) {
        if (std::isnan(n)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    if (n < 1) {
        sf_error("_gen_harmonic", SF_ERROR_DOMAIN,
                 "n >= 1 is required, but got n = %" PRId64, static_cast<int64_t>(n));
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (n == 1) {
        // IEEE: pow(1.0, _) is 1.0.
        return 1.0;
    }
    if (std::isnan(a)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (std::isinf(a)) {
        if (a > 0) {
            // a = +inf
            return 1.0;
        }
        else {
            // a = -inf
            return std::numeric_limits<double>::infinity();
        }
    }
    if (a == 0) {
        return static_cast<double>(n);
    }
    if (a <= 1) {
        return sum_powers(static_cast<T>(1), n, a);
    }
    else {
        // If here, we know a is finite and a > 1, and n > 1.
        //
        // For a > 1, we can use the formula zeta(a, 1) - zeta(a, n + 1),
        // where zeta(a, k) is the Hurwitiz zeta function.
        // But if zeta(a, 1) and zeta(a, n + 1) are close, precision is lost
        // in the subtraction, so we use the explicit sum instead. We consider
        // the values "close" if zeta(a, n + 1)/zeta(a, 1) > zeta_ratio_threshold.
        //
        double z1 = xsf::cephes::zeta(a, static_cast<double>(1));
        double znp1 =  xsf::cephes::zeta(a, static_cast<double>(n + 1));
        if (znp1 / z1 <= zeta_ratio_threshold) {
            return z1 - znp1;
        }
        else {
            return sum_powers(static_cast<T>(1), n, a);
        }
    }
}


//
// Computes
//
//      sum_{i=j}^{k} i**-a
//      -------------------
//      sum_{i=1}^{n} i**-a
//
// Requires 1 <= j <= k <= n and finite a.
//
// See the NOTE above regarding T.
//

template <typename T>
static inline double
normalized_sum_powers(T j, T k, T n, double a)
{
    double numer = 0.0;
    double denom = 0.0;
    if (a >= 0) {
        for (T i = n; i >= 1; --i) {
            double term = std::pow(i, -a);
            denom += term;
            if (i >= j && i <= k) {
                numer += term;
            }
        }
    }
    else {
        for (T i = 1; i <= n; ++i) {
            double term = std::pow(i, -a);
            denom += term;
            if (i >= j && i <= k) {
                numer += term;
            }
        }
    }
    return numer / denom;
}

//
// Computes
//
//      sum_{i=j}^{k} i**-a
//      -------------------
//      sum_{i=1}^{n} i**-a
//
// Requires 1 <= j <= k <= n; returns NAN if that is not true.
//
// See the NOTE above regarding T.
//
template <typename T>
static inline double
normalized_gen_harmonic(T j, T k, T n, double a)
{
    if constexpr(std::is_same_v<T, double>) {
        if (std::isnan(j) || std::isnan(k) || std::isnan(n)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    if (j < 1 || k < j || n < k) {
        sf_error("_normalized_gen_harmonic", SF_ERROR_DOMAIN,
                 "1 <= j <= k <= n is required, but got j = %" PRId64 ", "
                 "k = %" PRId64 ", and n = %" PRId64,
                 static_cast<int64_t>(j),
                 static_cast<int64_t>(k),
                 static_cast<int64_t>(n));
        return std::numeric_limits<double>::quiet_NaN();
    }
    //
    // Now we know 1 <= j <= k <= n
    //
    if (n == 1) {
        // IEEE: pow(1.0, _) is 1.0.
        // n == 1 implies j == k == 1.
        return 1.0;
    }
    if (std::isnan(a)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (std::isinf(a)) {
        if (a > 0) {
            // a = +inf
            if (j == 1) {
                // Numerator and denominator are both 1.
                return 1.0;
            }
            else {
                // Numerator is 0, denominator is 1.
                return 0.0;
            }
        }
        else {
            // a = -inf
            if (k == 1) {
                // Numerator is 1, denominator is +inf.
                return 0.0;
            }
            else {
                // Numerator and denominator are both +inf.
                return std::numeric_limits<double>::quiet_NaN();
            }
        }
    }
    if (a == 0) {
        return static_cast<double>(k - j + 1) / static_cast<double>(n);
    }
    if (a <= 1) {
        return normalized_sum_powers(j, k, n, a);
    }
    else {
        // If here, we know a is finite and a > 1, n > 1, and 1 <= j <= k <= n.
        // See the comments in _gen_harmonic() for an explanation of why
        // we check the ratios of the zeta functions before using them
        // to compute the result.
        double zj = xsf::cephes::zeta(a, static_cast<double>(j));
        double zkp1 =  xsf::cephes::zeta(a, static_cast<double>(k + 1));
        double z1  = xsf::cephes::zeta(a, static_cast<double>(1));
        double znp1 = xsf::cephes::zeta(a, static_cast<double>(n + 1));
        bool zeta_numer_ok = zkp1 / zj < zeta_ratio_threshold;
        bool zeta_denom_ok = znp1 / z1 < zeta_ratio_threshold;
        if (zeta_numer_ok) {
            if (zeta_denom_ok) {
                // OK to use the zeta formula for the numerator
                // and denominator.
                return (zj - zkp1) / (z1 - znp1);
            }
            else {
                return (zj - zkp1) / sum_powers(static_cast<T>(1), n, a);
            }
        }
        else {
            if (zeta_denom_ok) {
                return sum_powers(j, k, a) / (z1 - znp1);
            }
            else {
                return normalized_sum_powers(j, k, n, a);
            }
        }
    }
}
