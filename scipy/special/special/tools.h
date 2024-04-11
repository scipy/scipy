/* Building blocks for implementing special functions */

#pragma once

#include "config.h"
#include "error.h"
#include <utility>
#include <cstddef>

namespace special {
namespace detail {

    // Return NaN, handling both real and complex types.
    template <typename T>
    SPECFUN_HOST_DEVICE inline std::enable_if_t<std::is_floating_point_v<T>, T> maybe_complex_NaN() {
        return std::numeric_limits<T>::quiet_NaN();
    }

    template <typename T>
    SPECFUN_HOST_DEVICE inline std::enable_if_t<!std::is_floating_point_v<T>, T> maybe_complex_NaN() {
        using V = typename T::value_type;
        return {std::numeric_limits<V>::quiet_NaN(), std::numeric_limits<V>::quiet_NaN()};
    }

    // Series evaluators.
    template <typename Generator, typename T>
    SPECFUN_HOST_DEVICE T series_eval(Generator &g, T init_val, double tol, std::uint64_t max_terms,
                                           const char *func_name) {
        /* Sum an infinite series to a given precision.
         *
         * g : a generator of terms for the series.
         *
         * init_val : A starting value that terms are added to. This argument determines the
         *     type of the result.
         *
         * tol : relative tolerance for stopping criterion.
         *
         * max_terms : The maximum number of terms to add before giving up and declaring
         *     non-convergence.
         *
         * func_name : The name of the function within SciPy where this call to series_eval
         *     will ultimately be used. This is needed to pass to set_error in case
         *     of non-convergence.
         */
        T result = init_val;
        T term, previous;
        for (std::uint64_t i = 0; i < max_terms; ++i) {
            term = g();
            result += term;
            if (std::abs(term) < std::abs(result) * tol) {
                return result;
            }
        }
        // Exceeded max terms without converging. Return NaN.
        set_error(func_name, SF_ERROR_NO_RESULT, NULL);
        return maybe_complex_NaN<T>();
    }

    template <typename Generator, typename T>
    SPECFUN_HOST_DEVICE T series_eval_fixed_length(Generator &g, T init_val, std::uint64_t num_terms) {
        /* Sum a fixed number of terms from a series.
         *
         * g : a generator of terms for the series.
         *
         * init_val : A starting value that terms are added to. This argument determines the
         *     type of the result.
         *
         * max_terms : The number of terms from the series to sum.
         *
         */
        T result = init_val;
        for (std::uint64_t i = 0; i < num_terms; ++i) {
            result += g();
        }
        return result;
    }

    namespace maybe_complex_numeric_limits {
        // Handle numeric limits when type may be complex.
        template <typename T>
        SPECFUN_HOST_DEVICE inline std::enable_if_t<std::is_floating_point_v<T>, T> quiet_NaN() {
            return std::numeric_limits<T>::quiet_NaN();
        }

        template <typename T>
        SPECFUN_HOST_DEVICE inline std::enable_if_t<!std::is_floating_point_v<T>, T> quiet_NaN() {
            using V = typename T::value_type;
            return {std::numeric_limits<V>::quiet_NaN(), std::numeric_limits<V>::quiet_NaN()};
        }

        template <typename T>
        SPECFUN_HOST_DEVICE inline std::enable_if_t<std::is_floating_point_v<T>, T> min() {
            return std::numeric_limits<T>::min();
        }

        template <typename T>
        SPECFUN_HOST_DEVICE inline std::enable_if_t<!std::is_floating_point_v<T>, T> min() {
            using V = typename T::value_type;
            return std::numeric_limits<V>::min();
        }

    } // namespace maybe_complex_numeric_limits

    template <typename Generator, typename T>
    SPECFUN_HOST_DEVICE T continued_fraction_eval(Generator &g, T init_val, double tol, std::uint64_t max_terms,
                                                  const char *func_name) {
        /* Evaluate continued fraction with modified Lentz's algorithm.
         *
         * Evaluates the continued fraction b0 + a1 / (b1 + a1 / (b2 + a2 / ( ...
         *
         * g : Generator of pairs of values (a1, b1), (a2, b2), (a3, b3), ...
         *
         * init_val : Initial value b0
         *
         * tol : relative tolerance for stopping criterion.
         *
         * max_terms : The maximum number of iterations before giving up and declaring
         *     non-convergence.
         *
         * func_name : The name of the function within SciPy where this call to series_eval
         *     will ultimately be used. This is needed to pass to set_error in case
         *     of non-convergence.
         */
        std::pair<T, T> v;

        T tiny_value = 16.0 * maybe_complex_numeric_limits::min<T>();
        T f = (init_val == 0.0) ? tiny_value : init_val;

        double C = f;
        double D = 0.0;
        double delta;

        for (uint64_t i = 0; i < max_terms; i++) {
            v = g();
            D = v.second + v.first * D;
            if (D == 0.0) {
                D = tiny_value;
            }
            C = v.second + v.first / C;
            if (C == 0.0) {
                C = tiny_value;
            }
            D = 1.0 / D;
            delta = C * D;
            f *= delta;
            if (std::abs(delta - 1.0) <= tol) {
                return f;
            }
        }
        // Exceeded max terms without converging. Return NaN.
        set_error(func_name, SF_ERROR_NO_RESULT, NULL);
        return maybe_complex_numeric_limits::quiet_NaN<T>();
    }

    /* Evaluates a continued fraction using the modified Lentz algorithm.
     *
     * Parameters
     * ----------
     *
     *   init_val
     *       b0.  The type of this parameter (T) is used for intermediary
     *       computations as well as the return value.  It is expected to
     *       be a (possibly complex) floating point type.
     *
     *   cf
     *       Input iterator that yields the terms of the continued fraction.
     *       It must support two operations: `*cf` returns the current term
     *       (a[n], b[n]) as a pair; `++cf` advances to the next term. The
     *       iterator initially points to (a[1], b[1]).
     *
     *   tol
     *       Tolerance used to terminate the iteration.  Specifically, stop
     *       iteration as soon as ...
     *
     *   max_terms
     *       Maximum number of terms to evaluate.  This parameter is used to
     *       handle unexpected non-convergence.  It should normally be set
     *       large enough such that the continued fraction is mathematically
     *       guaranteed to have converged within that many terms.
     *
     *   tiny_val
     *       If any intermediary C or D is exactly zero, it is replaced by
     *       `tiny_val`.  It should be set to a "small" number; the default
     *       choice of boost is 16 times the smallest positive normalized
     *       value representable by T.  Setting `tiny_val` to zero disables
     *       this workaround.
     *
     * Return Value
     * ------------
     *
     * If the continued fraction converges (i.e. the tolerance condition is
     * satisfied) by f[n] with `n <= max_terms`, returns (f[n], n).
     * Otherwise, returns (f[max_terms], 0).
     *
     * Remarks
     * -------
     *
     * This is a low-level routine intended for expert use.  No error checking
     * is performed.  The caller must ensure that the parameters and terms are
     * finite and that intermediary calculations do not trigger floating point
     * exceptions such as overflow.
     *
     * The numerical accuracy of this method depends on the characteristics of
     * the continued fraction being evaluated.
     */
    template <typename T, typename InputIt>
    SPECFUN_HOST_DEVICE std::pair<T, std::size_t> continued_fraction_eval_lentz(
        T init_val, InputIt cf, T tol, std::size_t max_terms, T tiny_val) {

        T f = (init_val == 0)? tiny_val : init_val;
        T C = f;
        T D = 0;
        for (std::size_t n = 0; n < max_terms; ++n, ++cf) {
            T a = std::get<0>(*cf);  // current numerator
            T b = std::get<1>(*cf);  // current denominator

            D = b + a * D;
            if (D == 0) {
                D = tiny_val;
            }
            C = b + a / C;
            if (C == 0) {
                C = tiny_val;
            }
            D = 1 / D;
            T delta = C * D;
            f *= delta;
            if (abs(delta - 1) <= abs(tol)) {
                return {f, n+1};
            }
        }

        // Failed to converge within max_terms terms.
        return {f, 0};
    }

    /* Performs one step of Kahan summation. */
    template <typename T>
    void kahan_step(T *sum, T *comp, T x) {
        T y = x - *comp;
        T t = *sum + y;
        *comp = (t - *sum) - y;
        *sum = t;
    }

    /* Evaluates a continued fraction using the "series" method.
     *
     * Denote the continued fraction by
     *
     *              a[1]   a[2]   a[3]
     *   f = b[0] + ------ ------ ------ ...
     *              b[1] + b[2] + b[3] +
     *
     * Denote the n-th convergent by f[n], starting from f[0] := b[0].  The
     * series method (see [1])
     *
     * Reference
     * ---------
     *
     * [1] Gautschi, W. (1967). “Computational Aspects of Three-Term
     *     Recurrence Relations.” SIAM Review, 9(1):24-82.
     *
     * Parameters
     * ----------
     *
     *   init_val
     *       b0.  The type of this parameter (T) is used for intermediary
     *       computations as well as the return value.  It is expected to
     *       be a (possibly complex) floating point type.
     *
     *   cf
     *       Input iterator that yields the terms of the continued fraction.
     *       It must support two operations: `*cf` returns the current term
     *       (a[n], b[n]) as a pair; `++cf` advances to the next term. The
     *       iterator initially points to (a[1], b[1]).
     *
     *   tol
     *       Tolerance used to terminate the iteration.  Specifically, stop
     *       iteration as soon as `abs(f[n] - f[n-1]) <= abs(tol * f[n])`,
     *       where f[n], n >= 1 is the n-th convergent, and f[0] = b0.  The
     *       value of `tol` should be chosen to bound the truncation error.
     *
     *   max_terms
     *       Maximum number of terms to evaluate.  This parameter is used to
     *       handle unexpected non-convergence.  It should normally be set
     *       large enough such that the continued fraction is mathematically
     *       guaranteed to have converged within that many terms.
     *
     * Return Value
     * ------------
     *
     * If the continued fraction converges (i.e. the tolerance condition is
     * satisfied) by f[n] with `n <= max_terms`, returns (f[n], n).
     * Otherwise, returns (f[max_terms], 0).
     *
     * Remarks
     * -------
     *
     * This is a low-level routine intended for expert use.  No error checking
     * is performed.  The caller must ensure that the parameters and terms are
     * finite and that intermediary calculations do not trigger floating point
     * exceptions such as overflow.
     *
     * The numerical accuracy of this method depends on the characteristics of
     * the continued fraction being evaluated.
     */
    template <typename T, typename InputIt>
    SPECFUN_HOST_DEVICE std::pair<T, std::size_t> continued_fraction_eval_series(
        T init_val, InputIt cf, T tol, std::size_t max_terms) {

        T w = init_val;  // current convergent
        T c = 0;         // compensation for Kahan summation
        if (max_terms == 0) {
            return {w, 0};
        }

        // n = 1
        T a = std::get<0>(*cf);  // a == current numerator
        T b = std::get<1>(*cf);  // b == current denominator
        T v = a / b;             // v == f[n] - f[n-1]
        kahan_step(&w, &c, v);   // w == f[n]
        if (abs(v) <= abs(tol * w)) {
            return {w, 1};
        }

        // n = 2, 3, 4, ...
        T u = 1;
        T bb = b;  // last denominator
        for (std::size_t n = 1; n < max_terms; ++n) {
            ++cf;
            a = std::get<0>(*cf);
            b = std::get<1>(*cf);
            u = 1 / (1 + a / (b * bb) * u);
            v *= (u - 1);
            kahan_step(&w, &c, v);
            if (abs(v) <= abs(tol * w)) {
                return {w, n + 1u};
            }
            bb = b;
        }

        // Failed to converge within max_terms terms.
        return {w, 0};
    }

} // namespace detail
} // namespace special
