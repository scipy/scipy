/* Building blocks for implementing special functions */

#pragma once

#include "config.h"
#include "error.h"
#include <utility>
#include <cstddef>
#include <iterator>

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

    template <typename T, typename Iterator>
    std::pair<T, std::size_t> continued_fraction_eval_lentz(
        T init_val, Iterator cf, T tol, std::size_t max_terms) {

        std::pair<T, T> v;

        T tiny_value = 16.0 * maybe_complex_numeric_limits::min<T>();
        T f = (init_val == 0.0) ? tiny_value : init_val;

        double C = f;
        double D = 0.0;
        double delta;

        for (std::size_t i = 0; i < max_terms; i++) {
            v = *(cf++);
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
                return {f, i+1};
            }
        }

        // Exceeded max terms without converging.
        return {f, 0};
    }

    /* Evaluates a continued fraction using the "series" method.
     *
     * Denote the continued fraction by
     *
     *                a[1]     a[2]     a[3]
     *   f = b[0] + -------- -------- -------- ...
     *               b[1] +   b[2] +   b[3] +
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
     *   initial_value
     *       b0.  The type of this parameter is used as the "accumulator"
     *       as well as the result.  It may be different from the type of
     *       the terms of the continued fraction.  For example, a higher
     *       precision data type may be supplied to improve the accuracy.
     *
     *   cf
     *       Input iterator that yields the terms of the c.f.  The value_type
     *       of the iterator must be std::pair<T,T>, where T is the data type
     *       of the numerator and denominator.  The current term is given by
     *       ((*cf).first, (*cf).second).  The iterator initially points to
     *       (a1, b1), and advances to the next term each time `++cf` is
     *       called.
     *
     *   tol
     *       Tolerance used to terminate the iteration.  Specifically, stop
     *       iteration as soon as `abs(f[n] - f[n-1]) <= tol * abs(f[n])`,
     *       where f[n], n >= 1 is the n-th convergent, and f[0] = b0.
     *
     *   max_terms
     *       Maximum number of terms to evaluate.  Should be positive.
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
     * This is a low-level routine for internal use.  No error check of any
     * kind is performed.  The caller must ensure that all parameters and
     * terms are finite, and that all intermediary calculations are valid
     * (i.e. do not lead to overflow or nan).  `tol` should be suitably
     * chosen to bound the truncation error.  `max_terms` should be set
     * large enough and is meant to protect against unexpected errors.
     *
     * The numerical stability of this method depends on the characteristics
     * of the continued fraction being evaluated.  For example, it works well
     * if `f[n] - f[n-1]` are all positive and geometrically decreasing.
     */
    template <typename Result, typename Iterator, typename Tolerance>
    SPECFUN_HOST_DEVICE std::pair<Result, std::size_t> continued_fraction_eval_series(
        Result initial_value, Iterator cf, Tolerance tol, std::size_t max_terms) {

        // T is type of numerator and denominator.
        using T = typename std::iterator_traits<Iterator>::value_type::first_type;

        Result w = initial_value;  // current convergent
        if (max_terms == 0) {
            return {w, 0};
        }

        // n = 1
        std::pair<T, T> ab = *cf;  // ab == current fraction
        T a = ab.first;            // a == current numerator
        T b = ab.second;           // b == current denominator
        T v = a / b;               // v == f[n] - f[n-1]
        w += v;                    // w == f[n]
        if (std::abs(v) <= tol * std::abs(w)) {
            return {w, 1};
        }

        // n = 2, 3, 4, ...
        T u = 1;
        T bb = b;  // last denominator
        for (std::size_t n = 1; n < max_terms; ++n) {
            ab = *(++cf);
            a = ab.first;
            b = ab.second;
            u = 1 / (1 + a / (b * bb) * u);
            v *= (u - 1);
            w += v;
            if (std::abs(v) <= tol * std::abs(w)) {
                return {w, n + 1u};
            }
            bb = b;
        }

        // Failed to converge within max_terms terms.
        return {w, 0};
    }

    /* Performs Kahan summation. */
    template <typename T>
    class KahanSummer {

    public:
        KahanSummer() : _s(), _c() {} // value initialization

        KahanSummer & operator+=(T x) {
            T y = x - _c;
            T t = _s + y;
            _c = (t - _s) - y;
            _s = t;
            return *this;
        }

        operator T() const { return _s; }

    private:
        T _s; // current sum
        T _c; // current compensation
    };

} // namespace detail
} // namespace special
