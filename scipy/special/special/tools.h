/* Building blocks for implementing special functions */

#pragma once

#include "config.h"
#include "error.h"
#include <utility>
#include <cstddef>

namespace special {
namespace detail {

    template <typename T>
    struct real_type {
        using type = T;
    };

    template <typename T>
    struct real_type<std::complex<T>> {
        using type = T;
    };

    template <typename T>
    using real_type_t = typename real_type<T>::type;

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

    /* Performs one step of Kahan summation. */
    template <typename T>
    void kahan_step(T *sum, T *comp, T x) {
        T y = x - *comp;
        T t = *sum + y;
        *comp = (t - *sum) - y;
        *sum = t;
    }

    /* Evaluates an infinite series using Kahan summation.
     *
     * Denote the series by
     *
     *   S = a[0] + a[1] + a[2] + ...
     *
     * For n >= 0, denote its n-th partial sum by
     *
     *   S[n] = a[0] + a[1] + ... + a[n]
     *
     * This function computes S[0], S[1], ... until a[n] is sufficiently
     * small or if the maximum number of terms have been evaluated.
     *
     * Parameters
     * ----------
     *
     *   init_val
     *       a[0].  The type of this argument is used for intermediary
     *       computations as well as the return value.
     *
     *   g
     *       Input iterator that yields the terms a[1], a[2], a[3], ...
     *
     *   tol
     *       Relative tolerance for convergence.  Specifically, stop iteration
     *       as soon as `abs(a[n]) <= tol * abs(S[n])` for some n >= 1.
     *
     *   max_terms
     *       Maximum number of terms after a[0] to evaluate.
     *
     * Return Value
     * ------------
     *
     * If the convergence criterion is satisfied by some n <= max_terms,
     * returns `(S[n], n)`.  Otherwise, returns `(S[max_terms], 0)`.
     */
    template <typename T, typename InputIt>
    SPECFUN_HOST_DEVICE std::pair<T, std::size_t> series_eval_kahan(
        T init_val, InputIt g, real_type_t<T> tol, std::size_t max_terms) {

        T sum = init_val;
        T comp = 0;
        for (std::size_t i = 0; i < max_terms; ++i, ++g) {
            T term = *g;
            kahan_step(&sum, &comp, term);
            if (std::abs(term) <= tol * std::abs(sum)) {
                return {sum, i + 1};
            }
        }
        return {sum, 0};
    }

#if 0
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
#endif

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

    /* Forward iterator to generate the difference of successive convergents
     * of a continued fraction.
     *
     * Let f denote the continued fraction, and let f[n], n >= 0 denote its
     * n-th convergent.  Write
     *
     *              a[1]   a[2]   a[3]
     *   f = b[0] + ------ ------ ------ ...
     *              b[1] + b[2] + b[3] +
     *
     *                 a[1]   a[2]       a[n]
     *   f[n] = b[0] + ------ ------ ... ----
     *                 b[1] + b[2] +     b[n]
     *
     * with f[0] = b[0].  This iterator generates the sequence of values
     * f[1]-f[0], f[2]-f[1], f[3]-f[2], ...
     *
     * Type Arguments
     * --------------
     *
     *   T
     *       Type in which computations are performed and results are turned.
     *
     *   FwdIt
     *       Forward iterator that yields the terms of the continued fraction
     *       as (numerator, denominator) pairs, starting from (a[1], b[1]).
     *       It must support two operations: `*cf` to return the current term
     *       (a[n], b[n]) as a pair, and `++cf` to advance to the next term.
     *
     *       IMPORTANT: For performance reason, the generator always eagerly
     *       dereference the current term of the continued fraction.  That is,
     *       (a[1], b[1]) is dereferenced upon generator construction, and
     *       (a[n], b[n]) is dereferenced after (n-1) calls of `++`.  FwdIt
     *       must ensure it is always dereference-able.
     *
     * Remarks
     * -------
     *
     * The values are computed using the recurrence relation described in [1].
     *
     * No error checking is performed.  The caller must ensure that all terms
     * are finite and that intermediary computations do not trigger floating
     * point exceptions such as overflow.
     *
     * The numerical stability of this method depends on the characteristics
     * of the continued fraction being evaluated.
     *
     * Reference
     * ---------
     *
     * [1] Gautschi, W. (1967). “Computational Aspects of Three-Term
     *     Recurrence Relations.” SIAM Review, 9(1):24-82.
     */
    template <typename T, typename FwdIt>
    class ContinuedFractionDifferenceGenerator {
        using Self = ContinuedFractionDifferenceGenerator<T, FwdIt>;

    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = T;
        using difference_type = std::ptrdiff_t;
        using pointer = const T *;
        using reference = const T &;

        ContinuedFractionDifferenceGenerator(FwdIt cf) : _cf(cf) { _init(); }

        reference operator*() const { return _v; }

        Self& operator++() {
            _advance();
            return *this;
        }

        Self operator++(int) {
            Self self = *this;
            _advance();
            return self;
        }

    private:
        void _init() {
            auto [num, denom] = *_cf;
            T a = num;
            T b = denom;
            _u = 1;
            _v = a / b;
            _b = b;
        }

        void _advance() {
            auto [num, denom] = *(++_cf);
            T a = num;
            T b = denom;
            _u = 1 / (1 + (a * _u) / (b * _b));
            _v *= (_u - 1);
            _b = b;
        }

        FwdIt _cf;  // points to current fraction
        T _v;       // v[n] == f[n] - f[n-1], n >= 1
        T _u;       // u[1] = 1, u[n] = v[n]/v[n-1], n >= 2
        T _b;       // last denominator, i.e. b[n-1]
    };

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
    template <typename T, typename FwdIt>
    SPECFUN_HOST_DEVICE std::pair<T, std::size_t> continued_fraction_eval_series(
        T init_val, FwdIt cf, real_type_t<T> tol, std::size_t max_terms) {
        ContinuedFractionDifferenceGenerator<T, FwdIt> vs(cf);
        return series_eval_kahan(init_val, vs, tol, max_terms);
    }

} // namespace detail
} // namespace special
