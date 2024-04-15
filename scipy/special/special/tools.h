/* Building blocks for implementing special functions */

#pragma once

#include "config.h"
#include "error.h"

namespace special {
namespace detail {

    /* Returns the next value of a sequence.
     *
     * This overload assumes `g` to be a (stateful) callable object that
     * returns the next value of the sequence each time it is called.
     * `g` behaves like a bound `__next__` method in Python, except that
     * no `StopIteration` is supported.
     */
    template <typename Generator>
    std::enable_if_t<std::is_invocable_v<Generator>, std::invoke_result_t<Generator>>
    next(Generator &g) { return g(); }

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
    SPECFUN_HOST_DEVICE T series_eval(Generator &g, T init_val, real_type_t<T> tol, std::uint64_t max_terms,
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
            term = next(g);
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
            result += next(g);
        }
        return result;
    }

    /* Performs one step of Kahan summation. */
    template <typename T>
    SPECFUN_HOST_DEVICE void kahan_step(T *sum, T *comp, T x) {
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
     * And for n = 0, 1, 2, ..., denote its n-th partial sum by
     *
     *   S[n] = a[0] + a[1] + ... + a[n]
     *
     * This function computes S[0], S[1], ... until a[n] is sufficiently
     * small or if the maximum number of terms have been evaluated.
     *
     * Parameters
     * ----------
     *   init_val
     *       a[0].  The type of this parameter (T) is used for intermediary
     *       computations as well as the return value.
     *
     *   g
     *       Reference to generator that yields the sequence of values a[1],
     *       a[2], a[3], ...
     *
     *   tol
     *       Relative tolerance for convergence.  Specifically, stop iteration
     *       as soon as `abs(a[n]) <= tol * abs(S[n])` for some n >= 1.
     *
     *   max_terms
     *       Maximum number of terms after a[0] to evaluate.  This parameter
     *       is used to handle unexpected non-convergence.  It should be set
     *       large enough such that the convergence criterion is guaranteed
     *       to have been satisfied within that many terms if there is no
     *       rounding error.
     *
     * Return Value
     * ------------
     * If the convergence criterion is satisfied by some `n <= max_terms`,
     * returns `(S[n], n)`.  Otherwise, returns `(S[max_terms], 0)`.
     */
    template <typename T, typename Generator>
    SPECFUN_HOST_DEVICE std::pair<T, std::size_t> series_eval_kahan(
        T init_val, Generator &g, real_type_t<T> tol, std::size_t max_terms) {

        T sum = init_val;
        T comp = 0;
        for (std::size_t i = 0; i < max_terms; ++i) {
            T term = next(g);
            kahan_step(&sum, &comp, term);
            if (std::abs(term) <= tol * std::abs(sum)) {
                return {sum, i + 1};
            }
        }
        return {sum, 0};
    }

    /* Evaluates a continued fraction using the modified Lentz algorithm.
     *
     * Let f denote a continued fraction:
     *
     *              a[1]   a[2]   a[3]
     *   f = b[0] + ------ ------ ------ ...
     *              b[1] + b[2] + b[3] +
     *
     * Parameters
     * ----------
     *   init_val
     *       b0.  The type of this parameter (T) is used for intermediary
     *       computations as well as the return value.
     *
     *   cf
     *       Reference to generator that yields the terms of the continued
     *       fraction as (numerator, denominator) pairs, starting from
     *       (a[1], b[1]).
     *
     *   tol
     *       Relative tolerance for convergence.  Specifically, stop iteration
     *       as soon as `abs(f[n]/f[n-1] - 1) <= abs(tol)` for some `n >= 1`.
     *
     *   max_terms
     *       Maximum number of terms to evaluate.  This parameter is used to
     *       handle unexpected non-convergence.  It should be set large enough
     *       such that the convergence criterion is guaranteed to have been
     *       satisfied within that many terms if there is no rounding error.
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
     * If the convergence criterion is satisfied by some `n <= max_terms`,
     * returns `(f[n], n)`.  Otherwise, returns `(f[max_terms], 0)`.
     *
     * Remarks
     * -------
     * No error checking is performed.  The caller must ensure that all terms
     * are finite and that intermediary computations do not trigger floating
     * point exceptions such as overflow.
     *
     * The numerical accuracy of this method depends on the characteristics of
     * the continued fraction being evaluated.
     */
    template <typename T, typename Generator>
    SPECFUN_HOST_DEVICE std::pair<T, std::size_t> continued_fraction_eval_lentz(
        T init_val, Generator &cf, T tol, std::size_t max_terms,
        T tiny_val = 16 * std::numeric_limits<real_type_t<T>>::min()) {

        T f = (init_val == 0)? tiny_val : init_val;
        T C = f;
        T D = 0;
        for (std::size_t n = 0; n < max_terms; ++n) {
            auto [num, denom] = next(cf);
            T a = num;
            T b = denom;

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
            if (std::abs(delta - 1) <= tol) {
                return {f, n+1};
            }
        }

        // Failed to converge within max_terms terms.
        return {f, 0};
    }

    /* Generator that yields the difference of successive convergents of a
     * continued fraction.
     *
     * Let f denote a continued fraction, and let f[n], n = 0, 1, 2, 3, ...
     * denote its n-th convergent.  Write
     *
     *              a[1]   a[2]   a[3]
     *   f = b[0] + ------ ------ ------ ...
     *              b[1] + b[2] + b[3] +
     *
     *                 a[1]   a[2]       a[n]
     *   f[n] = b[0] + ------ ------ ... ----
     *                 b[1] + b[2] +     b[n]
     *
     * with f[0] = b[0].  This generator yields the sequence of values
     * f[1]-f[0], f[2]-f[1], f[3]-f[2], ...
     *
     * Type Arguments
     * --------------
     *   T
     *       Type in which computations are performed and results are turned.
     *
     *   Generator
     *       Reference to generator that yields the terms of the continued
     *       fraction as (numerator, denominator) pairs, starting from
     *       (a[1], b[1]).
     *
     *       IMPORTANT: For performance reason, the generator always eagerly
     *       retrieves the next term of the continued fraction.  Specifically,
     *       (a[1], b[1]) is retrieved upon generator construction, and
     *       (a[n], b[n]) is retrieved after (n-1) calls of `()`.
     *
     * Remarks
     * -------
     * The series is computed using the recurrence relation described in [1].
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
     * [1] Gautschi, W. (1967). “Computational Aspects of Three-Term
     *     Recurrence Relations.” SIAM Review, 9(1):24-82.
     */
    template <typename T, typename Generator>
    class ContinuedFractionSeriesGenerator {

    public:
        explicit ContinuedFractionSeriesGenerator(Generator &cf) : cf_(cf) {
            init();
        }

        double operator()() {
            double v = v_;
            advance();
            return v;
        }

    private:
        void init() {
            auto [num, denom] = next(cf_);
            T a = num;
            T b = denom;
            u_ = 1;
            v_ = a / b;
            b_ = b;
        }

        void advance() {
            auto [num, denom] = next(cf_);
            T a = num;
            T b = denom;
            u_ = 1 / (1 + (a * u_) / (b * b_));
            v_ *= (u_ - 1);
            b_ = b;
        }

        Generator& cf_; // reference to continued fraction generator
        T v_;           // v[n] == f[n] - f[n-1], n >= 1
        T u_;           // u[1] = 1, u[n] = v[n]/v[n-1], n >= 2
        T b_;           // last denominator, i.e. b[n-1]
    };

    /* Converts a continued fraction into a series whose terms are the
     * difference of its successive convergents.
     *
     * See ContinuedFractionSeriesGenerator for details.
     */
    template <typename T, typename Generator>
    SPECFUN_HOST_DEVICE ContinuedFractionSeriesGenerator<T, Generator>
    continued_fraction_series(Generator &cf) {
        return ContinuedFractionSeriesGenerator<T, Generator>(cf);
    }

} // namespace detail
} // namespace special
