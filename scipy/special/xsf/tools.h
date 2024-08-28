/* Building blocks for implementing special functions */

#pragma once

#include "config.h"
#include "error.h"

namespace xsf {
namespace detail {

    /* Result type of a "generator", a callable object that produces a value
     * each time it is called.
     */
    template <typename Generator>
    using generator_result_t = typename std::decay<typename std::invoke_result<Generator>::type>::type;

    /* Used to deduce the type of the numerator/denominator of a fraction. */
    template <typename Pair>
    struct pair_traits;

    template <typename T>
    struct pair_traits<std::pair<T, T>> {
        using value_type = T;
    };

    template <typename Pair>
    using pair_value_t = typename pair_traits<Pair>::value_type;

    /* Used to extract the "value type" of a complex type. */
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
    XSF_HOST_DEVICE inline typename std::enable_if<std::is_floating_point<T>::value, T>::type maybe_complex_NaN() {
        return std::numeric_limits<T>::quiet_NaN();
    }

    template <typename T>
    XSF_HOST_DEVICE inline typename std::enable_if<!std::is_floating_point<T>::value, T>::type maybe_complex_NaN() {
        using V = typename T::value_type;
        return {std::numeric_limits<V>::quiet_NaN(), std::numeric_limits<V>::quiet_NaN()};
    }

    // Series evaluators.
    template <typename Generator, typename T = generator_result_t<Generator>>
    XSF_HOST_DEVICE T
    series_eval(Generator &g, T init_val, real_type_t<T> tol, std::uint64_t max_terms, const char *func_name) {
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
        T term;
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

    template <typename Generator, typename T = generator_result_t<Generator>>
    XSF_HOST_DEVICE T series_eval_fixed_length(Generator &g, T init_val, std::uint64_t num_terms) {
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
    XSF_HOST_DEVICE void kahan_step(T &sum, T &comp, T x) {
        T y = x - comp;
        T t = sum + y;
        comp = (t - sum) - y;
        sum = t;
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
     *   g
     *       Reference to generator that yields the sequence of values a[1],
     *       a[2], a[3], ...
     *
     *   tol
     *       Relative tolerance for convergence.  Specifically, stop iteration
     *       as soon as `abs(a[n]) <= tol * abs(S[n])` for some n >= 1.
     *
     *   max_terms
     *       Maximum number of terms after a[0] to evaluate.  It should be set
     *       large enough such that the convergence criterion is guaranteed
     *       to have been satisfied within that many terms if there is no
     *       rounding error.
     *
     *   init_val
     *       a[0].  Default is zero.  The type of this parameter (T) is used
     *       for intermediary computations as well as the result.
     *
     * Return Value
     * ------------
     * If the convergence criterion is satisfied by some `n <= max_terms`,
     * returns `(S[n], n)`.  Otherwise, returns `(S[max_terms], 0)`.
     */
    template <typename Generator, typename T = generator_result_t<Generator>>
    XSF_HOST_DEVICE std::pair<T, std::uint64_t>
    series_eval_kahan(Generator &&g, real_type_t<T> tol, std::uint64_t max_terms, T init_val = T(0)) {

        T sum = init_val;
        T comp = 0;
        for (std::uint64_t i = 0; i < max_terms; ++i) {
            T term = g();
            kahan_step(sum, comp, term);
            if (std::abs(term) <= tol * std::abs(sum)) {
                return {sum, i + 1};
            }
        }
        return {sum, 0};
    }

    /* Generator that yields the difference of successive convergents of a
     * continued fraction.
     *
     * Let f[n] denote the n-th convergent of a continued fraction:
     *
     *                 a[1]   a[2]       a[n]
     *   f[n] = b[0] + ------ ------ ... ----
     *                 b[1] + b[2] +     b[n]
     *
     * with f[0] = b[0].  This generator yields the sequence of values
     * f[1]-f[0], f[2]-f[1], f[3]-f[2], ...
     *
     * Constructor Arguments
     * ---------------------
     *   cf
     *       Reference to generator that yields the terms of the continued
     *       fraction as (numerator, denominator) pairs, starting from
     *       (a[1], b[1]).
     *
     *       `cf` must outlive the ContinuedFractionSeriesGenerator object.
     *
     *       The constructed object always eagerly retrieves the next term
     *       of the continued fraction.  Specifically, (a[1], b[1]) is
     *       retrieved upon construction, and (a[n], b[n]) is retrieved after
     *       (n-1) calls of `()`.
     *
     * Type Arguments
     * --------------
     *   T
     *       Type in which computations are performed and results are turned.
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
    template <typename Generator, typename T = pair_value_t<generator_result_t<Generator>>>
    class ContinuedFractionSeriesGenerator {

      public:
        XSF_HOST_DEVICE explicit ContinuedFractionSeriesGenerator(Generator &cf) : cf_(cf) { init(); }

        XSF_HOST_DEVICE double operator()() {
            double v = v_;
            advance();
            return v;
        }

      private:
        XSF_HOST_DEVICE void init() {
            auto [num, denom] = cf_();
            T a = num;
            T b = denom;
            u_ = T(1);
            v_ = a / b;
            b_ = b;
        }

        XSF_HOST_DEVICE void advance() {
            auto [num, denom] = cf_();
            T a = num;
            T b = denom;
            u_ = T(1) / (T(1) + (a * u_) / (b * b_));
            v_ *= (u_ - T(1));
            b_ = b;
        }

        Generator &cf_; // reference to continued fraction generator
        T v_;           // v[n] == f[n] - f[n-1], n >= 1
        T u_;           // u[1] = 1, u[n] = v[n]/v[n-1], n >= 2
        T b_;           // last denominator, i.e. b[n-1]
    };

    /* Converts a continued fraction into a series whose terms are the
     * difference of its successive convergents.
     *
     * See ContinuedFractionSeriesGenerator for details.
     */
    template <typename Generator, typename T = pair_value_t<generator_result_t<Generator>>>
    XSF_HOST_DEVICE ContinuedFractionSeriesGenerator<Generator, T> continued_fraction_series(Generator &cf) {
        return ContinuedFractionSeriesGenerator<Generator, T>(cf);
    }

    /* Find initial bracket for a bracketing scalar root finder. A valid bracket is a pair of points a < b for
     * which the signs of f(a) and f(b) differ. x_left and x_right give the initial bracket, and x_min and
     * x_max define the search space. This is a private function intended specifically for the situation where
     * the goal is to invert a CDF function F for a parametrized family of distributions with respect to one
     * parameter, when the other parameters are known, and where F is monotonic with respect to the unknown parameter.
     *
     * Note that this takes a pointer to a function taking a tuple of args along with a scalar
     * double argument. A tuple of args for specializing func is also passed as the final argument.
     * It would be much cleaner to use std::function and capturing lambda's to specialize the
     * function we are finding a root for, but this was found unworkable in CuPy using NVRTC.
     * This should be revisited in the future in order to allow simplifying this code. */
    template <typename... Args>
    XSF_HOST_DEVICE inline std::tuple<double, double, double, double, int> bracket_root_for_cdf_inversion(
        double (*func)(double, std::tuple<Args...>), double x_left, double x_right, double xmin, double xmax,
        double factor, bool increasing, std::uint64_t maxiter, std::tuple<Args...> args
    ) {
        double y_left = func(x_left, args), y_right = func(x_right, args);
        double y_left_sgn = std::signbit(y_left), y_right_sgn = std::signbit(y_right);

        if (y_left_sgn != y_right_sgn || (y_left == 0 || y_right == 0)) {
            /* The initial bracket is valid. */
            std::tuple<double, double, double, double, int> result(x_left, x_right, y_left, y_right, 0);
            return result;
        }
        bool search_left;
        /* The frontier is the new leading endpoint of the expanding bracket. The
         * interior endpoint trails behind the frontier. In each step, the old frontier
         * endpoint becomes the new interior endpoint. */
        double interior, frontier, y_interior, y_frontier, y_interior_sgn, y_frontier_sgn, boundary;
        if ((increasing && y_right < 0) || (!increasing && y_right > 0)) {
            /* If func is increasing  and func(x_right) < 0 or if func is decreasing and
             *  f(y_right) > 0, we should expand the bracket to the right. */
            interior = x_left, frontier = x_right, y_interior = y_left, y_frontier = y_right;
            y_interior_sgn = y_left_sgn;
            y_frontier_sgn = y_right_sgn;
            search_left = false;
            boundary = xmax;
        } else {
            /* Otherwise we move and expand the bracket to the left. */
            interior = x_right, frontier = x_left, y_interior = y_right, y_frontier = y_left;
            y_interior_sgn = y_right_sgn;
            y_frontier_sgn = y_left_sgn;
            search_left = true;
            boundary = xmin;
        }

        bool stop = false;
        for (std::uint64_t i = 0; i < maxiter; i++) {
            double step = (frontier - interior) * factor;
            interior = frontier;
            y_interior = y_frontier;
            y_interior_sgn = y_frontier_sgn;
            frontier += step;
            if ((search_left && frontier <= boundary) || (!search_left && frontier >= boundary)) {
                /* If the frontier has reached the boundary, set it to the boundary and signal
                 * the algorithm to stop. We cannot search further. */
                frontier = boundary;
                stop = true;
            }
            y_frontier = func(frontier, args);
            y_frontier_sgn = std::signbit(y_frontier);
            if (y_frontier_sgn != y_interior_sgn || (y_frontier == 0.0)) {
                /* Stopping condition, func evaluated at endpoints of bracket has opposing signs,
                 * meeting requirement for bracketing root finder. (Or endpoint has reached a
                 * zero.) */
                if (search_left) {
                    /* Ensure we return an interval (a, b) with a < b. */
                    std::swap(interior, frontier);
                    std::swap(y_interior, y_frontier);
                }
                std::tuple<double, double, double, double, int> result(interior, frontier, y_interior, y_frontier, 0);
                return result;
            }

            if (stop) {
                /* We've reached a boundary point without finding a root . */
                std::tuple<double, double, double, double, int> result(
                    std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                    search_left ? 1 : 2
                );
                return result;
            }
        }
        /* Failed to converge within maxiter iterations. If maxiter is sufficiently high and
         * factor1 and factor2 are set appropriately, this should only happen due to a bug
         * in this function. Limiting the number of iterations is a defensive programming
         * measure. */
        std::tuple<double, double, double, double, int> result(
            std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), 3
        );
        return result;
    }

    /* Find root of a scalar function using Chandrupatla's algorithm */
    template <typename... Args>
    XSF_HOST_DEVICE inline std::pair<double, int> find_root_chandrupatla(
        double (*func)(double, std::tuple<Args...>), double x1, double x2, double f1, double f2, double rtol,
        double atol, std::uint64_t maxiter, std::tuple<Args...> args
    ) {
        if (f1 == 0) {
            return {x1, 0};
        }
        if (f2 == 0) {
            return {x2, 0};
        }
        double t = 0.5, x3, f3;
        for (uint64_t i = 0; i < maxiter; i++) {
            double x = x1 + t * (x2 - x1);
            double f = func(x, args);
            if (std::signbit(f) == std::signbit(f1)) {
                x3 = x1;
                x1 = x;
                f3 = f1;
                f1 = f;
            } else {
                x3 = x2;
                x2 = x1;
                x1 = x;
                f3 = f2;
                f2 = f1;
                f1 = f;
            }
            double xm, fm;
            if (std::abs(f2) < std::abs(f1)) {
                xm = x2;
                fm = f2;
            } else {
                xm = x1;
                fm = f1;
            }
            double tol = 2.0 * rtol * std::abs(xm) + 0.5 * atol;
            double tl = tol / std::abs(x2 - x1);
            if (tl > 0.5 || fm == 0) {
                return {xm, 0};
            }
            double xi = (x1 - x2) / (x3 - x2);
            double phi = (f1 - f2) / (f3 - f2);
            double fl = 1.0 - std::sqrt(1.0 - xi);
            double fh = std::sqrt(xi);

            if ((fl < phi) && (phi < fh)) {
                t = (f1 / (f2 - f1)) * (f3 / (f2 - f3)) + (f1 / (f3 - f1)) * (f2 / (f3 - f2)) * ((x3 - x1) / (x2 - x1));
            } else {
                t = 0.5;
            }
            t = std::fmin(std::fmax(t, tl), 1.0 - tl);
        }
        return {std::numeric_limits<double>::quiet_NaN(), 1};
    }

} // namespace detail
} // namespace xsf
