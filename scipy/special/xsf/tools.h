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
     * which the signs of f(a) and f(b) differ. Assumes function is monotonic and it is known
     * whether the function is increasing or decreasing. Passing in whether the function is
     * increasing or not allows for simplication of logic and avoidance of some branches.
     *
     * Note that this takes a pointer to a function taking a tuple of args along with a scalar
     * double argument. A tuple of args for specializing func is also passed as the final argument.
     * It would be much cleaner to use std::function and capturing lambda's to specialize the
     * function we are finding a root for, but this was found unworkable in CuPy using NVRTC.
     * This should be revisited in the future in order to allow simplifying this code. */
    template <typename... Args>
    XSF_HOST_DEVICE inline std::tuple<double, double, int> bracket_root(
        double (*func)(double, std::tuple<Args...>), double x_left, double x_right, double x_min, double x_max,
        double factor, bool increasing, std::tuple<Args...> args
    ) {
        double y_left = func(x_left, args), y_right = func(x_right, args);
        double y_left_sgn = std::signbit(y_left), y_right_sgn = std::signbit(y_right);

        if (y_left_sgn != y_right_sgn || (y_left == 0 || y_right == 0)) {
            /* Check if the initial bracket is valid. */
            std::tuple<double, double, int> result(x_left, x_right, 0);
            return result;
        }
        bool search_left;
        /* The frontier is the new leading endpoint of the expanding bracket. The
         * interior endpoint trails behind the frontier. In each step, except in
         * a special circumstance, the old frontier endpoint becomes the new interior
         * endpoint. */
        double interior, frontier, y_interior, y_frontier, boundary;
        if ((increasing && y_right < 0) || (!increasing && y_right > 0)) {
            /* If func is increasing  and func(x_right) < 0 or if func is decreasing and
             *  f(y_right) > 0, we should expand the bracket to the right. */
            interior = x_left, frontier = x_right, y_interior = y_left, y_frontier = y_right;
            search_left = false;
            boundary = x_max;
        } else {
            /* Otherwise we move and expand the bracket to the left. */
            interior = x_right, frontier = x_left, y_interior = y_right, y_frontier = y_left;
            search_left = true;
            boundary = x_min;
        }
        bool y_interior_sgn = std::signbit(y_interior);
        bool y_frontier_sgn = std::signbit(y_frontier);

        /* If f(x_left) = f(x_right), we say we are in a plateau. In this case, only
         * move the frontier endpoint while keeping the interior endpoint the same (in
         * contrast to setting the interior endpoint to the old frontier). This helps
         * avoid moving into a new plateau. */
        bool plateau = (y_left == y_right);
        bool stop = false;
        while (true) {
            double step = (frontier - interior) * factor;
            if (!plateau) {
                interior = frontier;
                y_interior = y_frontier;
                y_interior_sgn = y_frontier_sgn;
            }
            frontier += step;
            if ((search_left && frontier <= boundary) || (!search_left && frontier >= boundary)) {
                /* If the frontier has reached the boundary, set it to the boundary and signal
                 * the algorithm to stop. We cannot search further. */
                frontier = boundary;
                stop = true;
            }
            y_frontier = func(frontier, args);
            y_frontier_sgn = std::signbit(y_frontier);
            if (y_frontier == y_interior) {
                plateau = true;
            }
            if (y_frontier_sgn != y_interior_sgn || (y_frontier == 0.0)) {
                /* Stopping condition, func evaluated at endpoints of bracket has opposing signs,
                 * meeting requirement for bracketing root finder. (Or endpoint has reached a
                 * zero.) */
                if (search_left) {
                    /* Ensure we return an interval (a, b) with a < b. */
                    std::swap(interior, frontier);
                }
                std::tuple<double, double, int> result(interior, frontier, 0);
                return result;
            }

            if (stop) {
                /* We have reached the boundary and must stop the search. If f evaluated at the boundary is
                 * a zero, then we consider this a valid bracket because we have found a root, otherwise we
                 have not found a valid bracket. */
                if (y_frontier == 0.0) {
                    if (search_left) {
                        std::swap(interior, frontier);
                    }
                    std::tuple<double, double, int> result(interior, frontier, 0);
                    return result;
                }
                std::tuple<double, double, int> result(
                    std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                    search_left ? 1 : 2
                );
                return result;
            }
        }
    }

    /* Find root of a real valued continuous function of a single variable
     *
     * This is algorithm R from the paper "Two Efficient Algorithms with Guaranteed Convergence
     * for Finding a Zero of a Function" by Bus and Dekker. This algorithm used Barry Brown,
     * James Lovato, and Kathy Russell's cdflib library available at https://www.netlib.org/random/,
     * but it was translated from the Algol code in Bus & Dekker, not the Fortran from cdflib.
     *
     * The algorithm is similar to Brent's method, and uses a mix of linear interpolation,
     * (secant method), rational interpolation, and bisection.
     *
     * Note that this takes a pointer to a function taking a tuple of args along with a scalar
     * double argument. A tuple of args for specializing func is also passed as the final argument.
     * It would be much cleaner to use std::function and capturing lambda's to specialize the
     * function we are finding a root for, but this was found unworkable in CuPy using NVRTC.
     * This should be revisited in the future in order to allow simplifying this code. */
    template <typename... Args>
    XSF_HOST_DEVICE inline std::pair<double, int>
    find_root_bus_dekker_r(double (*func)(double, std::tuple<Args...>), double a, double b, std::tuple<Args...> args) {
        double fa = func(a, args), fb = func(b, args);
        // Handle cases where zero is on endpoint of initial bracket.
        if (fa == 0) {
            return {a, 0};
        }
        if (fb == 0) {
            return {b, 0};
        }
        if (std::signbit(fa) == std::signbit(fb)) {
            // Initial bracket is invalid.
            if (std::abs(fa) == std::abs(fb)) {
                return {std::numeric_limits<double>::quiet_NaN(), 3};
            } else if (std::abs(fa) > std::abs(fb)) {
                // Answer is likely beyond upper bound
                return {b, 2};
            }
            // Answer is likely before lower bound.
            return {a, 1};
        }
        bool first = true;
        /* Bus and Dekker distinguish between what they call intrapolation steps
         *  when a == c and extrapolation steps when a != c. For first iteration attempt
         * an intrapolation step. */
        int ext = 0;
        double c = a, fc = fa;
        /* d stores the previous value of a. We initialize to avoid compiler warnings,
         * but the initial values here don't actually matter. */
        double d = 0, fd = 0;
        while (true) {
            if (std::abs(fc) < std::abs(fb)) {
                // interchange to ensure f(b) is the smallest value seen so far.
                if (c != a) {
                    d = a;
                    fd = fa;
                }
                a = b;
                fa = fb;
                b = c;
                fb = fc;
                c = a;
                fc = fa;
            }
            double m = 0.5 * (b + c);
            double mb = m - b;
            double eps = std::nextafter(b, std::numeric_limits<double>::infinity()) - b;
            if (std::abs(mb) <= eps) {
                /* If we've converged, check doubles one step in each direction */
                double best_val = std::abs(fb);
                double b_r = b + eps;
                double b_l = std::nextafter(b, -std::numeric_limits<double>::infinity());
                double candidate_val = std::abs(func(b_l, args));
                if (candidate_val <= best_val) {
                    b = b_l, best_val = candidate_val;
                }
                if (best_val > 0.0) {
                    candidate_val = std::abs(func(b_r, args));
                    if (candidate_val < best_val) {
                        b = b_r;
                    }
                }
                return {b, 0};
            }
            // w is the step used to update b
            double w;
            if (ext > 3) {
                /* w = mb corresponds to bisection. After 3 successive extrapolations, run
                 * a bisection step. This ensures desired asymptotic behavior. See
                 * Bus and Dekker Section 4.2.*/
                w = mb;
            } else {
                //
                double p = (b - a) * fb;
                double q;
                if (first) {
                    // Only on first iteration use a linear interpolation step.
                    q = fa - fb;
                    first = false;
                } else {
                    // Otherwise, perform three point rational interpolation.
                    double fdb = (fd - fb) / (d - b);
                    double fda = (fd - fa) / (d - a);
                    p = fda * p;
                    q = fdb * fa - fda * fb;
                }
                if (p < 0) {
                    p = -p;
                    q = -q;
                }
                if (ext == 3) {
                    /* On attempt at 3rd extrapolation, step is doubled. This is to ensure desired
                     * asymptotic behavior. See Bus and Dekker, Section 4.2. */
                    p *= 2.0;
                }
                if (p < std::abs(q) * eps) {
                    /* If step size p / q is too small for b + p / q to be
                     * different from b, choose minimal step size that will make them different. */
                    w = std::signbit(mb) * eps;
                } else if (p < mb * q) {
                    // Use the calculated step size if it will keep b within the current interval.
                    w = p / q;
                } else {
                    // Otherwise bisection.
                    w = mb;
                }
            }
            d = a;
            fd = fa;
            a = b;
            fa = fb;
            b += w;
            fb = func(b, args);
            if (std::signbit(fb) == std::signbit(fc)) {
                /* If f(b) and f(c) have the same sign, next step is intrapolation.
                 * Set c to a and reset extrapolation count. */
                c = a;
                fc = fa;
                ext = 0;
            } else {
                if (w == mb) {
                    // Reset extrapolation count if bisection was performed.
                    ext = 0;
                } else {
                    /* Otherwise, increment extrapolation count. Bisection will be used after 4
                     * successive extrapolations */
                    ext += 1;
                }
            }
        }
    }

} // namespace detail
} // namespace xsf
