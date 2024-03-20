/* Building blocks for implementing special functions */

#pragma once

#include "config.h"
#include "error.h"

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

} // namespace detail
} // namespace special
