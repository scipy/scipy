/* Building blocks for implementing special functions */

#pragma once

#include "config.h"
#include "error.h"

namespace special {
namespace detail {

    // Series evaluators.
    template <typename Generator, typename Number>
    SPECFUN_HOST_DEVICE Number series_eval(Generator &g, Number init_val, double tol, std::uint64_t max_terms,
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
        Number result = init_val;
        Number term, previous;
        for (std::uint64_t i = 0; i < max_terms; ++i) {
            term = g();
            result += term;
            if (std::abs(term) < std::abs(result) * tol) {
                return result;
            }
        }
        // Exceeded max terms without converging. Return NaN.
        set_error(func_name, SF_ERROR_NO_RESULT, NULL);
        if constexpr (std::is_floating_point<Number>::value) {
            return std::numeric_limits<Number>::quiet_NaN();
        }
        // If result type is not a floating point type, assume it is complex.
        using FloatingType = typename Number::value_type;
        return Number(std::numeric_limits<FloatingType>::quiet_NaN(), std::numeric_limits<FloatingType>::quiet_NaN());
    }

    template <typename Generator, typename Number>
    SPECFUN_HOST_DEVICE Number series_eval_fixed_length(Generator &g, Number init_val, std::uint64_t num_terms) {
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
        Number result = init_val;
        for (std::uint64_t i = 0; i < num_terms; ++i) {
            result += g();
        }
        return result;
    }

} // namespace detail
} // namespace special
