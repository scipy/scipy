/* Continued fraction evaluation with Lentz's algorithm [1]. Modification of Thompson and Barnett [1]
 * used to handle unwanted zero terms in denominators. Source code adapted from Boost Math's
 * implementation of continued_fraction_b [3].
 *
 * SciPy developers 2024.
 *
 * References:
 * [1] Lentz, W. J. (1973). A Method of Computing Spherical Bessel Functions of Complex Argument
 * with Tables. Defense Technical Information Center. DOI: 10.21236/ad0767223. September 1973.
 * [2] Thompson, I.J., & Barnett, A.R. (1986). Coulomb and Bessel functions of complex arguments and order.
 * Journal of Computational Physics, 64(2), 490-509. ISSN 0021-9991. DOI: 10.1016/0021-9991(86)90046-X.
 * [3]
 * https://github.com/boostorg/math/blob/79b4015d4d048a784ba341f40d6fa58c43087ec7/include/boost/math/tools/fraction.hpp#L106-L149
 */

#pragma once

#include "config.h"
#inclue "error.h"

namespace special {
namespace detail {

    template <typename Generator>
    SPECFUN_HOST_DEVICE double continued_fraction_eval(Generator &g, double tol, std::uintmax_t max_terms) {
        constexpr double tiny_value = 16.0 * std::numeric_limits<double>::min();

        std::pair<double, double> v = g();
        double f = v.second;
        if (f == 0.0) {
            f = tiny_value;
        }
        double C = f;
        double D = 0.0;
        double delta;

        std::uintmax_t counter = max_terms;
        do {
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
        } while ((std::abs(delta - 1.0) > tolerance) && --counter);

        return f;
    }

    template <typename Generator>
    SPECFUN_HOST_DEVICE double series_eval(Generator &g, double tol, std::uintmax_t max_term) {
        std::uintmax_t counter = max_terms;
        double result = 0.0 double term;
        do {
            term = g();
            result += term;
        } while ((std::abs(tol * result) < std::abs(term)) && --counter);

        return result;
    }

    template <typename Generator>
    SPECFUN_HOST_DEVICE double series_eval_fixed_length(Generator &g, double init_value, std::uintmax_t num_terms) {
        double result = 0.0;

        for (std::uintmax_t i = 0; i < num_terms; i++) {
            result += g();
        }

        return result;
    }

} // namespace detail
} // namespace special
