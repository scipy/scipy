/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/* Compute the CDF of the Tukey-Lambda distribution
 * using a bracketing search with special checks
 *
 * The PPF of the Tukey-lambda distribution is
 * G(p) = (p**lam + (1-p)**lam) / lam
 *
 * Author:  Travis Oliphant
 */

#pragma once

#include "../config.h"

namespace xsf {
namespace cephes {

    namespace detail {

        constexpr double tukey_SMALLVAL = 1e-4;
        constexpr double tukey_EPS = 1.0e-14;
        constexpr int tukey_MAXCOUNT = 60;

    } // namespace detail

    XSF_HOST_DEVICE inline double tukeylambdacdf(double x, double lmbda) {
        double pmin, pmid, pmax, plow, phigh, xeval;
        int count;

        if (std::isnan(x) || std::isnan(lmbda)) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        xeval = 1.0 / lmbda;
        if (lmbda > 0.0) {
            if (x <= (-xeval)) {
                return 0.0;
            }
            if (x >= xeval) {
                return 1.0;
            }
        }

        if ((-detail::tukey_SMALLVAL < lmbda) && (lmbda < detail::tukey_SMALLVAL)) {
            if (x >= 0) {
                return 1.0 / (1.0 + std::exp(-x));
            } else {
                return exp(x) / (1.0 + exp(x));
            }
        }

        pmin = 0.0;
        pmid = 0.5;
        pmax = 1.0;
        plow = pmin;
        phigh = pmax;
        count = 0;

        while ((count < detail::tukey_MAXCOUNT) && (std::abs(pmid - plow) > detail::tukey_EPS)) {
            xeval = (std::pow(pmid, lmbda) - std::pow(1.0 - pmid, lmbda)) / lmbda;
            if (xeval == x) {
                return pmid;
            }
            if (xeval > x) {
                phigh = pmid;
                pmid = (pmid + plow) / 2.0;
            } else {
                plow = pmid;
                pmid = (pmid + phigh) / 2.0;
            }
            count++;
        }
        return pmid;
    }

} // namespace cephes
} // namespace xsf
