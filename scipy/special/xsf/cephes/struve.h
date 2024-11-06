/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*
 * Compute the Struve function.
 *
 * Notes
 * -----
 *
 * We use three expansions for the Struve function discussed in [1]:
 *
 * - power series
 * - expansion in Bessel functions
 * - asymptotic large-z expansion
 *
 * Rounding errors are estimated based on the largest terms in the sums.
 *
 * ``struve_convergence.py`` plots the convergence regions of the different
 * expansions.
 *
 * (i)
 *
 * Looking at the error in the asymptotic expansion, one finds that
 * it's not worth trying if z ~> 0.7 * v + 12 for v > 0.
 *
 * (ii)
 *
 * The Bessel function expansion tends to fail for |z| >~ |v| and is not tried
 * there.
 *
 * For Struve H it covers the quadrant v > z where the power series may fail to
 * produce reasonable results.
 *
 * (iii)
 *
 * The three expansions together cover for Struve H the region z > 0, v real.
 *
 * They also cover Struve L, except that some loss of precision may occur around
 * the transition region z ~ 0.7 |v|, v < 0, |v| >> 1 where the function changes
 * rapidly.
 *
 * (iv)
 *
 * The power series is evaluated in double-double precision. This fixes accuracy
 * issues in Struve H for |v| << |z| before the asymptotic expansion kicks in.
 * Moreover, it improves the Struve L behavior for negative v.
 *
 *
 * References
 * ----------
 * [1] NIST Digital Library of Mathematical Functions
 *     https://dlmf.nist.gov/11
 */

/*
 * Copyright (C) 2013  Pauli Virtanen
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * a. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * b. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * c. Neither the name of Enthought nor the names of the SciPy Developers
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */
#pragma once

#include "../bessel.h"
#include "../config.h"
#include "../error.h"

#include "dd_real.h"
#include "gamma.h"
#include "scipy_iv.h"

namespace xsf {
namespace cephes {

    namespace detail {

        constexpr int STRUVE_MAXITER = 10000;
        constexpr double STRUVE_SUM_EPS = 1e-16; /* be sure we are in the tail of the sum */
        constexpr double STRUVE_SUM_TINY = 1e-100;
        constexpr double STRUVE_GOOD_EPS = 1e-12;
        constexpr double STRUVE_ACCEPTABLE_EPS = 1e-7;
        constexpr double STRUVE_ACCEPTABLE_ATOL = 1e-300;

        /*
         * Large-z expansion for Struve H and L
         * https://dlmf.nist.gov/11.6.1
         */
        XSF_HOST_DEVICE inline double struve_asymp_large_z(double v, double z, int is_h, double *err) {
            int n, sgn, maxiter;
            double term, sum, maxterm;
            double m;

            if (is_h) {
                sgn = -1;
            } else {
                sgn = 1;
            }

            /* Asymptotic expansion divergenge point */
            m = z / 2;
            if (m <= 0) {
                maxiter = 0;
            } else if (m > STRUVE_MAXITER) {
                maxiter = STRUVE_MAXITER;
            } else {
                maxiter = (int) m;
            }
            if (maxiter == 0) {
                *err = std::numeric_limits<double>::infinity();
                return std::numeric_limits<double>::quiet_NaN();
            }

            if (z < v) {
                /* Exclude regions where our error estimation fails */
                *err = std::numeric_limits<double>::infinity();
                return std::numeric_limits<double>::quiet_NaN();
            }

            /* Evaluate sum */
            term = -sgn / std::sqrt(M_PI) * std::exp(-xsf::cephes::lgam(v + 0.5) + (v - 1) * std::log(z / 2)) *
                   xsf::cephes::gammasgn(v + 0.5);
            sum = term;
            maxterm = 0;

            for (n = 0; n < maxiter; ++n) {
                term *= sgn * (1 + 2 * n) * (1 + 2 * n - 2 * v) / (z * z);
                sum += term;
                if (std::abs(term) > maxterm) {
                    maxterm = std::abs(term);
                }
                if (std::abs(term) < STRUVE_SUM_EPS * std::abs(sum) || term == 0 || !std::isfinite(sum)) {
                    break;
                }
            }

            if (is_h) {
                sum += xsf::cyl_bessel_y(v, z);
            } else {
                sum += xsf::cephes::iv(v, z);
            }

            /*
             * This error estimate is strictly speaking valid only for
             * n > v - 0.5, but numerical results indicate that it works
             * reasonably.
             */
            *err = std::abs(term) + std::abs(maxterm) * STRUVE_SUM_EPS;

            return sum;
        }

        /*
         * Power series for Struve H and L
         * https://dlmf.nist.gov/11.2.1
         *
         * Starts to converge roughly at |n| > |z|
         */
        XSF_HOST_DEVICE inline double struve_power_series(double v, double z, int is_h, double *err) {
            int n, sgn;
            double term, sum, maxterm, scaleexp, tmp;
            double_double cterm, csum, cdiv, z2, c2v, ctmp;

            if (is_h) {
                sgn = -1;
            } else {
                sgn = 1;
            }

            tmp = -xsf::cephes::lgam(v + 1.5) + (v + 1) * std::log(z / 2);
            if (tmp < -600 || tmp > 600) {
                /* Scale exponent to postpone underflow/overflow */
                scaleexp = tmp / 2;
                tmp -= scaleexp;
            } else {
                scaleexp = 0;
            }

            term = 2 / std::sqrt(M_PI) * std::exp(tmp) * xsf::cephes::gammasgn(v + 1.5);
            sum = term;
            maxterm = 0;

            cterm = double_double(term);
            csum = double_double(sum);
            z2 = double_double(sgn * z * z);
            c2v = double_double(2 * v);

            for (n = 0; n < STRUVE_MAXITER; ++n) {
                /* cdiv = (3 + 2*n) * (3 + 2*n + 2*v)) */
                cdiv = double_double(3 + 2 * n);
                ctmp = double_double(3 + 2 * n);
                ctmp = ctmp + c2v;
                cdiv = cdiv * ctmp;

                /* cterm *= z2 / cdiv */
                cterm = cterm * z2;
                cterm = cterm / cdiv;

                csum = csum + cterm;

                term = static_cast<double>(cterm);
                sum = static_cast<double>(csum);

                if (std::abs(term) > maxterm) {
                    maxterm = std::abs(term);
                }
                if (std::abs(term) < STRUVE_SUM_TINY * std::abs(sum) || term == 0 || !std::isfinite(sum)) {
                    break;
                }
            }

            *err = std::abs(term) + std::abs(maxterm) * 1e-22;

            if (scaleexp != 0) {
                sum *= std::exp(scaleexp);
                *err *= std::exp(scaleexp);
            }

            if (sum == 0 && term == 0 && v < 0 && !is_h) {
                /* Spurious underflow */
                *err = std::numeric_limits<double>::infinity();
                return std::numeric_limits<double>::quiet_NaN();
                ;
            }

            return sum;
        }

        /*
         * Bessel series
         * https://dlmf.nist.gov/11.4.19
         */
        XSF_HOST_DEVICE inline double struve_bessel_series(double v, double z, int is_h, double *err) {
            int n;
            double term, cterm, sum, maxterm;

            if (is_h && v < 0) {
                /* Works less reliably in this region */
                *err = std::numeric_limits<double>::infinity();
                return std::numeric_limits<double>::quiet_NaN();
            }

            sum = 0;
            maxterm = 0;

            cterm = std::sqrt(z / (2 * M_PI));

            for (n = 0; n < STRUVE_MAXITER; ++n) {
                if (is_h) {
                    term = cterm * xsf::cyl_bessel_j(n + v + 0.5, z) / (n + 0.5);
                    cterm *= z / 2 / (n + 1);
                } else {
                    term = cterm * xsf::cephes::iv(n + v + 0.5, z) / (n + 0.5);
                    cterm *= -z / 2 / (n + 1);
                }
                sum += term;
                if (std::abs(term) > maxterm) {
                    maxterm = std::abs(term);
                }
                if (std::abs(term) < STRUVE_SUM_EPS * std::abs(sum) || term == 0 || !std::isfinite(sum)) {
                    break;
                }
            }

            *err = std::abs(term) + std::abs(maxterm) * 1e-16;

            /* Account for potential underflow of the Bessel functions */
            *err += 1e-300 * std::abs(cterm);

            return sum;
        }

        XSF_HOST_DEVICE inline double struve_hl(double v, double z, int is_h) {
            double value[4], err[4], tmp;
            int n;

            if (z < 0) {
                n = v;
                if (v == n) {
                    tmp = (n % 2 == 0) ? -1 : 1;
                    return tmp * struve_hl(v, -z, is_h);
                } else {
                    return std::numeric_limits<double>::quiet_NaN();
                }
            } else if (z == 0) {
                if (v < -1) {
                    return xsf::cephes::gammasgn(v + 1.5) * std::numeric_limits<double>::infinity();
                } else if (v == -1) {
                    return 2 / std::sqrt(M_PI) / xsf::cephes::Gamma(0.5);
                } else {
                    return 0;
                }
            }

            n = -v - 0.5;
            if (n == -v - 0.5 && n > 0) {
                if (is_h) {
                    return (n % 2 == 0 ? 1 : -1) * xsf::cyl_bessel_j(n + 0.5, z);
                } else {
                    return xsf::cephes::iv(n + 0.5, z);
                }
            }

            /* Try the asymptotic expansion */
            if (z >= 0.7 * v + 12) {
                value[0] = struve_asymp_large_z(v, z, is_h, &err[0]);
                if (err[0] < STRUVE_GOOD_EPS * std::abs(value[0])) {
                    return value[0];
                }
            } else {
                err[0] = std::numeric_limits<double>::infinity();
            }

            /* Try power series */
            value[1] = struve_power_series(v, z, is_h, &err[1]);
            if (err[1] < STRUVE_GOOD_EPS * std::abs(value[1])) {
                return value[1];
            }

            /* Try bessel series */
            if (std::abs(z) < std::abs(v) + 20) {
                value[2] = struve_bessel_series(v, z, is_h, &err[2]);
                if (err[2] < STRUVE_GOOD_EPS * std::abs(value[2])) {
                    return value[2];
                }
            } else {
                err[2] = std::numeric_limits<double>::infinity();
            }

            /* Return the best of the three, if it is acceptable */
            n = 0;
            if (err[1] < err[n])
                n = 1;
            if (err[2] < err[n])
                n = 2;
            if (err[n] < STRUVE_ACCEPTABLE_EPS * std::abs(value[n]) || err[n] < STRUVE_ACCEPTABLE_ATOL) {
                return value[n];
            }

            /* Maybe it really is an overflow? */
            tmp = -xsf::cephes::lgam(v + 1.5) + (v + 1) * std::log(z / 2);
            if (!is_h) {
                tmp = std::abs(tmp);
            }
            if (tmp > 700) {
                set_error("struve", SF_ERROR_OVERFLOW, NULL);
                return std::numeric_limits<double>::infinity() * xsf::cephes::gammasgn(v + 1.5);
            }

            /* Failure */
            set_error("struve", SF_ERROR_NO_RESULT, NULL);
            return std::numeric_limits<double>::quiet_NaN();
        }
    } // namespace detail

    XSF_HOST_DEVICE inline double struve_h(double v, double z) { return detail::struve_hl(v, z, 1); }

    XSF_HOST_DEVICE inline double struve_l(double v, double z) { return detail::struve_hl(v, z, 0); }

} // namespace cephes
} // namespace xsf
