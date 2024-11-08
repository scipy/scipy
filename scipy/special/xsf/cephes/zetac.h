/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     zetac.c
 *
 *     Riemann zeta function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, zetac();
 *
 * y = zetac( x );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 *
 *                inf.
 *                 -    -x
 *   zetac(x)  =   >   k   ,   x > 1,
 *                 -
 *                k=2
 *
 * is related to the Riemann zeta function by
 *
 *     Riemann zeta(x) = zetac(x) + 1.
 *
 * Extension of the function definition for x < 1 is implemented.
 * Zero is returned for x > log2(INFINITY).
 *
 * ACCURACY:
 *
 * Tabulated values have full machine accuracy.
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      1,50        10000       9.8e-16            1.3e-16
 *
 *
 */

/*
 * Cephes Math Library Release 2.1:  January, 1989
 * Copyright 1984, 1987, 1989 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */
#pragma once

#include "../config.h"

#include "const.h"
#include "lanczos.h"
#include "polevl.h"
#include "zeta.h"

namespace xsf {
namespace cephes {

    namespace detail {

        /* Riemann zeta(x) - 1
         * for integer arguments between 0 and 30.
         */
        constexpr double azetac[] = {-1.50000000000000000000E0, 0.0, /* Not used; zetac(1.0) is infinity. */
                                     6.44934066848226436472E-1, 2.02056903159594285400E-1,
                                     8.23232337111381915160E-2, 3.69277551433699263314E-2,
                                     1.73430619844491397145E-2, 8.34927738192282683980E-3,
                                     4.07735619794433937869E-3, 2.00839282608221441785E-3,
                                     9.94575127818085337146E-4, 4.94188604119464558702E-4,
                                     2.46086553308048298638E-4, 1.22713347578489146752E-4,
                                     6.12481350587048292585E-5, 3.05882363070204935517E-5,
                                     1.52822594086518717326E-5, 7.63719763789976227360E-6,
                                     3.81729326499983985646E-6, 1.90821271655393892566E-6,
                                     9.53962033872796113152E-7, 4.76932986787806463117E-7,
                                     2.38450502727732990004E-7, 1.19219925965311073068E-7,
                                     5.96081890512594796124E-8, 2.98035035146522801861E-8,
                                     1.49015548283650412347E-8, 7.45071178983542949198E-9,
                                     3.72533402478845705482E-9, 1.86265972351304900640E-9,
                                     9.31327432419668182872E-10};

        /* 2**x (1 - 1/x) (zeta(x) - 1) = P(1/x)/Q(1/x), 1 <= x <= 10 */
        constexpr double zetac_P[9] = {
            5.85746514569725319540E11, 2.57534127756102572888E11, 4.87781159567948256438E10,
            5.15399538023885770696E9,  3.41646073514754094281E8,  1.60837006880656492731E7,
            5.92785467342109522998E5,  1.51129169964938823117E4,  2.01822444485997955865E2,
        };

        constexpr double zetac_Q[8] = {
            /*  1.00000000000000000000E0, */
            3.90497676373371157516E11, 5.22858235368272161797E10, 5.64451517271280543351E9, 3.39006746015350418834E8,
            1.79410371500126453702E7,  5.66666825131384797029E5,  1.60382976810944131506E4, 1.96436237223387314144E2,
        };

        /* log(zeta(x) - 1 - 2**-x), 10 <= x <= 50 */
        constexpr double zetac_A[11] = {
            8.70728567484590192539E6,  1.76506865670346462757E8,   2.60889506707483264896E10,
            5.29806374009894791647E11, 2.26888156119238241487E13,  3.31884402932705083599E14,
            5.13778997975868230192E15, -1.98123688133907171455E15, -9.92763810039983572356E16,
            7.82905376180870586444E16, 9.26786275768927717187E16,
        };

        constexpr double zetac_B[10] = {
            /* 1.00000000000000000000E0, */
            -7.92625410563741062861E6,  -1.60529969932920229676E8,  -2.37669260975543221788E10,
            -4.80319584350455169857E11, -2.07820961754173320170E13, -2.96075404507272223680E14,
            -4.86299103694609136686E15, 5.34589509675789930199E15,  5.71464111092297631292E16,
            -1.79915597658676556828E16,
        };

        /* (1-x) (zeta(x) - 1), 0 <= x <= 1 */
        constexpr double zetac_R[6] = {
            -3.28717474506562731748E-1, 1.55162528742623950834E1, -2.48762831680821954401E2,
            1.01050368053237678329E3,   1.26726061410235149405E4, -1.11578094770515181334E5,
        };

        constexpr double zetac_S[5] = {
            /* 1.00000000000000000000E0, */
            1.95107674914060531512E1, 3.17710311750646984099E2, 3.03835500874445748734E3,
            2.03665876435770579345E4, 7.43853965136767874343E4,
        };

        constexpr double zetac_TAYLOR0[10] = {
            -1.0000000009110164892,    -1.0000000057646759799,
            -9.9999983138417361078e-1, -1.0000013011460139596,
            -1.000001940896320456,     -9.9987929950057116496e-1,
            -1.000785194477042408,     -1.0031782279542924256,
            -9.1893853320467274178e-1, -1.5,
        };

        constexpr int zetac_MAXL2 = 127;

        /*
         * Compute zetac for positive arguments
         */
        XSF_HOST_DEVICE inline double zetac_positive(double x) {
            int i;
            double a, b, s, w;

            if (x == 1.0) {
                return std::numeric_limits<double>::infinity();
            }

            if (x >= detail::zetac_MAXL2) {
                /* because first term is 2**-x */
                return 0.0;
            }

            /* Tabulated values for integer argument */
            w = std::floor(x);
            if (w == x) {
                i = x;
                if (i < 31) {
                    return (azetac[i]);
                }
            }

            if (x < 1.0) {
                w = 1.0 - x;
                a = xsf::cephes::polevl(x, zetac_R, 5) / (w * xsf::cephes::p1evl(x, zetac_S, 5));
                return a;
            }

            if (x <= 10.0) {
                b = std::pow(2.0, x) * (x - 1.0);
                w = 1.0 / x;
                s = (x * xsf::cephes::polevl(w, zetac_P, 8)) / (b * xsf::cephes::p1evl(w, zetac_Q, 8));
                return s;
            }

            if (x <= 50.0) {
                b = std::pow(2.0, -x);
                w = xsf::cephes::polevl(x, zetac_A, 10) / xsf::cephes::p1evl(x, zetac_B, 10);
                w = std::exp(w) + b;
                return w;
            }

            /* Basic sum of inverse powers */
            s = 0.0;
            a = 1.0;
            do {
                a += 2.0;
                b = std::pow(a, -x);
                s += b;
            } while (b / s > MACHEP);

            b = std::pow(2.0, -x);
            s = (s + b) / (1.0 - b);
            return s;
        }

        /*
         * Compute zetac for small negative x. We can't use the reflection
         * formula because to double precision 1 - x = 1 and zetac(1) = inf.
         */
        XSF_HOST_DEVICE inline double zetac_smallneg(double x) {
            return xsf::cephes::polevl(x, zetac_TAYLOR0, 9);
        }

        /*
         * Compute zetac using the reflection formula (see DLMF 25.4.2) plus
         * the Lanczos approximation for Gamma to avoid overflow.
         */
        XSF_HOST_DEVICE inline double zeta_reflection(double x) {
            double base, large_term, small_term, hx, x_shift;

            hx = x / 2;
            if (hx == std::floor(hx)) {
                /* Hit a zero of the sine factor */
                return 0;
            }

            /* Reduce the argument to sine */
            x_shift = std::fmod(x, 4);
            small_term = -SQRT2PI * sin(0.5 * M_PI * x_shift);
            small_term *= xsf::cephes::lanczos_sum_expg_scaled(x + 1) * xsf::cephes::zeta(x + 1, 1);

            /* Group large terms together to prevent overflow */
            base = (x + xsf::cephes::lanczos_g + 0.5) / (2 * M_PI * M_E);
            large_term = std::pow(base, x + 0.5);
            if (std::isfinite(large_term)) {
                return large_term * small_term;
            }
            /*
             * We overflowed, but we might be able to stave off overflow by
             * factoring in the small term earlier. To do this we compute
             *
             * (sqrt(large_term) * small_term) * sqrt(large_term)
             *
             * Since we only call this method for negative x bounded away from
             * zero, the small term can only be as small sine on that region;
             * i.e. about machine epsilon. This means that if the above still
             * overflows, then there was truly no avoiding it.
             */
            large_term = std::pow(base, 0.5 * x + 0.25);
            return (large_term * small_term) * large_term;
        }

    } // namespace detail

    /*
     * Riemann zeta function, minus one
     */
    XSF_HOST_DEVICE inline double zetac(double x) {
        if (std::isnan(x)) {
            return x;
        } else if (x == -std::numeric_limits<double>::infinity()) {
            return std::numeric_limits<double>::quiet_NaN();
        } else if (x < 0.0 && x > -0.01) {
            return detail::zetac_smallneg(x);
        } else if (x < 0.0) {
            return detail::zeta_reflection(-x) - 1;
        } else {
            return detail::zetac_positive(x);
        }
    }

    /*
     * Riemann zeta function
     */
    XSF_HOST_DEVICE inline double riemann_zeta(double x) {
        if (std::isnan(x)) {
            return x;
        } else if (x == -std::numeric_limits<double>::infinity()) {
            return std::numeric_limits<double>::quiet_NaN();
        } else if (x < 0.0 && x > -0.01) {
            return 1 + detail::zetac_smallneg(x);
        } else if (x < 0.0) {
            return detail::zeta_reflection(-x);
        } else {
            return 1 + detail::zetac_positive(x);
        }
    }

} // namespace cephes
} // namespace xsf
