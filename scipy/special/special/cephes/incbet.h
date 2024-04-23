/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     incbet.c
 *
 *     Incomplete beta integral
 *
 *
 * SYNOPSIS:
 *
 * double a, b, x, y, incbet();
 *
 * y = incbet( a, b, x );
 *
 *
 * DESCRIPTION:
 *
 * Returns incomplete beta integral of the arguments, evaluated
 * from zero to x.  The function is defined as
 *
 *                  x
 *     -            -
 *    | (a+b)      | |  a-1     b-1
 *  -----------    |   t   (1-t)   dt.
 *   -     -     | |
 *  | (a) | (b)   -
 *                 0
 *
 * The domain of definition is 0 <= x <= 1.  In this
 * implementation a and b are restricted to positive values.
 * The integral from x to 1 may be obtained by the symmetry
 * relation
 *
 *    1 - incbet( a, b, x )  =  incbet( b, a, 1-x ).
 *
 * The integral is evaluated by a continued fraction expansion
 * or, when b*x is small, by a power series.
 *
 * ACCURACY:
 *
 * Tested at uniformly distributed random points (a,b,x) with a and b
 * in "domain" and x between 0 and 1.
 *                                        Relative error
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,5         10000       6.9e-15     4.5e-16
 *    IEEE      0,85       250000       2.2e-13     1.7e-14
 *    IEEE      0,1000      30000       5.3e-12     6.3e-13
 *    IEEE      0,10000    250000       9.3e-11     7.1e-12
 *    IEEE      0,100000    10000       8.7e-10     4.8e-11
 * Outputs smaller than the IEEE gradual underflow threshold
 * were excluded from these statistics.
 *
 * ERROR MESSAGES:
 *   message         condition      value returned
 * incbet domain      x<0, x>1          0.0
 * incbet underflow                     0.0
 */

/*
 * Cephes Math Library, Release 2.3:  March, 1995
 * Copyright 1984, 1995 by Stephen L. Moshier
 */
#pragma once

#include "../config.h"
#include "../error.h"

#include "beta.h"
#include "const.h"

namespace special {
namespace cephes {

    namespace detail {

        constexpr double incbet_big = 4.503599627370496e15;
        constexpr double incbet_biginv = 2.22044604925031308085e-16;

        /* Power series for incomplete beta integral.
         * Use when b*x is small and x not too close to 1.  */

        SPECFUN_HOST_DEVICE inline double incbet_pseries(double a, double b, double x) {
            double s, t, u, v, n, t1, z, ai;

            ai = 1.0 / a;
            u = (1.0 - b) * x;
            v = u / (a + 1.0);
            t1 = v;
            t = u;
            n = 2.0;
            s = 0.0;
            z = MACHEP * ai;
            while (std::abs(v) > z) {
                u = (n - b) * x / n;
                t *= u;
                v = t / (a + n);
                s += v;
                n += 1.0;
            }
            s += t1;
            s += ai;

            u = a * std::log(x);
            if ((a + b) < MAXGAM && std::abs(u) < MAXLOG) {
                t = 1.0 / beta(a, b);
                s = s * t * std::pow(x, a);
            } else {
                t = -lbeta(a, b) + u + std::log(s);
                if (t < MINLOG) {
                    s = 0.0;
                } else {
                    s = exp(t);
                }
            }
            return (s);
        }

        /* Continued fraction expansion #1 for incomplete beta integral */
        SPECFUN_HOST_DEVICE inline double incbcf(double a, double b, double x) {
            double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
            double k1, k2, k3, k4, k5, k6, k7, k8;
            double r, t, ans, thresh;
            int n;

            k1 = a;
            k2 = a + b;
            k3 = a;
            k4 = a + 1.0;
            k5 = 1.0;
            k6 = b - 1.0;
            k7 = k4;
            k8 = a + 2.0;

            pkm2 = 0.0;
            qkm2 = 1.0;
            pkm1 = 1.0;
            qkm1 = 1.0;
            ans = 1.0;
            r = 1.0;
            n = 0;
            thresh = 3.0 * MACHEP;
            do {

                xk = -(x * k1 * k2) / (k3 * k4);
                pk = pkm1 + pkm2 * xk;
                qk = qkm1 + qkm2 * xk;
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;

                xk = (x * k5 * k6) / (k7 * k8);
                pk = pkm1 + pkm2 * xk;
                qk = qkm1 + qkm2 * xk;
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;

                if (qk != 0) {
                    r = pk / qk;
                }
                if (r != 0) {
                    t = std::abs((ans - r) / r);
                    ans = r;
                } else {
                    t = 1.0;
                }
                if (t < thresh) {
                    goto cdone;
                }

                k1 += 1.0;
                k2 += 1.0;
                k3 += 2.0;
                k4 += 2.0;
                k5 += 1.0;
                k6 -= 1.0;
                k7 += 2.0;
                k8 += 2.0;

                if ((std::abs(qk) + std::abs(pk)) > incbet_big) {
                    pkm2 *= incbet_biginv;
                    pkm1 *= incbet_biginv;
                    qkm2 *= incbet_biginv;
                    qkm1 *= incbet_biginv;
                }
                if ((std::abs(qk) < incbet_biginv) || (fabs(pk) < incbet_biginv)) {
                    pkm2 *= incbet_big;
                    pkm1 *= incbet_big;
                    qkm2 *= incbet_big;
                    qkm1 *= incbet_big;
                }
            } while (++n < 300);

        cdone:
            return (ans);
        }

        /* Continued fraction expansion #2 for incomplete beta integral */
        SPECFUN_HOST_DEVICE inline double incbd(double a, double b, double x) {
            double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
            double k1, k2, k3, k4, k5, k6, k7, k8;
            double r, t, ans, z, thresh;
            int n;

            k1 = a;
            k2 = b - 1.0;
            k3 = a;
            k4 = a + 1.0;
            k5 = 1.0;
            k6 = a + b;
            k7 = a + 1.0;
            ;
            k8 = a + 2.0;

            pkm2 = 0.0;
            qkm2 = 1.0;
            pkm1 = 1.0;
            qkm1 = 1.0;
            z = x / (1.0 - x);
            ans = 1.0;
            r = 1.0;
            n = 0;
            thresh = 3.0 * MACHEP;
            do {

                xk = -(z * k1 * k2) / (k3 * k4);
                pk = pkm1 + pkm2 * xk;
                qk = qkm1 + qkm2 * xk;
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;

                xk = (z * k5 * k6) / (k7 * k8);
                pk = pkm1 + pkm2 * xk;
                qk = qkm1 + qkm2 * xk;
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;

                if (qk != 0)
                    r = pk / qk;
                if (r != 0) {
                    t = std::abs((ans - r) / r);
                    ans = r;
                } else {
                    t = 1.0;
                }
                if (t < thresh) {
                    goto cdone;
                }

                k1 += 1.0;
                k2 -= 1.0;
                k3 += 2.0;
                k4 += 2.0;
                k5 += 1.0;
                k6 += 1.0;
                k7 += 2.0;
                k8 += 2.0;

                if ((std::abs(qk) + std::abs(pk)) > incbet_big) {
                    pkm2 *= incbet_biginv;
                    pkm1 *= incbet_biginv;
                    qkm2 *= incbet_biginv;
                    qkm1 *= incbet_biginv;
                }
                if ((std::abs(qk) < incbet_biginv) || (std::abs(pk) < incbet_biginv)) {
                    pkm2 *= incbet_big;
                    pkm1 *= incbet_big;
                    qkm2 *= incbet_big;
                    qkm1 *= incbet_big;
                }
            } while (++n < 300);
        cdone:
            return (ans);
        }

    } // namespace detail

    SPECFUN_HOST_DEVICE inline double incbet(double aa, double bb, double xx) {
        double a, b, t, x, xc, w, y;
        int flag;

        if (aa <= 0.0 || bb <= 0.0)
            goto domerr;

        if ((xx <= 0.0) || (xx >= 1.0)) {
            if (xx == 0.0)
                return (0.0);
            if (xx == 1.0)
                return (1.0);
        domerr:
            set_error("incbet", SF_ERROR_DOMAIN, NULL);
            return (std::numeric_limits<double>::quiet_NaN());
        }

        flag = 0;
        if ((bb * xx) <= 1.0 && xx <= 0.95) {
            t = detail::incbet_pseries(aa, bb, xx);
            goto done;
        }

        w = 1.0 - xx;

        /* Reverse a and b if x is greater than the mean. */
        if (xx > (aa / (aa + bb))) {
            flag = 1;
            a = bb;
            b = aa;
            xc = xx;
            x = w;
        } else {
            a = aa;
            b = bb;
            xc = w;
            x = xx;
        }

        if (flag == 1 && (b * x) <= 1.0 && x <= 0.95) {
            t = detail::incbet_pseries(a, b, x);
            goto done;
        }

        /* Choose expansion for better convergence. */
        y = x * (a + b - 2.0) - (a - 1.0);
        if (y < 0.0) {
            w = detail::incbcf(a, b, x);
        } else {
            w = detail::incbd(a, b, x) / xc;
        }

        /* Multiply w by the factor
         * a      b   _             _     _
         * x  (1-x)   | (a+b) / ( a | (a) | (b) ) .   */

        y = a * std::log(x);
        t = b * std::log(xc);
        if ((a + b) < detail::MAXGAM && std::abs(y) < detail::MAXLOG && std::abs(t) < detail::MAXLOG) {
            t = std::pow(xc, b);
            t *= std::pow(x, a);
            t /= a;
            t *= w;
            t *= 1.0 / beta(a, b);
            goto done;
        }
        /* Resort to logarithms.  */
        y += t - lbeta(a, b);
        y += std::log(w / a);
        if (y < detail::MINLOG) {
            t = 0.0;
        } else {
            t = exp(y);
        }
    done:
        if (flag == 1) {
            if (t <= detail::MACHEP) {
                t = 1.0 - detail::MACHEP;
            } else {
                t = 1.0 - t;
            }
        }
        return (t);
    }

} // namespace cephes
} // namespace special
