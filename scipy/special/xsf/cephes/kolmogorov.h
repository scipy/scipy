/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/* File altered for inclusion in cephes module for Python:
 * Main loop commented out.... */
/*  Travis Oliphant Nov. 1998 */

/* Re Kolmogorov statistics, here is Birnbaum and Tingey's (actually it was already present
 * in Smirnov's paper) formula for the
 * distribution of D+, the maximum of all positive deviations between a
 * theoretical distribution function P(x) and an empirical one Sn(x)
 * from n samples.
 *
 *     +
 *    D  =         sup     [P(x) - S (x)]
 *     n     -inf < x < inf         n
 *
 *
 *                  [n(1-d)]
 *        +            -                    v-1              n-v
 *    Pr{D   > d} =    >    C    d (d + v/n)    (1 - d - v/n)
 *        n            -   n v
 *                    v=0
 *
 * (also equals the following sum, but note the terms may be large and alternating in sign)
 * See Smirnov 1944, Dwass 1959
 *                         n
 *                         -                         v-1              n-v
 *                =  1 -   >         C    d (d + v/n)    (1 - d - v/n)
 *                         -        n v
 *                       v=[n(1-d)]+1
 *
 * [n(1-d)] is the largest integer not exceeding n(1-d).
 * nCv is the number of combinations of n things taken v at a time.

 * Sources:
 * [1] Smirnov, N.V. "Approximate laws of distribution of random variables from empirical data"
 *     Usp. Mat. Nauk, 1944. http://mi.mathnet.ru/umn8798
 * [2] Birnbaum, Z. W. and Tingey, Fred H.
 *     "One-Sided Confidence Contours for Probability Distribution Functions",
 *     Ann. Math. Statist. 1951.  https://doi.org/10.1214/aoms/1177729550
 * [3] Dwass, Meyer, "The Distribution of a Generalized $\mathrm{D}^+_n$ Statistic",
 *     Ann. Math. Statist., 1959.  https://doi.org/10.1214/aoms/1177706085
 * [4] van Mulbregt, Paul, "Computing the Cumulative Distribution Function and Quantiles of the One-sided
 Kolmogorov-Smirnov Statistic"
 *     http://arxiv.org/abs/1802.06966
 * [5] van Mulbregt, Paul,  "Computing the Cumulative Distribution Function and Quantiles of the limit of the Two-sided
 Kolmogorov-Smirnov Statistic"
 *     https://arxiv.org/abs/1803.00426
 *
 */

#pragma once

#include "../config.h"
#include "../error.h"

#include "const.h"
#include "dd_real.h"
#include "unity.h"

namespace xsf {
namespace cephes {

    namespace detail {
        /* ************************************************************************ */
        /* Algorithm Configuration */

        constexpr int KOLMOG_MAXITER = 500;

        /*
         * Kolmogorov Two-sided:
         * Switchover between the two series to compute K(x)
         *  0 <=  x <= KOLMOG_CUTOVER and
         *  KOLMOG_CUTOVER < x < infty
         */
        constexpr double KOLMOG_CUTOVER = 0.82;

        /*
         * Smirnov One-sided:
         * n larger than SMIRNOV_MAX_COMPUTE_N will result in an approximation
         */
        constexpr int SMIRNOV_MAX_COMPUTE_N = 1000000;

        /*
         * Use the upper sum formula, if the number of terms is at most SM_UPPER_MAX_TERMS,
         * and n is at least SM_UPPERSUM_MIN_N
         * Don't use the upper sum if lots of terms are involved as the series alternates
         *  sign and the terms get much bigger than 1.
         */
        constexpr int SM_UPPER_MAX_TERMS = 3;
        constexpr int SM_UPPERSUM_MIN_N = 10;

        /* ************************************************************************ */
        /* ************************************************************************ */

        /* exp() of anything below this returns 0 */
        constexpr int MIN_EXPABLE = (-708 - 38);

        /* Class to hold the CDF, SF and PDF, which are computed simultaneously */
        struct ThreeProbs {
            double sf;
            double cdf;
            double pdf;
        };

        constexpr double _xtol = std::numeric_limits<double>::epsilon();
        constexpr double _rtol = 2 * _xtol;

        XSF_HOST_DEVICE inline bool _within_tol(double x, double y, double atol, double rtol) {
            double diff = std::abs(x - y);
            bool result = (diff <= (atol + rtol * std::abs(y)));
            return result;
        }

        /* ************************************************************************ */
        /* Kolmogorov : Two-sided                      **************************** */
        /* ************************************************************************ */

        XSF_HOST_DEVICE inline ThreeProbs _kolmogorov(double x) {
            double P = 1.0;
            double D = 0;
            double sf, cdf, pdf;

            if (std::isnan(x)) {
                return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN()};
            }
            if (x <= 0) {
                return {1.0, 0.0, 0};
            }
            /* x <= 0.040611972203751713 */
            if (x <= M_PI / std::sqrt(-MIN_EXPABLE * 8)) {
                return {1.0, 0.0, 0};
            }

            P = 1.0;
            if (x <= KOLMOG_CUTOVER) {
                /*
                 *  u = e^(-pi^2/(8x^2))
                 *  w = sqrt(2pi)/x
                 *  P = w*u * (1 + u^8 + u^24 + u^48 + ...)
                 */
                double w = std::sqrt(2 * M_PI) / x;
                double logu8 = -M_PI * M_PI / (x * x); /* log(u^8) */
                double u = std::exp(logu8 / 8);
                if (u == 0) {
                    /*
                     * P = w*u, but u < 1e-308, and w > 1,
                     * so compute as logs, then exponentiate
                     */
                    double logP = logu8 / 8 + std::log(w);
                    P = std::exp(logP);
                } else {
                    /* Just unroll the loop, 3 iterations */
                    double u8 = std::exp(logu8);
                    double u8cub = std::pow(u8, 3);
                    P = 1 + u8cub * P;
                    D = 5 * 5 + u8cub * D;
                    P = 1 + u8 * u8 * P;
                    D = 3 * 3 + u8 * u8 * D;
                    P = 1 + u8 * P;
                    D = 1 * 1 + u8 * D;

                    D = M_PI * M_PI / 4 / (x * x) * D - P;
                    D *= w * u / x;
                    P = w * u * P;
                }
                cdf = P;
                sf = 1 - P;
                pdf = D;
            } else {
                /*
                 *  v = e^(-2x^2)
                 *  P = 2 (v - v^4 + v^9 - v^16 + ...)
                 *    = 2v(1 - v^3*(1 - v^5*(1 - v^7*(1 - ...)))
                 */
                double logv = -2 * x * x;
                double v = std::exp(logv);
                /*
                 * Want q^((2k-1)^2)(1-q^(4k-1)) / q(1-q^3) < epsilon to break out of loop.
                 * With KOLMOG_CUTOVER ~ 0.82, k <= 4.  Just unroll the loop, 4 iterations
                 */
                double vsq = v * v;
                double v3 = std::pow(v, 3);
                double vpwr;

                vpwr = v3 * v3 * v; /* v**7 */
                P = 1 - vpwr * P;   /* P <- 1 - (1-v**(2k-1)) * P */
                D = 3 * 3 - vpwr * D;

                vpwr = v3 * vsq;
                P = 1 - vpwr * P;
                D = 2 * 2 - vpwr * D;

                vpwr = v3;
                P = 1 - vpwr * P;
                D = 1 * 1 - vpwr * D;

                P = 2 * v * P;
                D = 8 * v * x * D;
                sf = P;
                cdf = 1 - sf;
                pdf = D;
            }
            pdf = std::fmax(0, pdf);
            cdf = std::clamp(cdf, 0.0, 1.0);
            sf = std::clamp(sf, 0.0, 1.0);
            return {sf, cdf, pdf};
        }

        /* Find x such kolmogorov(x)=psf, kolmogc(x)=pcdf */
        XSF_HOST_DEVICE inline double _kolmogi(double psf, double pcdf) {
            double x, t;
            double xmin = 0;
            double xmax = std::numeric_limits<double>::infinity();
            int iterations;
            double a = xmin, b = xmax;

            if (!(psf >= 0.0 && pcdf >= 0.0 && pcdf <= 1.0 && psf <= 1.0)) {
                set_error("kolmogi", SF_ERROR_DOMAIN, NULL);
                return (std::numeric_limits<double>::quiet_NaN());
            }
            if (std::abs(1.0 - pcdf - psf) > 4 * std::numeric_limits<double>::epsilon()) {
                set_error("kolmogi", SF_ERROR_DOMAIN, NULL);
                return (std::numeric_limits<double>::quiet_NaN());
            }
            if (pcdf == 0.0) {
                return 0.0;
            }
            if (psf == 0.0) {
                return std::numeric_limits<double>::infinity();
            }

            if (pcdf <= 0.5) {
                /* p ~ (sqrt(2pi)/x) *exp(-pi^2/8x^2).  Generate lower and upper bounds  */
                double logpcdf = std::log(pcdf);
                /* Now that 1 >= x >= sqrt(p) */
                /* Iterate twice: x <- pi/(sqrt(8) sqrt(log(sqrt(2pi)) - log(x) - log(pdf))) */
                a = M_PI / (2 * M_SQRT2 * std::sqrt(-(logpcdf + logpcdf / 2 - LOGSQRT2PI)));
                b = M_PI / (2 * M_SQRT2 * std::sqrt(-(logpcdf + 0 - LOGSQRT2PI)));
                a = M_PI / (2 * M_SQRT2 * std::sqrt(-(logpcdf + std::log(a) - LOGSQRT2PI)));
                b = M_PI / (2 * M_SQRT2 * std::sqrt(-(logpcdf + std::log(b) - LOGSQRT2PI)));
                x = (a + b) / 2.0;
            } else {
                /*
                 * Based on the approximation p ~ 2 exp(-2x^2)
                 * Found that needed to replace psf with a slightly smaller number in the second element
                 *  as otherwise _kolmogorov(b) came back as a very small number but with
                 *  the same sign as _kolmogorov(a)
                 *  kolmogi(0.5) = 0.82757355518990772
                 *  so (1-q^(-(4-1)*2*x^2)) = (1-exp(-6*0.8275^2) ~ (1-exp(-4.1)
                 */
                constexpr double jiggerb = 256 * std::numeric_limits<double>::epsilon();
                double pba = psf / (1.0 - std::exp(-4)) / 2, pbb = psf * (1 - jiggerb) / 2;
                double q0;
                a = std::sqrt(-0.5 * std::log(pba));
                b = std::sqrt(-0.5 * std::log(pbb));
                /*
                 * Use inversion of
                 *   p = q - q^4 + q^9 - q^16 + ...:
                 *   q = p + p^4 + 4p^7 - p^9 + 22p^10  - 13p^12 + 140*p^13 ...
                 */
                {
                    double p = psf / 2.0;
                    double p2 = p * p;
                    double p3 = p * p * p;
                    q0 = 1 + p3 * (1 + p3 * (4 + p2 * (-1 + p * (22 + p2 * (-13 + 140 * p)))));
                    q0 *= p;
                }
                x = std::sqrt(-std::log(q0) / 2);
                if (x < a || x > b) {
                    x = (a + b) / 2;
                }
            }
            XSF_ASSERT(a <= b);

            iterations = 0;
            do {
                double x0 = x;
                ThreeProbs probs = _kolmogorov(x0);
                double df = ((pcdf < 0.5) ? (pcdf - probs.cdf) : (probs.sf - psf));
                double dfdx;

                if (std::abs(df) == 0) {
                    break;
                }
                /* Update the bracketing interval */
                if (df > 0 && x > a) {
                    a = x;
                } else if (df < 0 && x < b) {
                    b = x;
                }

                dfdx = -probs.pdf;
                if (std::abs(dfdx) <= 0.0) {
                    x = (a + b) / 2;
                    t = x0 - x;
                } else {
                    t = df / dfdx;
                    x = x0 - t;
                }

                /*
                 * Check out-of-bounds.
                 * Not expecting this to happen often --- kolmogorov is convex near x=infinity and
                 * concave near x=0, and we should be approaching from the correct side.
                 * If out-of-bounds, replace x with a midpoint of the bracket.
                 */
                if (x >= a && x <= b) {
                    if (_within_tol(x, x0, _xtol, _rtol)) {
                        break;
                    }
                    if ((x == a) || (x == b)) {
                        x = (a + b) / 2.0;
                        /* If the bracket is already so small ... */
                        if (x == a || x == b) {
                            break;
                        }
                    }
                } else {
                    x = (a + b) / 2.0;
                    if (_within_tol(x, x0, _xtol, _rtol)) {
                        break;
                    }
                }

                if (++iterations > KOLMOG_MAXITER) {
                    set_error("kolmogi", SF_ERROR_SLOW, NULL);
                    break;
                }
            } while (1);
            return (x);
        }

    } // namespace detail

    XSF_HOST_DEVICE inline double kolmogorov(double x) {
        if (std::isnan(x)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return detail::_kolmogorov(x).sf;
    }

    XSF_HOST_DEVICE inline double kolmogc(double x) {
        if (std::isnan(x)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return detail::_kolmogorov(x).cdf;
    }

    XSF_HOST_DEVICE inline double kolmogp(double x) {
        if (std::isnan(x)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (x <= 0) {
            return -0.0;
        }
        return -detail::_kolmogorov(x).pdf;
    }

    /* Functional inverse of Kolmogorov survival statistic for two-sided test.
     * Finds x such that kolmogorov(x) = p.
     */
    XSF_HOST_DEVICE inline double kolmogi(double p) {
        if (std::isnan(p)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return detail::_kolmogi(p, 1 - p);
    }

    /* Functional inverse of Kolmogorov cumulative statistic for two-sided test.
     * Finds x such that kolmogc(x) = p = (or kolmogorov(x) = 1-p).
     */
    XSF_HOST_DEVICE inline double kolmogci(double p) {
        if (std::isnan(p)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return detail::_kolmogi(1 - p, p);
    }

    namespace detail {

        /* ************************************************************************ */
        /* ********** Smirnov : One-sided ***************************************** */
        /* ************************************************************************ */

        XSF_HOST_DEVICE inline double nextPowerOf2(double x) {
            double q = std::ldexp(x, 1 - std::numeric_limits<double>::digits);
            double L = std::abs(q + x);
            if (L == 0) {
                L = std::abs(x);
            } else {
                int Lint = (int) (L);
                if (Lint == L) {
                    L = Lint;
                }
            }
            return L;
        }

        XSF_HOST_DEVICE inline double modNX(int n, double x, int *pNXFloor, double *pNX) {
            /*
             * Compute floor(n*x) and remainder *exactly*.
             * If remainder is too close to 1 (E.g. (1, -std::numeric_limits<double>::epsilon()/2))
             *  round up and adjust   */
            double_double alphaD, nxD, nxfloorD;
            int nxfloor;
            double alpha;

            nxD = static_cast<double>(n) * double_double(x);
            nxfloorD = floor(nxD);
            alphaD = nxD - nxfloorD;
            alpha = alphaD.hi;
            nxfloor = static_cast<int>(nxfloorD);
            XSF_ASSERT(alpha >= 0);
            XSF_ASSERT(alpha <= 1);
            if (alpha == 1) {
                nxfloor += 1;
                alpha = 0;
            }
            XSF_ASSERT(alpha < 1.0);
            *pNX = static_cast<double>(nxD);
            *pNXFloor = nxfloor;
            return alpha;
        }

        /*
         * The binomial coefficient C overflows a 64 bit double, as the 11-bit
         * exponent is too small.
         * Store C as (Cman:double_double, Cexpt:int).
         *  I.e a Mantissa/significand, and an exponent.
         *  Cman lies between 0.5 and 1, and the exponent has >=32-bit.
         */
        XSF_HOST_DEVICE inline void updateBinomial(double_double *Cman, int *Cexpt, int n, int j) {
            int expt;
            double_double rat = double_double(n - j) / (j + 1.0);
            double_double man2 = *Cman * rat;
            man2 = frexp(man2, &expt);
            XSF_ASSERT(man2 != 0.0);
            *Cexpt += expt;
            *Cman = man2;
        }

        XSF_HOST_DEVICE double_double pow_D(const double_double &a, int m) {
            /*
             * Using dd_npwr() here would be quite time-consuming.
             * Tradeoff accuracy-time by using pow().
             */

            double ans, r, adj;
            if (m <= 0) {
                if (m == 0) {
                    return double_double(1.0);
                }
                return 1.0 / pow_D(a, -m);
            }
            if (a == 0.0) {
                return double_double(0.0);
            }
            ans = std::pow(a.hi, m);
            r = a.lo / a.hi;
            adj = m * r;
            if (std::abs(adj) > 1e-8) {
                if (std::abs(adj) < 1e-4) {
                    /* Take 1st two terms of Taylor Series for (1+r)^m */
                    adj += (m * r) * ((m - 1) / 2.0 * r);
                } else {
                    /* Take exp of scaled log */
                    adj = xsf::cephes::expm1(m * std::log1p(r));
                }
            }
            return double_double(ans) + ans * adj;
        }

        XSF_HOST_DEVICE inline double pow2(double a, double b, int m) {
            return static_cast<double>(pow_D(double_double(a) + b, m));
        }

        /*
         * Not 1024 as too big.  Want _MAX_EXPONENT < 1023-52 so as to keep both
         * elements of the double_double normalized
         */
        constexpr int SM_MAX_EXPONENT = 960;

        XSF_HOST_DEVICE double_double pow2Scaled_D(const double_double &a, int m, int *pExponent) {
            /* Compute a^m = significand*2^expt and return as (significand, expt) */
            double_double ans, y;
            int ansE, yE;
            int maxExpt = SM_MAX_EXPONENT;
            int q, r, y2mE, y2rE, y2mqE;
            double_double y2r, y2m, y2mq;

            if (m <= 0) {
                int aE1, aE2;
                if (m == 0) {
                    *pExponent = 0.0;
                    return double_double(1.0);
                }
                ans = pow2Scaled_D(a, -m, &aE1);
                ans = frexp(1.0 / ans, &aE2);
                ansE = -aE1 + aE2;
                *pExponent = ansE;
                return ans;
            }
            y = frexp(a, &yE);
            if (m == 1) {
                *pExponent = yE;
                return y;
            }
            /*
             *  y ^ maxExpt >= 2^{-960}
             *  =>  maxExpt = 960 / log2(y.x[0]) = 708 / log(y.x[0])
             *              = 665/((1-y.x[0] + y.x[0]^2/2 - ...)
             *              <= 665/(1-y.x[0])
             * Quick check to see if we might need to break up the exponentiation
             */
            if (m * (y.hi - 1) / y.hi < -SM_MAX_EXPONENT * M_LN2) {
                /* Now do it carefully, calling log() */
                double lg2y = std::log(y.hi) / M_LN2;
                double lgAns = m * lg2y;
                if (lgAns <= -SM_MAX_EXPONENT) {
                    maxExpt = static_cast<int>(nextPowerOf2(-SM_MAX_EXPONENT / lg2y + 1) / 2);
                }
            }
            if (m <= maxExpt) {
                double_double ans1 = pow_D(y, m);
                ans = frexp(ans1, &ansE);
                ansE += m * yE;
                *pExponent = ansE;
                return ans;
            }

            q = m / maxExpt;
            r = m % maxExpt;
            /* y^m = (y^maxExpt)^q * y^r */
            y2r = pow2Scaled_D(y, r, &y2rE);
            y2m = pow2Scaled_D(y, maxExpt, &y2mE);
            y2mq = pow2Scaled_D(y2m, q, &y2mqE);
            ans = frexp(y2r * y2mq, &ansE);
            y2mqE += y2mE * q;
            ansE += y2mqE + y2rE;
            ansE += m * yE;
            *pExponent = ansE;
            return ans;
        }

        XSF_HOST_DEVICE inline double_double pow4_D(double a, double b, double c, double d, int m) {
            /* Compute ((a+b)/(c+d)) ^ m */
            double_double A, C, X;
            if (m <= 0) {
                if (m == 0) {
                    return double_double(1.0);
                }
                return pow4_D(c, d, a, b, -m);
            }
            A = double_double(a) + b;
            C = double_double(c) + d;
            if (A == 0.0) {
                return (C == 0.0) ? quiet_NaN() : double_double(0.0);
            }
            if (C == 0.0) {
                return ((A < 0) ? -infinity() : infinity());
            }
            X = A / C;
            return pow_D(X, m);
        }

        XSF_HOST_DEVICE inline double pow4(double a, double b, double c, double d, int m) {
            double_double ret = pow4_D(a, b, c, d, m);
            return static_cast<double>(ret);
        }

        XSF_HOST_DEVICE inline double_double logpow4_D(double a, double b, double c, double d, int m) {
            /*
             * Compute log(((a+b)/(c+d)) ^ m)
             *    == m * log((a+b)/(c+d))
             *    == m * log( 1 + (a+b-c-d)/(c+d))
             */
            double_double ans;
            double_double A, C, X;
            if (m == 0) {
                return double_double(0.0);
            }
            A = double_double(a) + b;
            C = double_double(c) + d;
            if (A == 0.0) {
                return ((C == 0.0) ? double_double(0.0) : -infinity());
            }
            if (C == 0.0) {
                return infinity();
            }
            X = A / C;
            XSF_ASSERT(X.hi >= 0);
            if (0.5 <= X.hi && X.hi <= 1.5) {
                double_double A1 = A - C;
                double_double X1 = A1 / C;
                ans = log1p(X1);
            } else {
                ans = log(X);
            }
            ans = m * ans;
            return ans;
        }

        XSF_HOST_DEVICE inline double logpow4(double a, double b, double c, double d, int m) {
            double_double ans = logpow4_D(a, b, c, d, m);
            return static_cast<double>(ans);
        }

        /*
         *  Compute a single term in the summation, A_v(n, x):
         *  A_v(n, x) =  Binomial(n,v) * (1-x-v/n)^(n-v) * (x+v/n)^(v-1)
         */
        XSF_HOST_DEVICE inline void computeAv(int n, double x, int v, const double_double &Cman, int Cexpt,
                                                  double_double *pt1, double_double *pt2, double_double *pAv) {
            int t1E, t2E, ansE;
            double_double Av;
            double_double t2x = double_double(n - v) / n - x; /*  1 - x - v/n */
            double_double t2 = pow2Scaled_D(t2x, n - v, &t2E);
            double_double t1x = double_double(v) / n + x; /* x + v/n */
            double_double t1 = pow2Scaled_D(t1x, v - 1, &t1E);
            double_double ans = t1 * t2;
            ans = ans * Cman;
            ansE = Cexpt + t1E + t2E;
            Av = ldexp(ans, ansE);
            *pAv = Av;
            *pt1 = t1;
            *pt2 = t2;
        }

        XSF_HOST_DEVICE inline ThreeProbs _smirnov(int n, double x) {
            double nx, alpha;
            double_double AjSum = double_double(0.0);
            double_double dAjSum = double_double(0.0);
            double cdf, sf, pdf;

            int bUseUpperSum;
            int nxfl, n1mxfl, n1mxceil;

            if (!(n > 0 && x >= 0.0 && x <= 1.0)) {
                return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN()};
            }
            if (n == 1) {
                return {1 - x, x, 1.0};
            }
            if (x == 0.0) {
                return {1.0, 0.0, 1.0};
            }
            if (x == 1.0) {
                return {0.0, 1.0, 0.0};
            }

            alpha = modNX(n, x, &nxfl, &nx);
            n1mxfl = n - nxfl - (alpha == 0 ? 0 : 1);
            n1mxceil = n - nxfl;
            /*
             * If alpha is 0, don't actually want to include the last term
             * in either the lower or upper summations.
             */
            if (alpha == 0) {
                n1mxfl -= 1;
                n1mxceil += 1;
            }

            /* Special case:  x <= 1/n  */
            if (nxfl == 0 || (nxfl == 1 && alpha == 0)) {
                double t = pow2(1, x, n - 1);
                pdf = (nx + 1) * t / (1 + x);
                cdf = x * t;
                sf = 1 - cdf;
                /* Adjust if x=1/n *exactly* */
                if (nxfl == 1) {
                    XSF_ASSERT(alpha == 0);
                    pdf -= 0.5;
                }
                return {sf, cdf, pdf};
            }
            /* Special case:  x is so big, the sf underflows double64 */
            if (-2 * n * x * x < MINLOG) {
                return {0, 1, 0};
            }
            /* Special case:  x >= 1 - 1/n */
            if (nxfl >= n - 1) {
                sf = pow2(1, -x, n);
                cdf = 1 - sf;
                pdf = n * sf / (1 - x);
                return {sf, cdf, pdf};
            }
            /* Special case:  n is so big, take too long to compute */
            if (n > SMIRNOV_MAX_COMPUTE_N) {
                /* p ~ e^(-(6nx+1)^2 / 18n) */
                double logp = -std::pow(6.0 * n * x + 1, 2) / 18.0 / n;
                /* Maximise precision for small p-value. */
                if (logp < -M_LN2) {
                    sf = std::exp(logp);
                    cdf = 1 - sf;
                } else {
                    cdf = -xsf::cephes::expm1(logp);
                    sf = 1 - cdf;
                }
                pdf = (6.0 * n * x + 1) * 2 * sf / 3;
                return {sf, cdf, pdf};
            }
            {
                /*
                 * Use the upper sum if n is large enough, and x is small enough and
                 * the number of terms is going to be small enough.
                 * Otherwise it just drops accuracy, about 1.6bits * nUpperTerms
                 */
                int nUpperTerms = n - n1mxceil + 1;
                bUseUpperSum = (nUpperTerms <= 1 && x < 0.5);
                bUseUpperSum = (bUseUpperSum || ((n >= SM_UPPERSUM_MIN_N) && (nUpperTerms <= SM_UPPER_MAX_TERMS) &&
                                                 (x <= 0.5 / std::sqrt(n))));
            }
            {
                int start = 0, step = 1, nTerms = n1mxfl + 1;
                int j, firstJ = 0;
                int vmid = n / 2;
                double_double Cman = double_double(1.0);
                int Cexpt = 0;
                double_double Aj, dAj, t1, t2, dAjCoeff;
                double_double oneOverX = double_double(1.0) / x;

                if (bUseUpperSum) {
                    start = n;
                    step = -1;
                    nTerms = n - n1mxceil + 1;

                    t1 = pow4_D(1, x, 1, 0, n - 1);
                    t2 = double_double(1.0);
                    Aj = t1;

                    dAjCoeff = (n - 1) / (double_double(1.0) + x);
                    dAjCoeff = dAjCoeff + oneOverX;
                } else {
                    t1 = oneOverX;
                    t2 = pow4_D(1, -x, 1, 0, n);
                    Aj = t2 / x;

                    dAjCoeff = (-1 - double_double(n - 1) * x) / (double_double(1.0) - x);
                    dAjCoeff = dAjCoeff / x;
                    dAjCoeff = dAjCoeff + oneOverX;
                }

                dAj = Aj * dAjCoeff;
                AjSum = AjSum + Aj;
                dAjSum = dAjSum + dAj;

                updateBinomial(&Cman, &Cexpt, n, 0);
                firstJ++;

                for (j = firstJ; j < nTerms; j += 1) {
                    int v = start + j * step;

                    computeAv(n, x, v, Cman, Cexpt, &t1, &t2, &Aj);

                    if (isfinite(Aj) && (Aj != 0.0)) {
                        /* coeff = 1/x + (j-1)/(x+j/n) - (n-j)/(1-x-j/n) */
                        dAjCoeff = (n * (v - 1)) / (double_double(nxfl + v) + alpha) -
                                   ((n - v) * n) / (double_double(n - nxfl - v) - alpha);
                        dAjCoeff = dAjCoeff + oneOverX;
                        dAj = Aj * dAjCoeff;

                        XSF_ASSERT(isfinite(Aj));
                        AjSum = AjSum + Aj;
                        dAjSum = dAjSum + dAj;
                    }
                    /* Safe to terminate early? */
                    if (Aj != 0.0) {
                        if (((4 * (nTerms - j) * std::abs(static_cast<double>(Aj))) <
                             (std::numeric_limits<double>::epsilon() * static_cast<double>(AjSum))) &&
                            (j != nTerms - 1)) {
                            break;
                        }
                    } else if (j > vmid) {
                        XSF_ASSERT(Aj == 0.0);
                        break;
                    }
                    updateBinomial(&Cman, &Cexpt, n, j);
                }
                XSF_ASSERT(isfinite(AjSum));
                XSF_ASSERT(isfinite(dAjSum));
                {
                    double_double derivD = x * dAjSum;
                    double_double probD = x * AjSum;
                    double deriv = static_cast<double>(derivD);
                    double prob = static_cast<double>(probD);

                    XSF_ASSERT(nx != 1 || alpha > 0);
                    if (step < 0) {
                        cdf = prob;
                        sf = 1 - prob;
                        pdf = deriv;
                    } else {
                        cdf = 1 - prob;
                        sf = prob;
                        pdf = -deriv;
                    }
                }
            }
            pdf = std::fmax(0, pdf);
            cdf = std::clamp(cdf, 0.0, 1.0);
            sf = std::clamp(sf, 0.0, 1.0);
            return {sf, cdf, pdf};
        }

        /*
         * Functional inverse of Smirnov distribution
         * finds x such that smirnov(n, x) = psf; smirnovc(n, x) = pcdf).
         */
        XSF_HOST_DEVICE inline double _smirnovi(int n, double psf, double pcdf) {
            /*
             * Need to use a bracketing NR algorithm here and be very careful
             *  about the starting point.
             */
            double x, logpcdf;
            int iterations = 0;
            double a = 0, b = 1;
            double maxlogpcdf, psfrootn;
            double dx, dxold;

            if (!(n > 0 && psf >= 0.0 && pcdf >= 0.0 && pcdf <= 1.0 && psf <= 1.0)) {
                set_error("smirnovi", SF_ERROR_DOMAIN, NULL);
                return std::numeric_limits<double>::quiet_NaN();
            }
            if (std::abs(1.0 - pcdf - psf) > 4 * std::numeric_limits<double>::epsilon()) {
                set_error("smirnovi", SF_ERROR_DOMAIN, NULL);
                return (std::numeric_limits<double>::quiet_NaN());
            }
            /* STEP 1: Handle psf==0, or pcdf == 0 */
            if (pcdf == 0.0) {
                return 0.0;
            }
            if (psf == 0.0) {
                return 1.0;
            }
            /* STEP 2: Handle n=1 */
            if (n == 1) {
                return pcdf;
            }

            /* STEP 3 Handle psf *very* close to 0.  Correspond to (n-1)/n < x < 1  */
            psfrootn = std::pow(psf, 1.0 / n);
            /* xmin > 1 - 1.0 / n */
            if (n < 150 && n * psfrootn <= 1) {
                /* Solve exactly. */
                x = 1 - psfrootn;
                return x;
            }

            logpcdf = (pcdf < 0.5 ? std::log(pcdf) : std::log1p(-psf));

            /*
             * STEP 4 Find bracket and initial estimate for use in N-R
             * 4(a)  Handle 0 < x <= 1/n:   pcdf = x * (1+x)^*(n-1)
             */
            maxlogpcdf = logpow4(1, 0.0, n, 0, 1) + logpow4(n, 1, n, 0, n - 1);
            if (logpcdf <= maxlogpcdf) {
                double xmin = pcdf / SCIPY_El;
                double xmax = pcdf;
                double P1 = pow4(n, 1, n, 0, n - 1) / n;
                double R = pcdf / P1;
                double z0 = R;
                /*
                 * Do one iteration of N-R solving: z*e^(z-1) = R, with z0=pcdf/P1
                 *  z  <-  z - (z exp(z-1) - pcdf)/((z+1)exp(z-1))
                 *  If z_0 = R, z_1 = R(1-exp(1-R))/(R+1)
                 */
                if (R >= 1) {
                    /*
                     * R=1 is OK;
                     * R>1 can happen due to truncation error for x = (1-1/n)+-eps
                     */
                    R = 1;
                    x = R / n;
                    return x;
                }
                z0 = (z0 * z0 + R * std::exp(1 - z0)) / (1 + z0);
                x = z0 / n;
                a = xmin * (1 - 4 * std::numeric_limits<double>::epsilon());
                a = std::fmax(a, 0);
                b = xmax * (1 + 4 * std::numeric_limits<double>::epsilon());
                b = std::fmin(b, 1.0 / n);
                x = std::clamp(x, a, b);
            } else {
                /* 4(b) : 1/n < x < (n-1)/n */
                double xmin = 1 - psfrootn;
                double logpsf = (psf < 0.5 ? std::log(psf) : std::log1p(-pcdf));
                double xmax = std::sqrt(-logpsf / (2.0L * n));
                double xmax6 = xmax - 1.0L / (6 * n);
                a = xmin;
                b = xmax;
                /* Allow for a little rounding error */
                a *= 1 - 4 * std::numeric_limits<double>::epsilon();
                b *= 1 + 4 * std::numeric_limits<double>::epsilon();
                a = std::fmax(xmin, 1.0 / n);
                b = std::fmin(xmax, 1 - 1.0 / n);
                x = xmax6;
            }
            if (x < a || x > b) {
                x = (a + b) / 2;
            }
            XSF_ASSERT(x < 1);

            /*
             * Skip computing fa, fb as that takes cycles and the exact values
             * are not needed.
             */

            /* STEP 5 Run N-R.
             * smirnov should be well-enough behaved for NR starting at this location.
             * Use smirnov(n, x)-psf, or pcdf - smirnovc(n, x), whichever has smaller p.
             */
            dxold = b - a;
            dx = dxold;
            do {
                double dfdx, x0 = x, deltax, df;
                XSF_ASSERT(x < 1);
                XSF_ASSERT(x > 0);
                {
                    ThreeProbs probs = _smirnov(n, x0);
                    df = ((pcdf < 0.5) ? (pcdf - probs.cdf) : (probs.sf - psf));
                    dfdx = -probs.pdf;
                }
                if (df == 0) {
                    return x;
                }
                /* Update the bracketing interval */
                if (df > 0 && x > a) {
                    a = x;
                } else if (df < 0 && x < b) {
                    b = x;
                }

                if (dfdx == 0) {
                    /*
                     * x was not within tolerance, but now we hit a 0 derivative.
                     * This implies that x >> 1/sqrt(n), and even then |smirnovp| >= |smirnov|
                     * so this condition is unexpected.  Do a bisection step.
                     */
                    x = (a + b) / 2;
                    deltax = x0 - x;
                } else {
                    deltax = df / dfdx;
                    x = x0 - deltax;
                }
                /*
                 * Check out-of-bounds.
                 * Not expecting this to happen ofen --- smirnov is convex near x=1 and
                 * concave near x=0, and we should be approaching from the correct side.
                 * If out-of-bounds, replace x with a midpoint of the bracket.
                 * Also check fast enough convergence.
                 */
                if ((a <= x) && (x <= b) &&
                    (std::abs(2 * deltax) <= std::abs(dxold) ||
                     std::abs(dxold) < 256 * std::numeric_limits<double>::epsilon())) {
                    dxold = dx;
                    dx = deltax;
                } else {
                    dxold = dx;
                    dx = dx / 2;
                    x = (a + b) / 2;
                    deltax = x0 - x;
                }
                /*
                 * Note that if psf is close to 1, f(x) -> 1, f'(x) -> -1.
                 *  => abs difference |x-x0| is approx |f(x)-p| >= std::numeric_limits<double>::epsilon(),
                 *  => |x-x0|/x >= std::numeric_limits<double>::epsilon()/x.
                 *  => cannot use a purely relative criteria as it will fail for x close to 0.
                 */
                if (_within_tol(x, x0, (psf < 0.5 ? 0 : _xtol), _rtol)) {
                    break;
                }
                if (++iterations > KOLMOG_MAXITER) {
                    set_error("smirnovi", SF_ERROR_SLOW, NULL);
                    return (x);
                }
            } while (1);
            return x;
        }

    } // namespace detail

    XSF_HOST_DEVICE inline double smirnov(int n, double d) {
        if (std::isnan(d)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return detail::_smirnov(n, d).sf;
    }

    XSF_HOST_DEVICE inline double smirnovc(int n, double d) {
        if (std::isnan(d)) {
            return NAN;
        }
        return detail::_smirnov(n, d).cdf;
    }

    /*
     * Derivative of smirnov(n, d)
     *  One interior point of discontinuity at d=1/n.
     */
    XSF_HOST_DEVICE inline double smirnovp(int n, double d) {
        if (!(n > 0 && d >= 0.0 && d <= 1.0)) {
            return (std::numeric_limits<double>::quiet_NaN());
        }
        if (n == 1) {
            /* Slope is always -1 for n=1, even at d = 1.0 */
            return -1.0;
        }
        if (d == 1.0) {
            return -0.0;
        }
        /*
         * If d is 0, the derivative is discontinuous, but approaching
         * from the right the limit is -1
         */
        if (d == 0.0) {
            return -1.0;
        }
        return -detail::_smirnov(n, d).pdf;
    }

    XSF_HOST_DEVICE inline double smirnovi(int n, double p) {
        if (std::isnan(p)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return detail::_smirnovi(n, p, 1 - p);
    }

    XSF_HOST_DEVICE inline double smirnovci(int n, double p) {
        if (std::isnan(p)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return detail::_smirnovi(n, 1 - p, p);
    }

} // namespace cephes
} // namespace xsf
