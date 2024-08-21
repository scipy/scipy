/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 *
 * The parts of the qd double-double floating point package used in SciPy
 * have been reworked in a more modern C++ style using operator overloading.
 */

/*
 * include/double2.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract numbers DE-AC03-76SF00098 and
 * DE-AC02-05CH11231.
 *
 * Copyright (c) 2003-2009, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from U.S. Dept. of Energy) All rights reserved.
 *
 * By downloading or using this software you are agreeing to the modified
 * BSD license "BSD-LBNL-License.doc" (see LICENSE.txt).
 */
/*
 * Double-double precision (>= 106-bit significand) floating point
 * arithmetic package based on David Bailey's Fortran-90 double-double
 * package, with some changes. See
 *
 *   http://www.nersc.gov/~dhbailey/mpdist/mpdist.html
 *
 * for the original Fortran-90 version.
 *
 * Overall structure is similar to that of Keith Brigg's C++ double-double
 * package.  See
 *
 *   http://www-epidem.plansci.cam.ac.uk/~kbriggs/doubledouble.html
 *
 * for more details.  In particular, the fix for x86 computers is borrowed
 * from his code.
 *
 * Yozo Hida
 */

/*
 * This code taken from v2.3.18 of the qd package.
 */

#pragma once

#include "../config.h"

#include "unity.h"

namespace xsf {
namespace cephes {

    namespace detail {

        constexpr double __DD_SPLITTER = 134217729.0;               // = 2^27 + 1
        constexpr double __DD_SPLIT_THRESH = 6.69692879491417e+299; // = 2^996

        /*************************************************************************
         * The basic routines taking double arguments, returning 1 (or 2) doubles
         *************************************************************************/

        /* volatile is used below to prevent aggressive optimizations which may change
         * the result of the error calculations. These volatiles wer e included in the
         * original C code and may perhaps still be useful, e.g. if someone compiles with
         * --ffastmath.
         */

        /* Computes fl(a+b) and err(a+b).  Assumes |a| >= |b|. */
        XSF_HOST_DEVICE inline double quick_two_sum(double a, double b, double *err) {
            volatile double s = a + b;
            volatile double c = s - a;
            *err = b - c;
            return s;
        }

        /* Computes fl(a+b) and err(a+b).  */
        XSF_HOST_DEVICE inline double two_sum(double a, double b, double *err) {
            volatile double s = a + b;
            volatile double c = s - a;
            volatile double d = b - c;
            volatile double e = s - c;
            *err = (a - e) + d;
            return s;
        }

        /* Computes fl(a*b) and err(a*b). */
        XSF_HOST_DEVICE inline double two_prod(double a, double b, double *err) {
            volatile double p = a * b;
            *err = std::fma(a, b, -p);
            return p;
        }

        /* Computes fl(a*a) and err(a*a).  Faster than the above method. */
        XSF_HOST_DEVICE inline double two_sqr(double a, double *err) {
            volatile double p = a * a;
            *err = std::fma(a, a, -p);
            return p;
        }

        /* Computes the nearest integer to d. */
        XSF_HOST_DEVICE inline double two_nint(double d) {
            if (d == std::floor(d)) {
                return d;
            }
            return std::floor(d + 0.5);
        }

        struct double_double {
            double hi, lo;

            double_double() = default;
            double_double(double high, double low) : hi(high), lo(low) {}
            explicit double_double(double high) : hi(high), lo(0.0) {}

            XSF_HOST_DEVICE explicit operator double() const { return hi; }
            XSF_HOST_DEVICE explicit operator int() const { return static_cast<int>(hi); }
        };

        // Arithmetic operations

        XSF_HOST_DEVICE inline double_double operator-(const double_double &x) {
            return double_double(-x.hi, -x.lo);
        }

        XSF_HOST_DEVICE inline double_double operator+(const double_double &lhs, const double_double &rhs) {
            /* This one satisfies IEEE style error bound,
               due to K. Briggs and W. Kahan.                   */
            double s1, s2, t1, t2;

            s1 = two_sum(lhs.hi, rhs.hi, &s2);
            t1 = two_sum(lhs.lo, rhs.lo, &t2);
            s2 += t1;
            s1 = quick_two_sum(s1, s2, &s2);
            s2 += t2;
            s1 = quick_two_sum(s1, s2, &s2);
            return double_double(s1, s2);
        }

        XSF_HOST_DEVICE inline double_double operator+(const double_double &lhs, const double rhs) {
            double s1, s2;
            s1 = two_sum(lhs.hi, rhs, &s2);
            s2 += lhs.lo;
            s1 = quick_two_sum(s1, s2, &s2);
            return double_double(s1, s2);
        }

        XSF_HOST_DEVICE inline double_double operator+(const double lhs, const double_double &rhs) {
            double s1, s2;
            s1 = two_sum(lhs, rhs.hi, &s2);
            s2 += rhs.lo;
            s1 = quick_two_sum(s1, s2, &s2);
            return double_double(s1, s2);
        }

        XSF_HOST_DEVICE inline double_double operator-(const double_double &lhs, const double_double &rhs) {
            return lhs + (-rhs);
        }

        XSF_HOST_DEVICE inline double_double operator-(const double_double &lhs, const double rhs) {
            double s1, s2;
            s1 = two_sum(lhs.hi, -rhs, &s2);
            s2 += lhs.lo;
            s1 = quick_two_sum(s1, s2, &s2);
            return double_double(s1, s2);
        }

        XSF_HOST_DEVICE inline double_double operator-(const double lhs, const double_double &rhs) {
            double s1, s2;
            s1 = two_sum(lhs, -rhs.hi, &s2);
            s2 -= rhs.lo;
            s1 = quick_two_sum(s1, s2, &s2);
            return double_double(s1, s2);
        }

        XSF_HOST_DEVICE inline double_double operator*(const double_double &lhs, const double_double &rhs) {
            double p1, p2;
            p1 = two_prod(lhs.hi, rhs.hi, &p2);
            p2 += (lhs.hi * rhs.lo + lhs.lo * rhs.hi);
            p1 = quick_two_sum(p1, p2, &p2);
            return double_double(p1, p2);
        }

        XSF_HOST_DEVICE inline double_double operator*(const double_double &lhs, const double rhs) {
            double p1, p2, e1, e2;
            p1 = two_prod(lhs.hi, rhs, &e1);
            p2 = two_prod(lhs.lo, rhs, &e2);
            p1 = quick_two_sum(p1, e2 + p2 + e1, &e1);
            return double_double(p1, e1);
        }

        XSF_HOST_DEVICE inline double_double operator*(const double lhs, const double_double &rhs) {
            double p1, p2, e1, e2;
            p1 = two_prod(lhs, rhs.hi, &e1);
            p2 = two_prod(lhs, rhs.lo, &e2);
            p1 = quick_two_sum(p1, e2 + p2 + e1, &e1);
            return double_double(p1, e1);
        }

        XSF_HOST_DEVICE inline double_double operator/(const double_double &lhs, const double_double &rhs) {
            double q1, q2, q3;
            double_double r;

            q1 = lhs.hi / rhs.hi; /* approximate quotient */

            r = lhs - rhs * q1;

            q2 = r.hi / rhs.hi;
            r = r - rhs * q2;

            q3 = r.hi / rhs.hi;

            q1 = quick_two_sum(q1, q2, &q2);
            r = double_double(q1, q2) + q3;
            return r;
        }

        XSF_HOST_DEVICE inline double_double operator/(const double_double &lhs, const double rhs) {
            return lhs / double_double(rhs);
        }

        XSF_HOST_DEVICE inline double_double operator/(const double lhs, const double_double &rhs) {
            return double_double(lhs) / rhs;
        }

        XSF_HOST_DEVICE inline bool operator==(const double_double &lhs, const double_double &rhs) {
            return (lhs.hi == rhs.hi && lhs.lo == rhs.lo);
        }

        XSF_HOST_DEVICE inline bool operator==(const double_double &lhs, const double rhs) {
            return (lhs.hi == rhs && lhs.lo == 0.0);
        }

        XSF_HOST_DEVICE inline bool operator==(const double lhs, const double_double &rhs) {
            return (lhs == rhs.hi) && (rhs.lo == 0.0);
        }

        XSF_HOST_DEVICE inline bool operator!=(const double_double &lhs, const double_double &rhs) {
            return (lhs.hi != rhs.hi) || (lhs.lo != rhs.lo);
        }

        XSF_HOST_DEVICE inline bool operator!=(const double_double &lhs, const double rhs) {
            return (lhs.hi != rhs) || (lhs.lo != 0.0);
        }

        XSF_HOST_DEVICE inline bool operator!=(const double lhs, const double_double &rhs) {
            return (rhs.hi != lhs) || (rhs.lo != 0.0);
        }

        XSF_HOST_DEVICE inline bool operator<(const double_double &lhs, const double_double &rhs) {
            if (lhs.hi < rhs.hi) {
                return true;
            }
            if (lhs.hi > rhs.hi) {
                return false;
            }
            return lhs.lo < rhs.lo;
        }

        XSF_HOST_DEVICE inline bool operator<(const double_double &lhs, const double rhs) {
            if (lhs.hi < rhs) {
                return true;
            }
            if (lhs.hi > rhs) {
                return false;
            }
            return lhs.lo < 0.0;
        }

        template <typename T>
        XSF_HOST_DEVICE bool operator>(const double_double &lhs, const T &rhs) {
            return rhs < lhs;
        }

        XSF_HOST_DEVICE inline bool operator<(const double lhs, const double_double &rhs) { return rhs > lhs; }

        XSF_HOST_DEVICE inline bool operator>(const double lhs, const double_double &rhs) { return rhs < lhs; }

        XSF_HOST_DEVICE inline bool operator<=(const double_double &lhs, const double_double &rhs) {
            if (lhs.hi < rhs.hi) {
                return true;
            }
            if (lhs.hi > rhs.hi) {
                return false;
            }
            return lhs.lo <= rhs.lo;
        }

        XSF_HOST_DEVICE inline bool operator<=(const double_double &lhs, const double rhs) {
            if (lhs.hi < rhs) {
                return true;
            }
            if (lhs.hi > rhs) {
                return false;
            }
            return lhs.lo <= 0.0;
        }

        template <typename T>
        XSF_HOST_DEVICE bool operator>=(const double_double &lhs, const T &rhs) {
            return rhs <= lhs;
        }

        XSF_HOST_DEVICE inline bool operator>=(const double lhs, const double_double &rhs) { return rhs <= lhs; }

        XSF_HOST_DEVICE inline bool operator<=(const double lhs, const double_double &rhs) { return rhs >= lhs; }

        // Math functions

        XSF_HOST_DEVICE inline double_double mul_pwr2(const double_double &lhs, double rhs) {
            /* double-double * double,  where double is a power of 2. */
            return double_double(lhs.hi * rhs, lhs.lo * rhs);
        }

        XSF_HOST_DEVICE inline bool isfinite(const double_double &a) { return std::isfinite(a.hi); }

        XSF_HOST_DEVICE inline bool isinf(const double_double &a) { return std::isinf(a.hi); }

        XSF_HOST_DEVICE inline double_double round(const double_double &a) {
            double hi = two_nint(a.hi);
            double lo;

            if (hi == a.hi) {
                /* High word is an integer already.  Round the low word.*/
                lo = two_nint(a.lo);

                /* Renormalize. This is needed if a.hi = some integer, a.lo = 1/2.*/
                hi = quick_two_sum(hi, lo, &lo);
            } else {
                /* High word is not an integer. */
                lo = 0.0;
                if (std::abs(hi - a.hi) == 0.5 && a.lo < 0.0) {
                    /* There is a tie in the high word, consult the low word
                       to break the tie. */
                    hi -= 1.0; /* NOTE: This does not cause INEXACT. */
                }
            }
            return double_double(hi, lo);
        }

        XSF_HOST_DEVICE inline double_double floor(const double_double &a) {
            double hi = std::floor(a.hi);
            double lo = 0.0;

            if (hi == a.hi) {
                /* High word is integer already.  Round the low word. */
                lo = std::floor(a.lo);
                hi = quick_two_sum(hi, lo, &lo);
            }

            return double_double(hi, lo);
        }

        XSF_HOST_DEVICE inline double_double ceil(const double_double &a) {
            double hi = std::ceil(a.hi);
            double lo = 0.0;

            if (hi == a.hi) {
                /* High word is integer already.  Round the low word. */
                lo = std::ceil(a.lo);
                hi = quick_two_sum(hi, lo, &lo);
            }

            return double_double(hi, lo);
        }

        XSF_HOST_DEVICE inline double_double trunc(const double_double &a) {
            return (a.hi >= 0.0) ? floor(a) : ceil(a);
        }

        XSF_HOST_DEVICE inline double_double abs(const double_double &a) { return (a.hi < 0.0 ? -a : a); }

        XSF_HOST_DEVICE inline double_double fmod(const double_double &lhs, const double_double &rhs) {
            double_double n = trunc(lhs / rhs);
            return lhs - rhs * n;
        }

        XSF_HOST_DEVICE inline double_double remainder(const double_double &lhs, const double_double &rhs) {
            double_double n = round(lhs / rhs);
            return lhs - rhs * n;
        }

        XSF_HOST_DEVICE inline std::pair<double_double, double_double> divrem(const double_double &lhs,
                                                                                  const double_double &rhs) {
            double_double n = round(lhs / rhs);
            double_double remainder = lhs - n * rhs;
            return {n, remainder};
        }

        XSF_HOST_DEVICE inline double_double square(const double_double &a) {
            double p1, p2;
            double s1, s2;
            p1 = two_sqr(a.hi, &p2);
            p2 += 2.0 * a.hi * a.lo;
            p2 += a.lo * a.lo;
            s1 = quick_two_sum(p1, p2, &s2);
            return double_double(s1, s2);
        }

        XSF_HOST_DEVICE inline double_double square(const double a) {
            double p1, p2;
            p1 = two_sqr(a, &p2);
            return double_double(p1, p2);
        }

        XSF_HOST_DEVICE inline double_double ldexp(const double_double &a, int expt) {
            // float128 * (2.0 ^ expt)
            return double_double(std::ldexp(a.hi, expt), std::ldexp(a.lo, expt));
        }

        XSF_HOST_DEVICE inline double_double frexp(const double_double &a, int *expt) {
            //    r"""return b and l s.t. 0.5<=|b|<1 and 2^l == a
            //    0.5<=|b[0]|<1.0 or |b[0]| == 1.0 and b[0]*b[1]<0
            //    """
            int exponent;
            double man = std::frexp(a.hi, &exponent);
            double b1 = std::ldexp(a.lo, -exponent);
            if (std::abs(man) == 0.5 && man * b1 < 0) {
                man *= 2;
                b1 *= 2;
                exponent -= 1;
            }
            *expt = exponent;
            return double_double(man, b1);
        }

        // Numeric limits

        XSF_HOST_DEVICE inline double_double quiet_NaN() {
            return double_double(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        }

        XSF_HOST_DEVICE inline double_double infinity() {
            return double_double(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
        }

        const double_double inv_fact[] = {double_double(1.66666666666666657e-01, 9.25185853854297066e-18),
                                          double_double(4.16666666666666644e-02, 2.31296463463574266e-18),
                                          double_double(8.33333333333333322e-03, 1.15648231731787138e-19),
                                          double_double(1.38888888888888894e-03, -5.30054395437357706e-20),
                                          double_double(1.98412698412698413e-04, 1.72095582934207053e-22),
                                          double_double(2.48015873015873016e-05, 2.15119478667758816e-23),
                                          double_double(2.75573192239858925e-06, -1.85839327404647208e-22),
                                          double_double(2.75573192239858883e-07, 2.37677146222502973e-23),
                                          double_double(2.50521083854417202e-08, -1.44881407093591197e-24),
                                          double_double(2.08767569878681002e-09, -1.20734505911325997e-25),
                                          double_double(1.60590438368216133e-10, 1.25852945887520981e-26),
                                          double_double(1.14707455977297245e-11, 2.06555127528307454e-28),
                                          double_double(7.64716373181981641e-13, 7.03872877733453001e-30),
                                          double_double(4.77947733238738525e-14, 4.39920548583408126e-31),
                                          double_double(2.81145725434552060e-15, 1.65088427308614326e-31)};

        // Math constants
        const double_double E = double_double(2.718281828459045091e+00, 1.445646891729250158e-16);
        const double_double LOG2 = double_double(6.931471805599452862e-01, 2.319046813846299558e-17);
        const double EPS = 4.93038065763132e-32; // 2^-104

        /* Exponential.  Computes exp(x) in double-double precision. */
        XSF_HOST_DEVICE inline double_double exp(const double_double &a) {
            /* Strategy:  We first reduce the size of x by noting that

               exp(kr + m * log(2)) = 2^m * exp(r)^k

               where m and k are integers.  By choosing m appropriately
               we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is
               evaluated using the familiar Taylor series.  Reducing the
               argument substantially speeds up the convergence.       */

            constexpr double k = 512.0;
            constexpr double inv_k = 1.0 / k;
            double m;
            double_double r, s, t, p;
            int i = 0;

            if (a.hi <= -709.0) {
                return double_double(0.0);
            }

            if (a.hi >= 709.0) {
                return infinity();
            }

            if (a == 0.0) {
                return double_double(1.0);
            }

            if (a == 1.0) {
                return E;
            }

            m = std::floor(a.hi / LOG2.hi + 0.5);
            r = mul_pwr2(double_double(a) - LOG2 * m, inv_k);

            p = square(r);
            s = r + mul_pwr2(p, 0.5);
            p = p * r;
            t = p * inv_fact[0];
            do {
                s = s + t;
                p = p * r;
                ++i;
                t = p * inv_fact[i];
            } while ((std::abs(static_cast<double>(t)) > inv_k * EPS) && i < 5);

            s = s + t;

            for (int j = 0; j < 9; j++) {
                s = mul_pwr2(s, 2.0) + square(s);
            }
            s = s + 1.0;

            return ldexp(s, static_cast<int>(m));
        }

        /* Logarithm.  Computes log(x) in double-double precision.
           This is a natural logarithm (i.e., base e).            */
        XSF_HOST_DEVICE inline double_double log(const double_double &a) {
            /* Strategy.  The Taylor series for log converges much more
               slowly than that of exp, due to the lack of the factorial
               term in the denominator.  Hence this routine instead tries
               to determine the root of the function

               f(x) = exp(x) - a

               using Newton iteration.  The iteration is given by

               x' = x - f(x)/f'(x)
               = x - (1 - a * exp(-x))
               = x + a * exp(-x) - 1.

               Only one iteration is needed, since Newton's iteration
               approximately doubles the number of digits per iteration. */
            double_double x;

            if (a == 1.0) {
                return double_double(0.0);
            }

            if (a.hi <= 0.0) {
                return quiet_NaN();
            }

            x = double_double(std::log(a.hi)); /* Initial approximation */

            /* x = x + a * exp(-x) - 1.0; */
            x = x + a * exp(-x) - 1.0;
            return x;
        }

        XSF_HOST_DEVICE inline double_double log1p(const double_double &a) {
            double_double ans;
            double la, elam1, ll;
            if (a.hi <= -1.0) {
                return -infinity();
            }
            la = std::log1p(a.hi);
            elam1 = xsf::cephes::expm1(la);
            ll = std::log1p(a.lo / (1 + a.hi));
            if (a.hi > 0) {
                ll -= (elam1 - a.hi) / (elam1 + 1);
            }
            ans = double_double(la) + ll;
            return ans;
        }
    } // namespace detail

} // namespace cephes
} // namespace xsf
