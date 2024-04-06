#pragma once

#include "../config.h"

#include "unity.h"

namespace special {
namespace dd_real {

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
        SPECFUN_HOST_DEVICE inline double quick_two_sum(double a, double b, double *err) {
            volatile double s = a + b;
            volatile double c = s - a;
            *err = b - c;
            return s;
        }

        /* Computes fl(a+b) and err(a+b).  */
        SPECFUN_HOST_DEVICE inline double two_sum(double a, double b, double *err) {
            volatile double s = a + b;
            volatile double c = s - a;
            volatile double d = b - c;
            volatile double e = s - c;
            *err = (a - e) + d;
            return s;
        }

        /* Computes fl(a*b) and err(a*b). */
        SPECFUN_HOST_DEVICE inline double two_prod(double a, double b, double *err) {
            volatile double p = a * b;
            *err = std::fma(a, b, -p);
            return p;
        }

        /* Computes fl(a*a) and err(a*a).  Faster than the above method. */
        SPECFUN_HOST_DEVICE inline double two_sqr(double a, double *err) {
            volatile double p = a * a;
            *err = std::fma(a, a, -p);
            return p;
        }

        /* Computes the nearest integer to d. */
        SPECFUN_HOST_DEVICE inline double two_nint(double d) {
            if (d == std::floor(d)) {
                return d;
            }
            return std::floor(d + 0.5);
        }

    } // namespace detail

    class DoubleDouble {
      public:
        double hi, lo;

        DoubleDouble() : hi(0.0), lo(0.0) {}
        DoubleDouble(double high, double low) : hi(high), lo(low) {}
        explicit DoubleDouble(double high) : hi(high), lo(0.0) {}

        SPECFUN_HOST_DEVICE inline DoubleDouble operator-() const { return DoubleDouble(-hi, -lo); }

        SPECFUN_HOST_DEVICE inline DoubleDouble operator+(const double rhs) const {
            double s1, s2;
            s1 = detail::two_sum(hi, rhs, &s2);
            s2 += lo;
            s1 = detail::quick_two_sum(s1, s2, &s2);
            return DoubleDouble(s1, s2);
        }

        SPECFUN_HOST_DEVICE inline DoubleDouble operator+(const DoubleDouble &rhs) const {
            /* This one satisfies IEEE style error bound,
               due to K. Briggs and W. Kahan.                   */
            double s1, s2, t1, t2;

            s1 = detail::two_sum(hi, rhs.hi, &s2);
            t1 = detail::two_sum(lo, rhs.lo, &t2);
            s2 += t1;
            s1 = detail::quick_two_sum(s1, s2, &s2);
            s2 += t2;
            s1 = detail::quick_two_sum(s1, s2, &s2);
            return DoubleDouble(s1, s2);
        }

        SPECFUN_HOST_DEVICE inline DoubleDouble operator-(const double rhs) const {
            double s1, s2;
            s1 = detail::two_sum(hi, -rhs, &s2);
            s2 += lo;
            s1 = detail::quick_two_sum(s1, s2, &s2);
            return DoubleDouble(s1, s2);
        }

        SPECFUN_HOST_DEVICE inline DoubleDouble operator-(const DoubleDouble &rhs) const { return *this + (-rhs); }

        SPECFUN_HOST_DEVICE inline DoubleDouble operator*(const double rhs) const {
            double p1, p2, e1, e2;
            p1 = detail::two_prod(hi, rhs, &e1);
            p2 = detail::two_prod(lo, rhs, &e2);
            p1 = detail::quick_two_sum(p1, e2 + p2 + e1, &e1);
            return DoubleDouble(p1, e1);
        }

        SPECFUN_HOST_DEVICE inline DoubleDouble operator*(const DoubleDouble &rhs) const {
            double p1, p2;
            p1 = detail::two_prod(hi, rhs.hi, &p2);
            p2 += (hi * rhs.lo + lo * rhs.hi);
            p1 = detail::quick_two_sum(p1, p2, &p2);
            return DoubleDouble(p1, p2);
        }

        SPECFUN_HOST_DEVICE inline DoubleDouble operator/(const DoubleDouble &rhs) const {
            double q1, q2, q3;
            DoubleDouble r;

            q1 = hi / rhs.hi; /* approximate quotient */

            r = *this - rhs * q1;

            q2 = r.hi / rhs.hi;
            r = r - rhs * q2;

            q3 = r.hi / rhs.hi;

            q1 = detail::quick_two_sum(q1, q2, &q2);
            r = DoubleDouble(q1, q2) + q3;
            return r;
        }

        SPECFUN_HOST_DEVICE inline DoubleDouble operator/(const double rhs) const { return *this / DoubleDouble(rhs); }

        SPECFUN_HOST_DEVICE inline bool operator==(const DoubleDouble &rhs) const {
            return (hi == rhs.hi && lo == rhs.lo);
        }

        SPECFUN_HOST_DEVICE inline bool operator==(const double rhs) const { return (hi == rhs && lo == 0.0); }

        SPECFUN_HOST_DEVICE inline bool operator!=(const DoubleDouble &rhs) const {
            return (hi != rhs.hi) || (lo != rhs.lo);
        }

        SPECFUN_HOST_DEVICE inline bool operator!=(const double rhs) const { return (hi != rhs) || (lo != 0.0); }

        SPECFUN_HOST_DEVICE inline bool operator<(const DoubleDouble &rhs) const {
            if (hi < rhs.hi) {
                return true;
            }
            if (hi > rhs.hi) {
                return false;
            }
            return lo < rhs.lo;
        }

        SPECFUN_HOST_DEVICE inline bool operator<(const double rhs) const {
            if (hi < rhs) {
                return true;
            }
            if (hi > rhs) {
                return false;
            }
            return lo < 0.0;
        }

        SPECFUN_HOST_DEVICE inline bool operator<=(const DoubleDouble &rhs) const {
            if (hi < rhs.hi) {
                return true;
            }
            if (hi > rhs.hi) {
                return false;
            }
            return lo <= rhs.lo;
        }

        SPECFUN_HOST_DEVICE inline bool operator<=(const double rhs) const {
            if (hi < rhs) {
                return true;
            }
            if (hi > rhs) {
                return false;
            }
            return lo <= 0.0;
        }

        template <typename T>
        SPECFUN_HOST_DEVICE inline bool operator>(const T &rhs) const {
            return rhs < *this;
        }

        template <typename T>
        SPECFUN_HOST_DEVICE inline bool operator>=(const T &rhs) const {
            return rhs <= *this;
        }

        SPECFUN_HOST_DEVICE inline explicit operator double() const { return hi; }

        SPECFUN_HOST_DEVICE inline explicit operator int() const { return static_cast<int>(hi); }
    };

    SPECFUN_HOST_DEVICE inline DoubleDouble operator+(const double lhs, const DoubleDouble &rhs) {
        double s1, s2;
        s1 = detail::two_sum(lhs, rhs.hi, &s2);
        s2 += rhs.lo;
        s1 = detail::quick_two_sum(s1, s2, &s2);
        return DoubleDouble(s1, s2);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble operator-(const double lhs, const DoubleDouble &rhs) {
        double s1, s2;
        s1 = detail::two_sum(lhs, -rhs.hi, &s2);
        s2 -= rhs.lo;
        s1 = detail::quick_two_sum(s1, s2, &s2);
        return DoubleDouble(s1, s2);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble operator*(const double lhs, const DoubleDouble &rhs) {
        double p1, p2, e1, e2;
        p1 = detail::two_prod(lhs, rhs.hi, &e1);
        p2 = detail::two_prod(lhs, rhs.lo, &e2);
        p1 = detail::quick_two_sum(p1, e2 + p2 + e1, &e1);
        return DoubleDouble(p1, e1);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble operator/(const double lhs, const DoubleDouble &rhs) {
        return DoubleDouble(lhs) / rhs;
    }

    SPECFUN_HOST_DEVICE inline bool operator==(const double lhs, const DoubleDouble &rhs) {
        return (lhs == rhs.hi) && (rhs.lo == 0.0);
    }

    SPECFUN_HOST_DEVICE inline bool operator!=(const double lhs, const DoubleDouble &rhs) {
        return (rhs.hi != lhs) || (rhs.lo != 0.0);
    }

    SPECFUN_HOST_DEVICE inline bool operator<(const double lhs, const DoubleDouble &rhs) { return rhs > lhs; }

    SPECFUN_HOST_DEVICE inline bool operator<=(const double lhs, const DoubleDouble &rhs) { return rhs >= lhs; }

    SPECFUN_HOST_DEVICE inline bool operator>(const double lhs, const DoubleDouble &rhs) { return rhs < lhs; }

    SPECFUN_HOST_DEVICE inline bool operator>=(const double lhs, const DoubleDouble &rhs) { return rhs <= lhs; }

    SPECFUN_HOST_DEVICE inline DoubleDouble mul_pwr2(const DoubleDouble &lhs, double rhs) {
        /* double-double * double,  where double is a power of 2. */
        return DoubleDouble(lhs.hi * rhs, lhs.lo * rhs);
    }

    SPECFUN_HOST_DEVICE inline bool isfinite(const DoubleDouble &a) { return std::isfinite(a.hi); }

    SPECFUN_HOST_DEVICE inline bool isinf(const DoubleDouble &a) { return std::isinf(a.hi); }

    SPECFUN_HOST_DEVICE inline DoubleDouble round(const DoubleDouble &a) {
        double hi = detail::two_nint(a.hi);
        double lo;

        if (hi == a.hi) {
            /* High word is an integer already.  Round the low word.*/
            lo = detail::two_nint(a.lo);

            /* Renormalize. This is needed if a.hi = some integer, a.lo = 1/2.*/
            hi = detail::quick_two_sum(hi, lo, &lo);
        } else {
            /* High word is not an integer. */
            lo = 0.0;
            if (std::abs(hi - a.hi) == 0.5 && a.lo < 0.0) {
                /* There is a tie in the high word, consult the low word
                   to break the tie. */
                hi -= 1.0; /* NOTE: This does not cause INEXACT. */
            }
        }
        return DoubleDouble(hi, lo);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble floor(const DoubleDouble &a) {
        double hi = std::floor(a.hi);
        double lo = 0.0;

        if (hi == a.hi) {
            /* High word is integer already.  Round the low word. */
            lo = std::floor(a.lo);
            hi = detail::quick_two_sum(hi, lo, &lo);
        }

        return DoubleDouble(hi, lo);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble ceil(const DoubleDouble &a) {
        double hi = std::ceil(a.hi);
        double lo = 0.0;

        if (hi == a.hi) {
            /* High word is integer already.  Round the low word. */
            lo = std::ceil(a.lo);
            hi = detail::quick_two_sum(hi, lo, &lo);
        }

        return DoubleDouble(hi, lo);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble trunc(const DoubleDouble &a) { return (a.hi >= 0.0) ? floor(a) : ceil(a); }

    SPECFUN_HOST_DEVICE inline DoubleDouble abs(const DoubleDouble &a) { return (a.hi < 0.0 ? -a : a); }

    SPECFUN_HOST_DEVICE inline DoubleDouble fmod(const DoubleDouble &lhs, const DoubleDouble &rhs) {
        DoubleDouble n = trunc(lhs / rhs);
        return lhs - rhs * n;
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble remainder(const DoubleDouble &lhs, const DoubleDouble &rhs) {
        DoubleDouble n = round(lhs / rhs);
        return lhs - rhs * n;
    }

    SPECFUN_HOST_DEVICE inline std::pair<DoubleDouble, DoubleDouble> divrem(const DoubleDouble &lhs,
                                                                            const DoubleDouble &rhs) {
        DoubleDouble n = round(lhs / rhs);
        DoubleDouble remainder = lhs - n * rhs;
        return {n, remainder};
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble square(const DoubleDouble &a) {
        double p1, p2;
        double s1, s2;
        p1 = detail::two_sqr(a.hi, &p2);
        p2 += 2.0 * a.hi * a.lo;
        p2 += a.lo * a.lo;
        s1 = detail::quick_two_sum(p1, p2, &s2);
        return DoubleDouble(s1, s2);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble square(const double a) {
        double p1, p2;
        p1 = detail::two_sqr(a, &p2);
        return DoubleDouble(p1, p2);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble ldexp(const DoubleDouble &a, int expt) {
        // float128 * (2.0 ^ expt)
        return DoubleDouble(std::ldexp(a.hi, expt), std::ldexp(a.lo, expt));
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble frexp(const DoubleDouble &a, int *expt) {
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
        return DoubleDouble(man, b1);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble quiet_NaN() {
        return DoubleDouble(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble infinity() {
        return DoubleDouble(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
    }

    namespace detail {

        const DoubleDouble inv_fact[] = {DoubleDouble(1.66666666666666657e-01, 9.25185853854297066e-18),
                                         DoubleDouble(4.16666666666666644e-02, 2.31296463463574266e-18),
                                         DoubleDouble(8.33333333333333322e-03, 1.15648231731787138e-19),
                                         DoubleDouble(1.38888888888888894e-03, -5.30054395437357706e-20),
                                         DoubleDouble(1.98412698412698413e-04, 1.72095582934207053e-22),
                                         DoubleDouble(2.48015873015873016e-05, 2.15119478667758816e-23),
                                         DoubleDouble(2.75573192239858925e-06, -1.85839327404647208e-22),
                                         DoubleDouble(2.75573192239858883e-07, 2.37677146222502973e-23),
                                         DoubleDouble(2.50521083854417202e-08, -1.44881407093591197e-24),
                                         DoubleDouble(2.08767569878681002e-09, -1.20734505911325997e-25),
                                         DoubleDouble(1.60590438368216133e-10, 1.25852945887520981e-26),
                                         DoubleDouble(1.14707455977297245e-11, 2.06555127528307454e-28),
                                         DoubleDouble(7.64716373181981641e-13, 7.03872877733453001e-30),
                                         DoubleDouble(4.77947733238738525e-14, 4.39920548583408126e-31),
                                         DoubleDouble(2.81145725434552060e-15, 1.65088427308614326e-31)};

        // Math constants
        const DoubleDouble E = DoubleDouble(2.718281828459045091e+00, 1.445646891729250158e-16);
        const DoubleDouble LOG2 = DoubleDouble(6.931471805599452862e-01, 2.319046813846299558e-17);
        const double EPS = 4.93038065763132e-32; // 2^-104
    }                                            // namespace detail

    /* Exponential.  Computes exp(x) in double-double precision. */
    SPECFUN_HOST_DEVICE inline DoubleDouble exp(const DoubleDouble &a) {
        /* Strategy:  We first reduce the size of x by noting that

           exp(kr + m * log(2)) = 2^m * exp(r)^k

           where m and k are integers.  By choosing m appropriately
           we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is
           evaluated using the familiar Taylor series.  Reducing the
           argument substantially speeds up the convergence.       */

        constexpr double k = 512.0;
        constexpr double inv_k = 1.0 / k;
        double m;
        DoubleDouble r, s, t, p;
        int i = 0;

        if (a.hi <= -709.0) {
            return DoubleDouble(0.0);
        }

        if (a.hi >= 709.0) {
            return infinity();
        }

        if (a == 0.0) {
            return DoubleDouble(1.0);
        }

        if (a == 1.0) {
            return detail::E;
        }

        m = std::floor(a.hi / detail::LOG2.hi + 0.5);
        r = mul_pwr2(DoubleDouble(a) - detail::LOG2 * m, inv_k);

        p = square(r);
        s = r + mul_pwr2(p, 0.5);
        p = p * r;
        t = p * detail::inv_fact[0];
        do {
            s = s + t;
            p = p * r;
            ++i;
            t = p * detail::inv_fact[i];
        } while ((std::abs(static_cast<double>(t)) > inv_k * detail::EPS) && i < 5);

        s = s + t;

        for (int j = 0; j < 9; j++) {
            s = mul_pwr2(s, 2.0) + square(s);
        }
        s = s + 1.0;

        return ldexp(s, static_cast<int>(m));
    }

    /* Logarithm.  Computes log(x) in double-double precision.
       This is a natural logarithm (i.e., base e).            */
    SPECFUN_HOST_DEVICE inline DoubleDouble log(const DoubleDouble &a) {
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
        DoubleDouble x;

        if (a == 1.0) {
            return DoubleDouble(0.0);
        }

        if (a.hi <= 0.0) {
            return quiet_NaN();
        }

        x = DoubleDouble(std::log(a.hi)); /* Initial approximation */

        /* x = x + a * exp(-x) - 1.0; */
        x = x + a * exp(-x) - 1.0;
        return x;
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble log1p(const DoubleDouble &a) {
        DoubleDouble ans;
        double la, elam1, ll;
        if (a.hi <= -1.0) {
            return -infinity();
        }
        la = std::log1p(a.hi);
        elam1 = special::cephes::expm1(la);
        ll = std::log1p(a.lo / (1 + a.hi));
        if (a.hi > 0) {
            ll -= (elam1 - a.hi) / (elam1 + 1);
        }
        ans = DoubleDouble(la) + ll;
        return ans;
    }

} // namespace dd_real
} // namespace special
