#pragma once


#include "../config.h"


namespace special {
namespace dd_real {

    namespace detail {

        constexpr double __DD_SPLITTER = 134217729.0; // = 2^27 + 1
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

        /* Computes fl(a-b) and err(a-b).  Assumes |a| >= |b| */
        SPECFUN_HOST_DEVICE inline double quick_two_diff(double a, double b, double *err) {
            volatile double s = a - b;
            volatile double c = a - s;
            *err = c - b;
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

        /* Computes fl(a-b) and err(a-b).  */
        SPECFUN_HOST_DEVICE inline double two_diff(double a, double b, double *err) {
            volatile double s = a - b;
            volatile double c = s - a;
            volatile double d = b + c;
            volatile double e = s - c;
            *err = (a - e) - d;
            return s;
        }

        /* Computes high word and lo word of a */
        SPECFUN_HOST_DEVICE inline void two_split(double a, double *hi, double *lo) {
            volatile double temp, tempma;
            if (a > __DD_SPLIT_THRESH || a < -__DD_SPLIT_THRESH) {
                a *= 3.7252902984619140625e-09; // 2^-28
                temp = __DD_SPLITTER * a;
                tempma = temp - a;
                *hi = temp - tempma;
                *lo = a - *hi;
                *hi *= 268435456.0; // 2^28
                *lo *= 268435456.0; // 2^28
            } else {
                temp = __DD_SPLITTER * a;
                tempma = temp - a;
                *hi = temp - tempma;
                *lo = a - *hi;
            }
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

        SPECFUN_HOST_DEVICE inline double two_div(double a, double b, double *err) {
            volatile double q1, q2;
            double p1, p2;
            double s, e;

            q1 = a / b;

            /* Compute  a - q1 * b */
            p1 = two_prod(q1, b, &p2);
            s = two_diff(a, p1, &e);
            e -= p2;

            /* get next approximation */
            q2 = (s + e) / b;

            return quick_two_sum(q1, q2, err);
        }

        /* Computes the nearest integer to d. */
        SPECFUN_HOST_DEVICE inline double two_nint(double d) {
            if (d == std::floor(d)) {
                return d;
            }
            return std::floor(d + 0.5);
        }

        /* Computes the truncated integer. */
        SPECFUN_HOST_DEVICE inline double two_aint(double d) {
            return (d >= 0.0 ? std::floor(d) : std::ceil(d));
        }


        /* Compare a and b */
        SPECFUN_HOST_DEVICE inline int two_comp(const double a, const double b) {
            /* Works for non-NAN inputs */
            return (a < b ? -1 : (a > b ? 1 : 0));
        }

    }

    class DoubleDouble {
    public:
        double hi, lo;

        DoubleDouble() : hi(0.0), lo(0.0) {}
        DoubleDouble(double high, double low) : hi(high), lo(low) {}
        explicit DoubleDouble(double high) : hi(high), lo(0.0) {}

        SPECFUN_HOST_DEVICE inline DoubleDouble operator-() const {
            return DoubleDouble(-hi, -lo);
        }

        SPECFUN_HOST_DEVICE inline DoubleDouble operator+(const double rhs) const {
            double s1, s2;
            s1 = detail::two_sum(hi, rhs, &s2);
            s2 += lo;
            s1 = detail::quick_two_sum(s1, s2, &s2);
            return DoubleDouble(s1, s2);
        }

        SPECFUN_HOST_DEVICE inline DoubleDouble operator+(const DoubleDouble& rhs) const {
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

        SPECFUN_HOST_DEVICE inline DoubleDouble operator-(const DoubleDouble& rhs) const {
            return *this + (-rhs);
        }

        SPECFUN_HOST_DEVICE inline DoubleDouble operator*(const double rhs) const {
            double p1, p2, e1, e2;
            p1 = detail::two_prod(hi, rhs, &e1);
            p2 = detail::two_prod(lo, rhs, &e2);
            p1 = detail::quick_two_sum(p1, e2 + p2 + e1, &e1);
            return DoubleDouble(p1, e1);
        }

        SPECFUN_HOST_DEVICE inline DoubleDouble operator*(const DoubleDouble& rhs) const {
            double p1, p2;
            p1 = detail::two_prod(hi, rhs.hi, &p2);
            p2 += (hi * rhs.lo + lo * rhs.hi);
            p1 = detail::quick_two_sum(p1, p2, &p2);
            return DoubleDouble(p1, p2);
        }

        SPECFUN_HOST_DEVICE inline DoubleDouble operator/(const DoubleDouble& rhs) const {
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

        SPECFUN_HOST_DEVICE inline DoubleDouble operator/(const double rhs) const {
            return *this / DoubleDouble(rhs);
        }

        SPECFUN_HOST_DEVICE inline bool operator==(const DoubleDouble& rhs) const {
            return (hi == rhs.hi && lo == rhs.lo);

        }

        SPECFUN_HOST_DEVICE inline bool operator==(const double rhs) const {
            return (hi == rhs && lo == 0.0);
        }

        SPECFUN_HOST_DEVICE inline bool operator!=(const DoubleDouble& rhs) const {
            return (hi != rhs.hi) || (lo != rhs.lo);
        }

        SPECFUN_HOST_DEVICE inline bool operator!=(const double rhs) const {
            return (hi != rhs) || (lo != 0.0);
        }

        SPECFUN_HOST_DEVICE inline bool operator<(const DoubleDouble& rhs) const {
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

        SPECFUN_HOST_DEVICE inline bool operator<=(const DoubleDouble& rhs) const {
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


        template<typename T>
        SPECFUN_HOST_DEVICE inline bool operator>(const T& rhs) const {
            return rhs < *this;
        }

        template<typename T>
        SPECFUN_HOST_DEVICE inline bool operator>=(const T& rhs) const {
            return rhs <= *this;
        }

        SPECFUN_HOST_DEVICE inline explicit operator double() const {
            return hi;
        }

        SPECFUN_HOST_DEVICE inline explicit operator int() const {
            return static_cast<int>(hi);
        }
    };

    SPECFUN_HOST_DEVICE inline DoubleDouble operator+(const double lhs, const DoubleDouble& rhs) {
        double s1, s2;
        s1 = detail::two_sum(lhs, rhs.hi, &s2);
        s2 += rhs.lo;
        s1 = detail::quick_two_sum(s1, s2, &s2);
        return DoubleDouble(s1, s2);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble operator-(const double lhs, const DoubleDouble& rhs) {
        double s1, s2;
        s1 = detail::two_sum(lhs, -rhs.hi, &s2);
        s2 -= rhs.lo;
        s1 = detail::quick_two_sum(s1, s2, &s2);
        return DoubleDouble(s1, s2);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble operator*(const double lhs, const DoubleDouble& rhs) {
        double p1, p2, e1, e2;
        p1 = detail::two_prod(lhs, rhs.hi, &e1);
        p2 = detail::two_prod(lhs, rhs.lo, &e2);
        p1 = detail::quick_two_sum(p1, e2 + p2 + e1, &e1);
        return DoubleDouble(p1, e1);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble operator/(const double lhs, const DoubleDouble& rhs) {
        return DoubleDouble(lhs) / rhs;
    }

    SPECFUN_HOST_DEVICE inline bool operator==(const double lhs, const DoubleDouble& rhs) {
        return (lhs == rhs.hi) && (rhs.lo == 0.0);
    }

    SPECFUN_HOST_DEVICE inline bool operator!=(const double lhs, const DoubleDouble& rhs) {
        return (rhs.hi != lhs) || (rhs.lo != 0.0);
    }

    SPECFUN_HOST_DEVICE inline bool operator<(const double lhs, const DoubleDouble& rhs) {
        return rhs > lhs;
    }

    SPECFUN_HOST_DEVICE inline bool operator<=(const double lhs, const DoubleDouble& rhs) {
        return rhs >= lhs;
    }

    SPECFUN_HOST_DEVICE inline bool operator>(const double lhs, const DoubleDouble& rhs) {
        return rhs < lhs;
    }

    SPECFUN_HOST_DEVICE inline bool operator>=(const double lhs, const DoubleDouble& rhs) {
        return rhs <= lhs;
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble add(const double lhs, const double rhs) {
        double s, e;
        s = detail::two_sum(lhs, rhs, &e);
        return DoubleDouble(s, e);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble sub(const double lhs, const double rhs) {
        double s, e;
        s = detail::two_diff(lhs, rhs, &e);
        return DoubleDouble(s, e);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble mul(const double lhs, const double rhs) {
        double p, e;
        p = detail::two_prod(lhs, rhs, &e);
        return DoubleDouble(p, e);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble mul_pwr2(const DoubleDouble& lhs, double rhs) {
        /* double-double * double,  where double is a power of 2. */
        return DoubleDouble(lhs.hi * rhs, lhs.lo * rhs);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble div(const double lhs, const double rhs) {
        return DoubleDouble(lhs) / DoubleDouble(rhs);
    }


    /* Compare a and b */
    SPECFUN_HOST_DEVICE inline int comp(const double a, const double b) {
        /* Works for non-NAN inputs */
        return (a < b ? -1 : (a > b ? 1 : 0));
    }

    SPECFUN_HOST_DEVICE inline int comp(const DoubleDouble& a, const DoubleDouble& b) {
        int cmp = comp(a.hi, b.hi);
        if (cmp == 0) {
            cmp = comp(a.lo, b.lo);
        }
        return cmp;
    }

    SPECFUN_HOST_DEVICE inline int comp(const DoubleDouble& a, const double b) {
        int cmp = comp(a.hi, b);
        if (cmp == 0) {
            cmp = comp(a.lo, 0);
        }
        return cmp;
    }

    SPECFUN_HOST_DEVICE inline int comp(const double a, const DoubleDouble& b) {
        int cmp = comp(a, b.hi);
        if (cmp == 0) {
            cmp = comp(0.0, b.lo);
        }
        return cmp;
    }

    SPECFUN_HOST_DEVICE inline bool isfinite(const DoubleDouble& a) {
        return std::isfinite(a.hi);
    }

    SPECFUN_HOST_DEVICE inline bool isinf(const DoubleDouble& a) {
        return std::isinf(a.hi);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble round(const DoubleDouble& a) {
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

    SPECFUN_HOST_DEVICE inline DoubleDouble floor(const DoubleDouble& a) {
        double hi = std::floor(a.hi);
        double lo = 0.0;

        if (hi == a.hi) {
            /* High word is integer already.  Round the low word. */
            lo = std::floor(a.hi);
            hi = detail::quick_two_sum(hi, lo, &lo);
        }

        return DoubleDouble(hi, lo);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble ceil(const DoubleDouble& a) {
        double hi = std::ceil(a.hi);
        double lo = 0.0;

        if (hi == a.hi) {
            /* High word is integer already.  Round the low word. */
            lo = std::ceil(a.lo);
            hi = detail::quick_two_sum(hi, lo, &lo);
        }

        return DoubleDouble(hi, lo);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble trunc(const DoubleDouble& a) {
        return (a.hi >= 0.0) ? floor(a) : ceil(a);
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble abs(const DoubleDouble& a) {
        return (a.hi < 0.0 ? -a : a);
    }
        

    SPECFUN_HOST_DEVICE inline DoubleDouble fmod(const DoubleDouble& lhs, const DoubleDouble& rhs) {
        DoubleDouble n = trunc(lhs / rhs);
        return lhs - rhs * n;
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble remainder(const DoubleDouble& lhs, const DoubleDouble& rhs) {
        DoubleDouble n = round(lhs / rhs);
        return lhs - rhs * n;
    }
    
    SPECFUN_HOST_DEVICE inline std::pair<DoubleDouble, DoubleDouble> divrem(const DoubleDouble& lhs, const DoubleDouble& rhs) {
        DoubleDouble n = round(lhs / rhs);
        DoubleDouble remainder = lhs - n * rhs;
        return {n, remainder};
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble square(const DoubleDouble& a) {
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

    SPECFUN_HOST_DEVICE inline DoubleDouble ldexp(const DoubleDouble& a, int expt) {
        // float128 * (2.0 ^ expt)
        return DoubleDouble(std::ldexp(a.hi, expt), std::ldexp(a.lo, expt));
    }

    SPECFUN_HOST_DEVICE inline DoubleDouble frexp(const DoubleDouble& a, int *expt) {
        //    r"""return b and l s.t. 0.5<=|b|<1 and 2^l == a
        //    0.5<=|b[0]|<1.0 or |b[0]| == 1.0 and b[0]*b[1]<0
        //    """
        int exponent;
        double man = std::frexp(a.hi, &exponent);
        double b1 = std::ldexp(a.lo, -exponent);
        if (std::abs(man) == 0.5 && man * b1 < 0) {
            man *=2;
            b1 *= 2;
            exponent -= 1;
        }
        *expt = exponent;
        return DoubleDouble(man, b1);
    }

}
}
