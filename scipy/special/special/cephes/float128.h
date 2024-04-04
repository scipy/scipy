#pragma once


#include "../config.h"


namespace special {
namespace float128 {

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
            if (a > _DD_SPLIT_THRESH || a < -_DD_SPLIT_THRESH) {
                a *= 3.7252902984619140625e-09; // 2^-28
                temp = _DD_SPLITTER * a;
                tempma = temp - a;
                *hi = temp - tempma;
                *lo = a - *hi;
                *hi *= 268435456.0; // 2^28
                *lo *= 268435456.0; // 2^28
            } else {
                temp = _DD_SPLITTER * a;
                tempma = temp - a;
                *hi = temp - tempma;
                *lo = a - *hi;
            }
        }

        /* Computes fl(a*b) and err(a*b). */
        SPECFUN_HOST_DEVICE inline double two_prod(double a, double b, double *err) {
            volatile double p = a * b;
            *err = DD_FMS(a, b, p);
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

    SPECFUN_HOST_DEVICE inline int comp(const Float128& a, const Float128& b);
    SPECFUN_HOST_DEVICE inline int comp(const double a, const Float128& b);
    SPECFUN_HOST_DEVICE inline int comp(const Float128& a, const double b);
        

    class Float128 {
    public:
        double hi, lo;

    Float128() : hi(0.0), lo(0.0) {}
    Float128(double high, double low) : hi(high), lo(low) {}
        explicit Float128(double high) : hi(high), lo(0.0) {}

        SPECFUN_HOST_DEVICE inline Float128 operator-() const {
            return Float128(-hi, -lo);
        }

        SPECFUN_HOST_DEVICE inline Float128 operator+(const double rhs) const {
            double s1, s2;
            s1 = detail::two_sum(hi, rhs, &s2);
            s2 += lo;
            s1 = detail::quick_two_sum(s1, s2, &s2);
            return Float128(s1, s2);
        }

        SPECFUN_HOST_DEVICE inline Float128 operator+(const Float128& rhs) const {
            /* This one satisfies IEEE style error bound,
               due to K. Briggs and W. Kahan.                   */
            double s1, s2, t1, t2;

            s1 = detail::two_sum(hi, rhs.hi, &s2);
            t1 = detail::two_sum(lo, rhs.lo, &t2);
            s2 += t1;
            s1 = detail::quick_two_sum(s1, s2, &s2);
            s2 += t2;
            s1 = detail::quick_two_sum(s1, s2, &s2);
            return Float128(s1, s2);
        }

        SPECFUN_HOST_DEVICE inline Float128 operator-(const double rhs) const {
            double s1, s2;
            s1 = detail::two_sum(hi, -rhs, &s2);
            s2 += lo;
            s1 = quick_two_sum(s1, s2, &s2);
            return Float128(s1, s2);
        }

        SPECFUN_HOST_DEVICE inline Float128 operator-(const Float128& rhs) const {
            return *this + (-rhs);
        }

        SPECFUN_HOST_DEVICE inline Float128 operator*(const double rhs) const {
            double p1, p2, e1, e2;
            p1 = detail::two_prod(hi, rhs, &e1);
            p2 = detail::two_prod(lo, rhs, &e2);
            p1 = detail::quick_two_sum(p1, e2 + p2 + e1, &e1);
            return Float128(p1, e1);
        }

        SPECFUN_HOST_DEVICE inline Float128 operator*(const Float128& rhs) const {
            double p1, p2;
            p1 = detail::two_prod(hi, rhs.hi, &p2);
            p2 += (hi * rhs.lo + lo * rhs.hi);
            p1 = detail::quick_two_sum(p1, p2, &p2);
            return Float128(p1, p2);
        }

        SPECFUN_HOST_DEVICE inline Float128 operator/(const Float128& rhs) const {
            double q1, q2, q3;
            Float128 r;

            q1 = hi / rhs.hi; /* approximate quotient */

            r = *this - rhs * q1;

            q2 = r.hi / rhs.hi;
            r = r - rhs * q2;

            q3 = r.hi / rhs.hi;

            q1 = detail::quick_two_sum(q1, q2, &q2);
            r = Float128(q1, q2) + q3;
            return r;
        }

        SPECFUN_HOST_DEVICE inline Float128 operator/(const double rhs) const {
            return *this / Float128(rhs);
        }

        SPECFUN_HOST_DEVICE inline bool operator==(const Float128& rhs) const {
            return (hi == rhs.hi && lo == rhs.lo);

        }

        SPECFUN_HOST_DEVICE inline bool operator==(const double rhs) const {
            return (hi == rhs && lo == 0.0);
        }

        SPECFUN_HOST_DEVICE inline bool operator<(const Float128& rhs) const {
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

        SPECFUN_HOST_DEVICE inline bool operator<=(const Float128& rhs) const {
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
    };

    SPECFUN_HOST_DEVICE inline Float128 operator+(const double lhs, const Float128& rhs) {
        double s1, s2;
        s1 = detail::two_sum(lhs, rhs.hi, &s2);
        s2 += rhs.lo;
        s1 = detail::quick_two_sum(s1, s2, &s2);
        return Float128(s1, s2);
    }

    SPECFUN_HOST_DEVICE inline Float128 operator-(const double lhs, const Float128& rhs) {
        double s1, s2;
        s1 = detail::two_sum(lhs, -rhs.hi, &s2);
        s2 -= rhs.lo;
        s1 = quick_two_sum(s1, s2, &s2);
        return Float128(s1, s2);
    }

    SPECFUN_HOST_DEVICE inline Float128 operator*(const double lhs, const Float128& rhs) {
        double p1, p2, e1, e2;
        p1 = detail::two_prod(lhs, rhs.hi, &e1);
        p2 = two_prod(lhs, rhs.lo, &e2);
        p1 = quick_two_sum(p1, e2 + p2 + e1, &e1);
        return Float128(p1, e1);
    }

    SPECFUN_HOST_DEVICE inline Float128 operator/(const double lhs, const Float128& rhs) {
        return Float128(lhs) / rhs;
    }

    SPECFUN_HOST_DEVICE inline bool operator==(const double lhs, const Float128& rhs) {
        return (lhs == rhs.hi) && (rhs.lo == 0.0);

    SPECFUN_HOST_DEVICE inline bool operator<(const double lhs, const Float128& rhs) {
        return rhs > lhs;
    }

    SPECFUN_HOST_DEVICE inline bool operator<=(const double lhs, const Float128& rhs) {
        return rhs >= lhs;
    }

    SPECFUN_HOST_DEVICE inline bool operator>(const double lhs, const Float128& rhs) {
        return rhs < lhs;
    }

    SPECFUN_HOST_DEVICE inline bool operator>=(const double lhs, const Float128& rhs) {
        return rhs <= lhs;
    }

    SPECFUN_HOST_DEVICE inline Float128 add(const double lhs, const double rhs) {
        double s, e;
        s = detail::two_sum(lhs, rhs, &e);
        return Float128(s, e);
    }

    SPECFUN_HOST_DEVICE inline Float128 sub(const double lhs, const double rhs) {
        double s, e;
        s = detail::two_diff(lhs, rhs, &e);
        return Float128(s, e);
    }

    SPECFUN_HOST_DEVICE inline Float128 mul(const double lhs, const double rhs) {
        double p, e;
        p = two_prod(lhs, rhs, &e);
        return Float128(p, e);
    }

    SPECFUN_HOST_DEVICE inline Float128 mul_pwr2(const Float128& lhs, double rhs) {
        /* double-double * double,  where double is a power of 2. */
        return Float128(lhs.hi * rhs, lhs.lo * rhs);
    }

    SPECFUN_HOST_DEVICE inline Float128 div(const double lhs, const double rhs) {
        return Float128(lhs) / Float128(rhs);
    }


    /* Compare a and b */
    SPECFUN_HOST_DEVICE inline int comp(const double a, const double b) {
        /* Works for non-NAN inputs */
        return (a < b ? -1 : (a > b ? 1 : 0));
    }

    SPECFUN_HOST_DEVICE inline int comp(const Float128& a, const Float128& b) {
        int cmp = comp(a.hi, b.hi);
        if (cmp == 0) {
            cmp = comp(a.lo, b.lo);
        }
        return cmp;
    }

    SPECFUN_HOST_DEVICE inline int comp(const Float128& a, const double b) {
        int cmp = comp(a.hi, b);
        if (cmp == 0) {
            cmp = comp(a.lo, 0);
        }
        return cmp;
    }

    SPECFUN_HOST_DEVICE inline int comp(const double a, const Float128& b) {
        int cmp = comp(a, b.hi);
        if (cmp == 0) {
            cmp = comp(0.0, b.lo);
        }
        return cmp;
    }

    SPECFUN_HOST_DEVICE inline bool isfinite(const Float128& a) {
        return std::isfinite(a.hi);
    }

    SPECFUN_HOST_DEVICE inline bool isinf(const Float128& a) {
        return std::isinf(a.hi);
    }

    SPECFUN_HOST_DEVICE inline Float128 round(const Float128& a) {
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
        return Float128(hi, lo);
    }

    SPECFUN_HOST_DEVICE inline Float128 floor(const Float128& a) {
        double hi = std::floor(a.hi);
        double lo = 0.0;

        if (hi == a.hi) {
            /* High word is integer already.  Round the low word. */
            lo = std::floor(a.hi);
            hi = detail::quick_two_sum(hi, lo, &lo);
        }

        return Float128(hi, lo);
    }

    SPECFUN_HOST_DEVICE inline Float128 ceil(const Float128& a) {
        double hi = std::ceil(a.hi);
        double lo = 0.0;

        if (hi == a.hi) {
            /* High word is integer already.  Round the low word. */
            lo = std::ceil(a.lo);
            hi = detail::quick_two_sum(hi, lo, &lo);
        }

        return Float128(hi, lo);
    }

    SPECFUN_HOST_DEVICE inline Float128 trunc(const Float128& a) {
        return (a.hi >= 0.0) ? floor(a) : ceil(a);
    }

    SPECFUN_HOST_DEVICE inline Float128 abs(const Float128& a) {
        return (a.hi < 0.0 ? -a : a);
    }
        

    SPECFUN_HOST_DEVICE inline Float128 fmod(const Float128& lhs, const Float128& rhs) {
        double2 n = trunc(lhs, rhs);
        return lhs - rhs * n;
    }

    SPECFUN_HOST_DEVICE inline Float128 remainder(const Float128& lhs, const Float128& rhs) {
        double2 n = round(lhs / rhs);
        return lhs - rhs * n;
    }
    
    SPECFUN_HOST_DEVICE inline std::pair<Float128, Float128> divrem(const Float128& lhs, const Float128& rhs) {
        Float128 n = round(lhs / rhs);
        Float128 remainder = lhs - n * rhs;
        return {n, remainder};
    }

    SPECFUN_HOST_DEVICE inline Float128 square(const Float128& a) {
        double p1, p2;
        double s1, s2;
        p1 = detail::two_sqr(a.hi, &p2);
        p2 += 2.0 * a.hi * a.lo;
        p2 += a.lo * a.lo;
        s1 = detail::quick_two_sum(p1, p2, &s2);
        return Float128(s1, s2);
    }

    SPECFUN_HOST_DEVICE inline Float128 square(const double a) {
        double p1, p2;
        p1 = detail::two_sqr(a, &p2);
        return Float128(p1, p2);
    }

    SPECFUN_HOST_DEVICE inline Float128 ldexp(const Float128& a, int expt) {
        // float128 * (2.0 ^ expt)
        return Float128(std::ldexp(a.hi, expt), ldexp(a.lo, expt));
    }

    SPECFUN_HOST_DEVICE inline Float128 frexp(const Float128&, int *expt) {
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
        return Float128(man, b1);
    }

}
}
