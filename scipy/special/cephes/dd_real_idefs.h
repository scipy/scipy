/*
 * include/dd_inline.h
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
 * Contains small functions (suitable for inlining) in the double-double
 * arithmetic package.
 */

#ifndef  _DD_REAL_IDEFS_H_
#define  _DD_REAL_IDEFS_H_ 1

#include <float.h>
#include <limits.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "dd_idefs.h"

/*
 ************************************************************************
  Now for the double2 routines
 ************************************************************************
*/

static NPY_INLINE double
dd_hi(const double2 a)
{
    return a.x[0];
}

static NPY_INLINE double
dd_lo(const double2 a)
{
    return a.x[1];
}

static NPY_INLINE int
dd_isnan(const double2 a)
{
    return (DD_ISNAN(a.x[0]) || DD_ISNAN(a.x[1]));
}

static NPY_INLINE int
dd_isfinite(const double2 a)
{
    return DD_ISFINITE(a.x[0]);
}

static NPY_INLINE int
dd_isinf(const double2 a)
{
    return DD_ISINF(a.x[0]);
}

static NPY_INLINE int
dd_is_zero(const double2 a)
{
    return (a.x[0] == 0.0);
}

static NPY_INLINE int
dd_is_one(const double2 a)
{
    return (a.x[0] == 1.0 && a.x[1] == 0.0);
}

static NPY_INLINE int
dd_is_positive(const double2 a)
{
    return (a.x[0] > 0.0);
}

static NPY_INLINE int
dd_is_negative(const double2 a)
{
    return (a.x[0] < 0.0);
}

/* Cast to double. */
static NPY_INLINE double
dd_to_double(const double2 a)
{
    return a.x[0];
}

/* Cast to int. */
static NPY_INLINE int
dd_to_int(const double2 a)
{
    return DD_STATIC_CAST(int, a.x[0]);
}

/*********** Equality and Other Comparisons ************/
static NPY_INLINE int
dd_comp(const double2 a, const double2 b)
{
    int cmp = two_comp(a.x[0], b.x[0]);
    if (cmp == 0) {
        cmp = two_comp(a.x[1], b.x[1]);
    }
    return cmp;
}

static NPY_INLINE int
dd_comp_dd_d(const double2 a, double b)
{
    int cmp = two_comp(a.x[0], b);
    if (cmp == 0) {
        cmp = two_comp(a.x[1], 0);
    }
    return cmp;
}

static NPY_INLINE int
dd_comp_d_dd(double a, const double2 b)
{
    int cmp = two_comp(a, b.x[0]);
    if (cmp == 0) {
        cmp = two_comp(0.0, b.x[1]);
    }
    return cmp;
}


/*********** Creation ************/
static NPY_INLINE double2
dd_create(double hi, double lo)
{
    double2 ret = {{hi, lo}};
    return ret;
}

static NPY_INLINE double2
dd_zero(void)
{
    return DD_C_ZERO;
}

static NPY_INLINE double2
dd_create_d(double hi)
{
    double2 ret = {{hi, 0.0}};
    return ret;
}

static NPY_INLINE double2
dd_create_i(int hi)
{
    double2 ret = {{DD_STATIC_CAST(double, hi), 0.0}};
    return ret;
}

static NPY_INLINE double2
dd_create_dp(const double *d)
{
    double2 ret = {{d[0], d[1]}};
    return ret;
}


/*********** Unary Minus ***********/
static NPY_INLINE double2
dd_neg(const double2 a)
{
    double2 ret = {{-a.x[0], -a.x[1]}};
    return ret;
}

/*********** Rounding ************/
/* Round to Nearest integer */
static NPY_INLINE double2
dd_nint(const double2 a)
{
    double hi = two_nint(a.x[0]);
    double lo;

    if (hi == a.x[0]) {
        /* High word is an integer already.  Round the low word.*/
        lo = two_nint(a.x[1]);

        /* Renormalize. This is needed if x[0] = some integer, x[1] = 1/2.*/
        hi = quick_two_sum(hi, lo, &lo);
    }
    else {
        /* High word is not an integer. */
        lo = 0.0;
        if (fabs(hi - a.x[0]) == 0.5 && a.x[1] < 0.0) {
            /* There is a tie in the high word, consult the low word
               to break the tie. */
            hi -= 1.0; /* NOTE: This does not cause INEXACT. */
        }
    }

    return dd_create(hi, lo);
}

static NPY_INLINE double2
dd_floor(const double2 a)
{
    double hi = floor(a.x[0]);
    double lo = 0.0;

    if (hi == a.x[0]) {
        /* High word is integer already.  Round the low word. */
        lo = floor(a.x[1]);
        hi = quick_two_sum(hi, lo, &lo);
    }

    return dd_create(hi, lo);
}

static NPY_INLINE double2
dd_ceil(const double2 a)
{
    double hi = ceil(a.x[0]);
    double lo = 0.0;

    if (hi == a.x[0]) {
        /* High word is integer already.  Round the low word. */
        lo = ceil(a.x[1]);
        hi = quick_two_sum(hi, lo, &lo);
    }

    return dd_create(hi, lo);
}

static NPY_INLINE double2
dd_aint(const double2 a)
{
    return (a.x[0] >= 0.0) ? dd_floor(a) : dd_ceil(a);
}

/* Absolute value */
static NPY_INLINE double2
dd_abs(const double2 a)
{
    return (a.x[0] < 0.0 ? dd_neg(a) : a);
}

static NPY_INLINE double2
dd_fabs(const double2 a)
{
    return dd_abs(a);
}


/*********** Normalizing ***********/
/* double-double * (2.0 ^ expt) */
static NPY_INLINE double2
dd_ldexp(const double2 a, int expt)
{
    return dd_create(ldexp(a.x[0], expt), ldexp(a.x[1], expt));
}

static NPY_INLINE double2
dd_frexp(const double2 a, int *expt)
{
//    r"""return b and l s.t. 0.5<=|b|<1 and 2^l == a
//    0.5<=|b[0]|<1.0 or |b[0]| == 1.0 and b[0]*b[1]<0
//    """
    int exponent;
    double man = frexp(a.x[0], &exponent);
    double b1 = ldexp(a.x[1], -exponent);
    if (fabs(man) == 0.5 && man * b1 < 0)
    {
        man *=2;
        b1 *= 2;
        exponent -= 1;
    }
    *expt = exponent;
    return dd_create(man, b1);
}


/*********** Additions ************/
static NPY_INLINE double2
dd_add_d_d(double a, double b)
{
    double s, e;
    s = two_sum(a, b, &e);
    return dd_create(s, e);
}

static NPY_INLINE double2
dd_add_dd_d(const double2 a, double b)
{
    double s1, s2;
    s1 = two_sum(a.x[0], b, &s2);
    s2 += a.x[1];
    s1 = quick_two_sum(s1, s2, &s2);
    return dd_create(s1, s2);
}

static NPY_INLINE double2
dd_add_d_dd(double a, const double2 b)
{
    double s1, s2;
    s1 = two_sum(a, b.x[0], &s2);
    s2 += b.x[1];
    s1 = quick_two_sum(s1, s2, &s2);
    return dd_create(s1, s2);
}

static NPY_INLINE double2
dd_ieee_add(const double2 a, const double2 b)
{
    /* This one satisfies IEEE style error bound,
       due to K. Briggs and W. Kahan.                   */
    double s1, s2, t1, t2;

    s1 = two_sum(a.x[0], b.x[0], &s2);
    t1 = two_sum(a.x[1], b.x[1], &t2);
    s2 += t1;
    s1 = quick_two_sum(s1, s2, &s2);
    s2 += t2;
    s1 = quick_two_sum(s1, s2, &s2);
    return dd_create(s1, s2);
}

static NPY_INLINE double2
dd_sloppy_add(const double2 a, const double2 b)
{
    /* This is the less accurate version ... obeys Cray-style
       error bound. */
    double s, e;

    s = two_sum(a.x[0], b.x[0], &e);
    e += (a.x[1] + b.x[1]);
    s = quick_two_sum(s, e, &e);
    return dd_create(s, e);
}

static NPY_INLINE double2
dd_add(const double2 a, const double2 b)
{
    /* Always require IEEE-style error bounds */
    return dd_ieee_add(a, b);
}

/*********** Subtractions ************/
/* double-double = double - double */
static NPY_INLINE double2
dd_sub_d_d(double a, double b)
{
    double s, e;
    s = two_diff(a, b, &e);
    return dd_create(s, e);
}

static NPY_INLINE double2
dd_sub(const double2 a, const double2 b)
{
    return dd_ieee_add(a, dd_neg(b));
}

static NPY_INLINE double2
dd_sub_dd_d(const double2 a, double b)
{
    double s1, s2;
    s1 = two_sum(a.x[0], -b, &s2);
    s2 += a.x[1];
    s1 = quick_two_sum(s1, s2, &s2);
    return dd_create(s1, s2);
}

static NPY_INLINE double2
dd_sub_d_dd(double a, const double2 b)
{
    double s1, s2;
    s1 = two_sum(a, -b.x[0], &s2);
    s2 -= b.x[1];
    s1 = quick_two_sum(s1, s2, &s2);
    return dd_create(s1, s2);
}


/*********** Multiplications ************/
/* double-double = double * double */
static NPY_INLINE double2
dd_mul_d_d(double a, double b)
{
    double p, e;
    p = two_prod(a, b, &e);
    return dd_create(p, e);
}

/* double-double * double,  where double is a power of 2. */
static NPY_INLINE double2
dd_mul_pwr2(const double2 a, double b)
{
    return dd_create(a.x[0] * b, a.x[1] * b);
}

static NPY_INLINE double2
dd_mul(const double2 a, const double2 b)
{
    double p1, p2;
    p1 = two_prod(a.x[0], b.x[0], &p2);
    p2 += (a.x[0] * b.x[1] + a.x[1] * b.x[0]);
    p1 = quick_two_sum(p1, p2, &p2);
    return dd_create(p1, p2);
}

static NPY_INLINE double2
dd_mul_dd_d(const double2 a, double b)
{
    double p1, p2, e1, e2;
    p1 = two_prod(a.x[0], b, &e1);
    p2 = two_prod(a.x[1], b, &e2);
    p1 = quick_two_sum(p1, e2 + p2 + e1, &e1);
    return dd_create(p1, e1);
}

static NPY_INLINE double2
dd_mul_d_dd(double a, const double2 b)
{
    double p1, p2, e1, e2;
    p1 = two_prod(a, b.x[0], &e1);
    p2 = two_prod(a, b.x[1], &e2);
    p1 = quick_two_sum(p1, e2 + p2 + e1, &e1);
    return dd_create(p1, e1);
}


/*********** Divisions ************/
static NPY_INLINE double2
dd_sloppy_div(const double2 a, const double2 b)
{
    double s1, s2;
    double q1, q2;
    double2 r;

    q1 = a.x[0] / b.x[0]; /* approximate quotient */

    /* compute  this - q1 * dd */
    r = dd_sub(a, dd_mul_dd_d(b, q1));
    s1 = two_diff(a.x[0], r.x[0], &s2);
    s2 -= r.x[1];
    s2 += a.x[1];

    /* get next approximation */
    q2 = (s1 + s2) / b.x[0];

    /* renormalize */
    r.x[0] = quick_two_sum(q1, q2, &r.x[1]);
    return r;
}

static NPY_INLINE double2
dd_accurate_div(const double2 a, const double2 b)
{
    double q1, q2, q3;
    double2 r;

    q1 = a.x[0] / b.x[0]; /* approximate quotient */

    r = dd_sub(a, dd_mul_dd_d(b, q1));

    q2 = r.x[0] / b.x[0];
    r = dd_sub(r, dd_mul_dd_d(b, q2));

    q3 = r.x[0] / b.x[0];

    q1 = quick_two_sum(q1, q2, &q2);
    r = dd_add_dd_d(dd_create(q1, q2), q3);
    return r;
}

static NPY_INLINE double2
dd_div(const double2 a, const double2 b)
{
    return dd_accurate_div(a, b);
}

static NPY_INLINE double2
dd_div_d_d(double a, double b)
{
    return dd_accurate_div(dd_create_d(a), dd_create_d(b));
}

static NPY_INLINE double2
dd_div_dd_d(const double2 a, double b)
{
    return dd_accurate_div(a, dd_create_d(b));
}

static NPY_INLINE double2
dd_div_d_dd(double a, const double2 b)
{
    return dd_accurate_div(dd_create_d(a), b);
}

static NPY_INLINE double2
dd_inv(const double2 a)
{
    return dd_div(DD_C_ONE, a);
}


/********** Remainder **********/
static NPY_INLINE double2
dd_drem(const double2 a, const double2 b)
{
    double2 n = dd_nint(dd_div(a, b));
    return dd_sub(a, dd_mul(n, b));
}

static NPY_INLINE double2
dd_divrem(const double2 a, const double2 b, double2 *r)
{
    double2 n = dd_nint(dd_div(a, b));
    *r = dd_sub(a, dd_mul(n, b));
    return n;
}

static NPY_INLINE double2
dd_fmod(const double2 a, const double2 b)
{
    double2 n = dd_aint(dd_div(a, b));
    return dd_sub(a, dd_mul(b, n));
}

/*********** Squaring **********/
static NPY_INLINE double2
dd_sqr(const double2 a)
{
    double p1, p2;
    double s1, s2;
    p1 = two_sqr(a.x[0], &p2);
    p2 += 2.0 * a.x[0] * a.x[1];
    p2 += a.x[1] * a.x[1];
    s1 = quick_two_sum(p1, p2, &s2);
    return dd_create(s1, s2);
}

static NPY_INLINE double2
dd_sqr_d(double a)
{
    double p1, p2;
    p1 = two_sqr(a, &p2);
    return dd_create(p1, p2);
}

#ifdef __cplusplus
}
#endif

#endif  /*  _DD_REAL_IDEFS_H_ */
