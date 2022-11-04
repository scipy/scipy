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

#ifndef  _DD_IDEFS_H_
#define  _DD_IDEFS_H_ 1

#include <float.h>
#include <limits.h>
#include <math.h>

#include <numpy/utils.h>

#ifdef __cplusplus
extern "C" {
#endif

#define _DD_SPLITTER 134217729.0               // = 2^27 + 1
#define _DD_SPLIT_THRESH 6.69692879491417e+299 // = 2^996

/*
 ************************************************************************
  The basic routines taking double arguments, returning 1 (or 2) doubles
 ************************************************************************
*/

/* Computes fl(a+b) and err(a+b).  Assumes |a| >= |b|. */
static inline double
quick_two_sum(double a, double b, double *err)
{
    volatile double s = a + b;
    volatile double c = s - a;
    *err = b - c;
    return s;
}

/* Computes fl(a-b) and err(a-b).  Assumes |a| >= |b| */
static inline double
quick_two_diff(double a, double b, double *err)
{
    volatile double s = a - b;
    volatile double c = a - s;
    *err = c - b;
    return s;
}

/* Computes fl(a+b) and err(a+b).  */
static inline double
two_sum(double a, double b, double *err)
{
    volatile double s = a + b;
    volatile double c = s - a;
    volatile double d = b - c;
    volatile double e = s - c;
    *err = (a - e) + d;
    return s;
}

/* Computes fl(a-b) and err(a-b).  */
static inline double
two_diff(double a, double b, double *err)
{
    volatile double s = a - b;
    volatile double c = s - a;
    volatile double d = b + c;
    volatile double e = s - c;
    *err = (a - e) - d;
    return s;
}

/* Computes high word and lo word of a */
static inline void
two_split(double a, double *hi, double *lo)
{
    volatile double temp, tempma;
    if (a > _DD_SPLIT_THRESH || a < -_DD_SPLIT_THRESH) {
        a *= 3.7252902984619140625e-09; // 2^-28
        temp = _DD_SPLITTER * a;
        tempma = temp - a;
        *hi = temp - tempma;
        *lo = a - *hi;
        *hi *= 268435456.0; // 2^28
        *lo *= 268435456.0; // 2^28
    }
    else {
        temp = _DD_SPLITTER * a;
        tempma = temp - a;
        *hi = temp - tempma;
        *lo = a - *hi;
    }
}

/* Computes fl(a*b) and err(a*b). */
static inline double
two_prod(double a, double b, double *err)
{
#ifdef DD_FMS
    volatile double p = a * b;
    *err = DD_FMS(a, b, p);
    return p;
#else
    double a_hi, a_lo, b_hi, b_lo;
    double p = a * b;
    volatile double c, d;
    two_split(a, &a_hi, &a_lo);
    two_split(b, &b_hi, &b_lo);
    c = a_hi * b_hi - p;
    d = c + a_hi * b_lo + a_lo * b_hi;
    *err = d + a_lo * b_lo;
    return p;
#endif  /* DD_FMA */
}

/* Computes fl(a*a) and err(a*a).  Faster than the above method. */
static inline double
two_sqr(double a, double *err)
{
#ifdef DD_FMS
    volatile double p = a * a;
    *err = DD_FMS(a, a, p);
    return p;
#else
    double hi, lo;
    volatile double c;
    double q = a * a;
    two_split(a, &hi, &lo);
    c = hi * hi - q;
    *err = (c + 2.0 * hi * lo) + lo * lo;
    return q;
#endif /* DD_FMS */
}

static inline double
two_div(double a, double b, double *err)
{
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
static inline double
two_nint(double d)
{
    if (d == floor(d)) {
        return d;
    }
    return floor(d + 0.5);
}

/* Computes the truncated integer. */
static inline double
two_aint(double d)
{
    return (d >= 0.0 ? floor(d) : ceil(d));
}


/* Compare a and b */
static inline int
two_comp(const double a, const double b)
{
    /* Works for non-NAN inputs */
    return (a < b ? -1 : (a > b ? 1 : 0));
}


#ifdef __cplusplus
}
#endif

#endif  /*  _DD_IDEFS_H_ */
