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

#ifndef _DD_REAL_H
#define _DD_REAL_H

#include <float.h>
#include <limits.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Some configuration defines */

/* If fast fused multiply-add is available, define to the correct macro for
   using it.  It is invoked as DD_FMA(a, b, c) to compute fl(a * b + c).
   If correctly rounded multiply-add is not available (or if unsure),
   keep it undefined. */
#ifndef DD_FMA
#ifdef FP_FAST_FMA
#define DD_FMA(A, B, C) fma((A), (B), (C))
#endif
#endif

/* Same with fused multiply-subtract */
#ifndef DD_FMS
#ifdef FP_FAST_FMA
#define DD_FMS(A, B, C) fma((A), (B), (-C))
#endif
#endif

#ifdef __cplusplus
#define DD_STATIC_CAST(T, X) (static_cast<T>(X))
#else
#define DD_STATIC_CAST(T, X) ((T)(X))
#endif

/* double2 struct definition, some external always-present double2 constants.
*/
typedef struct double2
{
    double x[2];
} double2;

extern const double DD_C_EPS;
extern const double DD_C_MIN_NORMALIZED;
extern const double2 DD_C_MAX;
extern const double2 DD_C_SAFE_MAX;
extern const int DD_C_NDIGITS;

extern const double2 DD_C_2PI;
extern const double2 DD_C_PI;
extern const double2 DD_C_3PI4;
extern const double2 DD_C_PI2;
extern const double2 DD_C_PI4;
extern const double2 DD_C_PI16;
extern const double2 DD_C_E;
extern const double2 DD_C_LOG2;
extern const double2 DD_C_LOG10;
extern const double2 DD_C_ZERO;
extern const double2 DD_C_ONE;
extern const double2 DD_C_NEGONE;

/* NAN definition in AIX's math.h doesn't make it qualify as constant literal. */
#if defined(__STDC__) && defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L) && defined(NAN) && !defined(_AIX)
#define DD_C_NAN_IS_CONST
extern const double2 DD_C_NAN;
extern const double2 DD_C_INF;
extern const double2 DD_C_NEGINF;
#else
#define DD_C_NAN (dd_create(NAN, NAN))
#define DD_C_INF (dd_create(INFINITY, INFINITY))
#define DD_C_NEGINF (dd_create(-INFINITY, -INFINITY))
#endif


/* Include the inline definitions of functions */
#include "dd_real_idefs.h"

/* Non-inline functions */

/********** Exponentiation **********/
double2 dd_npwr(const double2 a, int n);

/*********** Transcendental Functions ************/
double2 dd_exp(const double2 a);
double2 dd_log(const double2 a);
double2 dd_expm1(const double2 a);
double2 dd_log1p(const double2 a);
double2 dd_log10(const double2 a);
double2 dd_log_d(double a);

/* Returns the exponent of the double precision number.
   Returns INT_MIN is x is zero, and INT_MAX if x is INF or NaN. */
int get_double_expn(double x);

/*********** Polynomial Functions ************/
double2 dd_polyeval(const double2 *c, int n, const double2 x);

/*********** Random number generator ************/
extern double2 dd_rand(void);


#ifdef __cplusplus
}
#endif


#endif /* _DD_REAL_H */
