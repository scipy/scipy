/*
 * include/double2.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 *
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
#include <stdbool.h>

#include "_c99compat.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Some configuration defines */

/* If fast fused multiply-add is available, define to the correct macro for
   using it.  It is invoked as DD_FMA(a, b, c) to compute fl(a * b + c).
   If correctly rounded multiply-add is not available (or if unsure),
   keep it undefined.*/
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

/* Define this macro to be the isfinite(x) function. */
#define DD_ISFINITE isfinite

/* Define this macro to be the isinf(x) function. */
#define DD_ISINF isinf

/* Define this macro to be the isnan(x) function. */
#define DD_ISNAN isnan


/* Set the following to 1 to set commonly used function
   to be inlined.  This should be set to 1 unless the compiler
   does not support the "inline" keyword, or if building for
   debugging purposes. */
#define DD_DO_INLINE 1

#ifdef DD_DO_INLINE
#define DD_INLINE inline
#else
#define DD_INLINE
#endif

typedef struct double2
{
    double x[2];
} double2;

extern const double2 DD_C_2PI;
extern const double2 DD_C_PI;
extern const double2 DD_C_3PI4;
extern const double2 DD_C_PI2;
extern const double2 DD_C_PI4;
extern const double2 DD_C_PI16;
extern const double2 DD_C_E;
extern const double2 DD_C_LOG2;
extern const double2 DD_C_LOG10;
extern const double2 DD_C_NAN;
extern const double2 DD_C_INF;
extern const double2 DD_C_NEGINF;
extern const double2 DD_C_ZERO;
extern const double2 DD_C_ONE;

extern const double DD_C_EPS;
extern const double DD_C_MIN_NORMALIZED;
extern const double2 DD_C_MAX;
extern const double2 DD_C_SAFE_MAX;
extern const int DD_C_NDIGITS;

#ifdef __cplusplus
#define DD_STATIC_CAST(T, X) (static_cast<T>(X))
#else
#define DD_STATIC_CAST(T, X) ((T)(X))
#endif

#ifdef DD_DO_INLINE
#include "dd_real_inline.h"

#else

/*********** Inline basic functions ************/

DD_INLINE double dd_hi(const double2 a);
DD_INLINE double dd_lo(const double2 a);
DD_INLINE bool dd_isnan(const double2 a);
DD_INLINE bool dd_isfinite(const double2 a);
DD_INLINE bool dd_isinf(const double2 a);
DD_INLINE bool dd_is_zero(const double2 a);
DD_INLINE bool dd_is_one(const double2 a);
DD_INLINE bool dd_is_positive(const double2 a);
DD_INLINE bool dd_is_negative(const double2 a);
DD_INLINE double dd_to_double(const double2 a);
DD_INLINE int dd_to_int(const double2 a);

/*********** Equality Comparisons ************/

DD_INLINE int dd_comp(const double2 a, const double2 b);
DD_INLINE int dd_comp_dd_d(const double2 a, double b);
DD_INLINE int dd_comp_d_dd(double a, const double2 b);

/*********** Creation ************/
DD_INLINE double2 dd_create(double hi, double lo);
DD_INLINE double2 dd_zero(void);
DD_INLINE double2 dd_create_d(double h);
DD_INLINE double2 dd_create_i(int h);
DD_INLINE double2 dd_create_dp(const double *d);

/*********** Unary Minus ***********/
DD_INLINE double2 dd_neg(double2 a);

/*********** Rounding to integer ***********/

DD_INLINE double2 dd_nint(const double2 a);
DD_INLINE double2 dd_floor(const double2 a);
DD_INLINE double2 dd_ceil(const double2 a);
DD_INLINE double2 dd_aint(const double2 a);

DD_INLINE double2 dd_abs(const double2 a);
DD_INLINE double2 dd_fabs(const double2 a);

/*********** Normalizing ***********/

DD_INLINE double2 dd_ldexp(const double2 a, int expt);
DD_INLINE double2 dd_frexp(const double2 a, int *expt);

/*********** Additions/Subtractions ************/
DD_INLINE double2 dd_add_d_d(double a, double b);
DD_INLINE double2 dd_ieee_add(const double2 a, const double2 b);
DD_INLINE double2 dd_sloppy_add(const double2 a, const double2 b);
DD_INLINE double2 dd_add(const double2 a, const double2 b);
DD_INLINE double2 dd_add_dd_d(const double2 a, double b);
DD_INLINE double2 dd_add_d_dd(double a, const double2 b);

DD_INLINE double2 dd_sub_d_d(double a, double b);
DD_INLINE double2 dd_sub(const double2 a, const double2 b);
DD_INLINE double2 dd_sub_dd_d(const double2 a, double b);
DD_INLINE double2 dd_sub_d_dd(double a, const double2 b);

/*********** Multiplications ************/
DD_INLINE double2 dd_mul_d_d(double a, double b);
DD_INLINE double2 dd_mul_pwr2(const double2 a, double b);
DD_INLINE double2 dd_mul_dd_d(const double2 a, double b);
DD_INLINE double2 dd_mul_d_dd(double a, const double2 b);
DD_INLINE double2 dd_mul(const double2 a, const double2 b);

/*********** Divisions ************/

DD_INLINE double2 dd_sloppy_div(const double2 a, const double2 b);
DD_INLINE double2 dd_accurate_div(const double2 a, const double2 b);
DD_INLINE double2 dd_div(const double2 a, const double2 b);
DD_INLINE double2 dd_inv(const double2 a);

DD_INLINE double2 dd_div_d_dd(double a, const double2 b);
DD_INLINE double2 dd_div_dd_d(double2 a, double b);
DD_INLINE double2 dd_div_d_d(double a, double b);

/********** Remainder **********/
DD_INLINE double2 dd_drem(const double2 a, const double2 b);
DD_INLINE double2 dd_divrem(const double2 a, const double2 b, double2 *r);
DD_INLINE double2 dd_fmod(const double2 a, const double2 b);

/*********** Squaring **********/
DD_INLINE double2 dd_sqr(const double2 a);
DD_INLINE double2 dd_sqr_d(double a);

/*********** Random number generator ************/
DD_INLINE double2 dd_rand(void);

#endif  /* DD_INLINE */


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


#ifdef __cplusplus
}
#endif


#endif /* _DD_REAL_H */
