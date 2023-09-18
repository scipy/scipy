/* Copyright (C) 1997, 1998, 1999, 2000, 2001, 2003, 2004, 2005, 2007
   Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

/*
 *	ISO C99 Standard: 7.22 Type-generic math	<tgmath.h>
 */

#ifndef _TGMATH_H
#define _TGMATH_H	1

/* Include the needed headers.  */
#include <math.h>
#include <complex.h>


/* Since `complex' is currently not really implemented in most C compilers
   and if it is implemented, the implementations differ.  This makes it
   quite difficult to write a generic implementation of this header.  We
   do not try this for now and instead concentrate only on GNU CC.  Once
   we have more information support for other compilers might follow.  */

#if __GNUC_PREREQ (2, 7)

# ifdef __NO_LONG_DOUBLE_MATH
#  define __tgml(fct) fct
# else
#  define __tgml(fct) fct ## l
# endif

/* This is ugly but unless gcc gets appropriate builtins we have to do
   something like this.  Don't ask how it works.  */

/* 1 if 'type' is a floating type, 0 if 'type' is an integer type.
   Allows for _Bool.  Expands to an integer constant expression.  */
# if __GNUC_PREREQ (3, 1)
#  define __floating_type(type) \
  (__builtin_classify_type ((type) 0) == 8 \
   || (__builtin_classify_type ((type) 0) == 9 \
       && __builtin_classify_type (__real__ ((type) 0)) == 8))
# else
#  define __floating_type(type) (((type) 0.25) && ((type) 0.25 - 1))
# endif

/* The tgmath real type for T, where E is 0 if T is an integer type and
   1 for a floating type.  */
# define __tgmath_real_type_sub(T, E) \
  __typeof__ (*(0 ? (__typeof__ (0 ? (double *) 0 : (void *) (E))) 0	      \
		  : (__typeof__ (0 ? (T *) 0 : (void *) (!(E)))) 0))

/* The tgmath real type of EXPR.  */
# define __tgmath_real_type(expr) \
  __tgmath_real_type_sub (__typeof__ ((__typeof__ (expr)) 0),		      \
			  __floating_type (__typeof__ (expr)))


/* We have two kinds of generic macros: to support functions which are
   only defined on real valued parameters and those which are defined
   for complex functions as well.  */
# define __TGMATH_UNARY_REAL_ONLY(Val, Fct) \
     (__extension__ ((sizeof (Val) == sizeof (double)			      \
		      || __builtin_classify_type (Val) != 8)		      \
		     ? (__tgmath_real_type (Val)) Fct (Val)		      \
		     : (sizeof (Val) == sizeof (float))			      \
		     ? (__tgmath_real_type (Val)) Fct##f (Val)		      \
		     : (__tgmath_real_type (Val)) __tgml(Fct) (Val)))

# define __TGMATH_UNARY_REAL_RET_ONLY(Val, RetType, Fct) \
     (__extension__ ((sizeof (Val) == sizeof (double)			      \
		      || __builtin_classify_type (Val) != 8)		      \
		     ? (RetType) Fct (Val)				      \
		     : (sizeof (Val) == sizeof (float))			      \
		     ? (RetType) Fct##f (Val)				      \
		     : (RetType) __tgml(Fct) (Val)))

# define __TGMATH_BINARY_FIRST_REAL_ONLY(Val1, Val2, Fct) \
     (__extension__ ((sizeof (Val1) == sizeof (double)			      \
		      || __builtin_classify_type (Val1) != 8)		      \
		     ? (__tgmath_real_type (Val1)) Fct (Val1, Val2)	      \
		     : (sizeof (Val1) == sizeof (float))		      \
		     ? (__tgmath_real_type (Val1)) Fct##f (Val1, Val2)	      \
		     : (__tgmath_real_type (Val1)) __tgml(Fct) (Val1, Val2)))

# define __TGMATH_BINARY_REAL_ONLY(Val1, Val2, Fct) \
     (__extension__ (((sizeof (Val1) > sizeof (double)			      \
		       || sizeof (Val2) > sizeof (double))		      \
		      && __builtin_classify_type ((Val1) + (Val2)) == 8)      \
		     ? (__typeof ((__tgmath_real_type (Val1)) 0		      \
				   + (__tgmath_real_type (Val2)) 0))	      \
		       __tgml(Fct) (Val1, Val2)				      \
		     : (sizeof (Val1) == sizeof (double)		      \
			|| sizeof (Val2) == sizeof (double)		      \
			|| __builtin_classify_type (Val1) != 8		      \
			|| __builtin_classify_type (Val2) != 8)		      \
		     ? (__typeof ((__tgmath_real_type (Val1)) 0		      \
				   + (__tgmath_real_type (Val2)) 0))	      \
		       Fct (Val1, Val2)					      \
		     : (__typeof ((__tgmath_real_type (Val1)) 0		      \
				   + (__tgmath_real_type (Val2)) 0))	      \
		       Fct##f (Val1, Val2)))

# define __TGMATH_TERNARY_FIRST_SECOND_REAL_ONLY(Val1, Val2, Val3, Fct) \
     (__extension__ (((sizeof (Val1) > sizeof (double)			      \
		       || sizeof (Val2) > sizeof (double))		      \
		      && __builtin_classify_type ((Val1) + (Val2)) == 8)      \
		     ? (__typeof ((__tgmath_real_type (Val1)) 0		      \
				   + (__tgmath_real_type (Val2)) 0))	      \
		       __tgml(Fct) (Val1, Val2, Val3)			      \
		     : (sizeof (Val1) == sizeof (double)		      \
			|| sizeof (Val2) == sizeof (double)		      \
			|| __builtin_classify_type (Val1) != 8		      \
			|| __builtin_classify_type (Val2) != 8)		      \
		     ? (__typeof ((__tgmath_real_type (Val1)) 0		      \
				   + (__tgmath_real_type (Val2)) 0))	      \
		       Fct (Val1, Val2, Val3)				      \
		     : (__typeof ((__tgmath_real_type (Val1)) 0		      \
				   + (__tgmath_real_type (Val2)) 0))	      \
		       Fct##f (Val1, Val2, Val3)))

# define __TGMATH_TERNARY_REAL_ONLY(Val1, Val2, Val3, Fct) \
     (__extension__ (((sizeof (Val1) > sizeof (double)			      \
		       || sizeof (Val2) > sizeof (double)		      \
		       || sizeof (Val3) > sizeof (double))		      \
		      && __builtin_classify_type ((Val1) + (Val2) + (Val3))   \
			 == 8)						      \
		     ? (__typeof ((__tgmath_real_type (Val1)) 0		      \
				   + (__tgmath_real_type (Val2)) 0	      \
				   + (__tgmath_real_type (Val3)) 0))	      \
		       __tgml(Fct) (Val1, Val2, Val3)			      \
		     : (sizeof (Val1) == sizeof (double)		      \
			|| sizeof (Val2) == sizeof (double)		      \
			|| sizeof (Val3) == sizeof (double)		      \
			|| __builtin_classify_type (Val1) != 8		      \
			|| __builtin_classify_type (Val2) != 8		      \
			|| __builtin_classify_type (Val3) != 8)		      \
		     ? (__typeof ((__tgmath_real_type (Val1)) 0		      \
				   + (__tgmath_real_type (Val2)) 0	      \
				   + (__tgmath_real_type (Val3)) 0))	      \
		       Fct (Val1, Val2, Val3)				      \
		     : (__typeof ((__tgmath_real_type (Val1)) 0		      \
				   + (__tgmath_real_type (Val2)) 0	      \
				   + (__tgmath_real_type (Val3)) 0))	      \
		       Fct##f (Val1, Val2, Val3)))

/* XXX This definition has to be changed as soon as the compiler understands
   the imaginary keyword.  */
# define __TGMATH_UNARY_REAL_IMAG(Val, Fct, Cfct) \
     (__extension__ ((sizeof (__real__ (Val)) == sizeof (double)	      \
		      || __builtin_classify_type (__real__ (Val)) != 8)	      \
		     ? ((sizeof (__real__ (Val)) == sizeof (Val))	      \
			? (__tgmath_real_type (Val)) Fct (Val)		      \
			: (__tgmath_real_type (Val)) Cfct (Val))	      \
		     : (sizeof (__real__ (Val)) == sizeof (float))	      \
		     ? ((sizeof (__real__ (Val)) == sizeof (Val))	      \
			? (__tgmath_real_type (Val)) Fct##f (Val)	      \
			: (__tgmath_real_type (Val)) Cfct##f (Val))	      \
		     : ((sizeof (__real__ (Val)) == sizeof (Val))	      \
			? (__tgmath_real_type (Val)) __tgml(Fct) (Val)	      \
			: (__tgmath_real_type (Val)) __tgml(Cfct) (Val))))

# define __TGMATH_UNARY_IMAG(Val, Cfct) \
     (__extension__ ((sizeof (__real__ (Val)) == sizeof (double)	      \
		      || __builtin_classify_type (__real__ (Val)) != 8)	      \
		     ? (__typeof__ ((__tgmath_real_type (Val)) 0	      \
				    + _Complex_I)) Cfct (Val)		      \
		     : (sizeof (__real__ (Val)) == sizeof (float))	      \
		     ? (__typeof__ ((__tgmath_real_type (Val)) 0	      \
				    + _Complex_I)) Cfct##f (Val)	      \
		     : (__typeof__ ((__tgmath_real_type (Val)) 0	      \
				    + _Complex_I)) __tgml(Cfct) (Val)))

/* XXX This definition has to be changed as soon as the compiler understands
   the imaginary keyword.  */
# define __TGMATH_UNARY_REAL_IMAG_RET_REAL(Val, Fct, Cfct) \
     (__extension__ ((sizeof (__real__ (Val)) == sizeof (double)	      \
		      || __builtin_classify_type (__real__ (Val)) != 8)	      \
		     ? ((sizeof (__real__ (Val)) == sizeof (Val))	      \
			? (__typeof__ (__real__ (__tgmath_real_type (Val)) 0))\
			  Fct (Val)					      \
			: (__typeof__ (__real__ (__tgmath_real_type (Val)) 0))\
			  Cfct (Val))					      \
		     : (sizeof (__real__ (Val)) == sizeof (float))	      \
		     ? ((sizeof (__real__ (Val)) == sizeof (Val))	      \
			? (__typeof__ (__real__ (__tgmath_real_type (Val)) 0))\
			  Fct##f (Val)					      \
			: (__typeof__ (__real__ (__tgmath_real_type (Val)) 0))\
			  Cfct##f (Val))				      \
		     : ((sizeof (__real__ (Val)) == sizeof (Val))	      \
			? (__typeof__ (__real__ (__tgmath_real_type (Val)) 0))\
			  __tgml(Fct) (Val)				      \
			: (__typeof__ (__real__ (__tgmath_real_type (Val)) 0))\
			  __tgml(Cfct) (Val))))

/* XXX This definition has to be changed as soon as the compiler understands
   the imaginary keyword.  */
# define __TGMATH_BINARY_REAL_IMAG(Val1, Val2, Fct, Cfct) \
     (__extension__ (((sizeof (__real__ (Val1)) > sizeof (double)	      \
		       || sizeof (__real__ (Val2)) > sizeof (double))	      \
		      && __builtin_classify_type (__real__ (Val1)	      \
						  + __real__ (Val2)) == 8)    \
		     ? ((sizeof (__real__ (Val1)) == sizeof (Val1)	      \
			 && sizeof (__real__ (Val2)) == sizeof (Val2))	      \
			? (__typeof ((__tgmath_real_type (Val1)) 0	      \
				   + (__tgmath_real_type (Val2)) 0))	      \
			  __tgml(Fct) (Val1, Val2)			      \
			: (__typeof ((__tgmath_real_type (Val1)) 0	      \
				   + (__tgmath_real_type (Val2)) 0))	      \
			  __tgml(Cfct) (Val1, Val2))			      \
		     : (sizeof (__real__ (Val1)) == sizeof (double)	      \
			|| sizeof (__real__ (Val2)) == sizeof (double)	      \
			|| __builtin_classify_type (__real__ (Val1)) != 8     \
			|| __builtin_classify_type (__real__ (Val2)) != 8)    \
		     ? ((sizeof (__real__ (Val1)) == sizeof (Val1)	      \
			 && sizeof (__real__ (Val2)) == sizeof (Val2))	      \
			? (__typeof ((__tgmath_real_type (Val1)) 0	      \
				   + (__tgmath_real_type (Val2)) 0))	      \
			  Fct (Val1, Val2)				      \
			: (__typeof ((__tgmath_real_type (Val1)) 0	      \
				   + (__tgmath_real_type (Val2)) 0))	      \
			  Cfct (Val1, Val2))				      \
		     : ((sizeof (__real__ (Val1)) == sizeof (Val1)	      \
			 && sizeof (__real__ (Val2)) == sizeof (Val2))	      \
			? (__typeof ((__tgmath_real_type (Val1)) 0	      \
				   + (__tgmath_real_type (Val2)) 0))	      \
			  Fct##f (Val1, Val2)				      \
			: (__typeof ((__tgmath_real_type (Val1)) 0	      \
				   + (__tgmath_real_type (Val2)) 0))	      \
			  Cfct##f (Val1, Val2))))
#else
# error "Unsupported compiler; you cannot use <tgmath.h>"
#endif


/* Unary functions defined for real and complex values.  */


/* Trigonometric functions.  */

/* Arc cosine of X.  */
#define acos(Val) __TGMATH_UNARY_REAL_IMAG (Val, acos, cacos)
/* Arc sine of X.  */
#define asin(Val) __TGMATH_UNARY_REAL_IMAG (Val, asin, casin)
/* Arc tangent of X.  */
#define atan(Val) __TGMATH_UNARY_REAL_IMAG (Val, atan, catan)
/* Arc tangent of Y/X.  */
#define atan2(Val1, Val2) __TGMATH_BINARY_REAL_ONLY (Val1, Val2, atan2)

/* Cosine of X.  */
#define cos(Val) __TGMATH_UNARY_REAL_IMAG (Val, cos, ccos)
/* Sine of X.  */
#define sin(Val) __TGMATH_UNARY_REAL_IMAG (Val, sin, csin)
/* Tangent of X.  */
#define tan(Val) __TGMATH_UNARY_REAL_IMAG (Val, tan, ctan)


/* Hyperbolic functions.  */

/* Hyperbolic arc cosine of X.  */
#define acosh(Val) __TGMATH_UNARY_REAL_IMAG (Val, acosh, cacosh)
/* Hyperbolic arc sine of X.  */
#define asinh(Val) __TGMATH_UNARY_REAL_IMAG (Val, asinh, casinh)
/* Hyperbolic arc tangent of X.  */
#define atanh(Val) __TGMATH_UNARY_REAL_IMAG (Val, atanh, catanh)

/* Hyperbolic cosine of X.  */
#define cosh(Val) __TGMATH_UNARY_REAL_IMAG (Val, cosh, ccosh)
/* Hyperbolic sine of X.  */
#define sinh(Val) __TGMATH_UNARY_REAL_IMAG (Val, sinh, csinh)
/* Hyperbolic tangent of X.  */
#define tanh(Val) __TGMATH_UNARY_REAL_IMAG (Val, tanh, ctanh)


/* Exponential and logarithmic functions.  */

/* Exponential function of X.  */
#define exp(Val) __TGMATH_UNARY_REAL_IMAG (Val, exp, cexp)

/* Break VALUE into a normalized fraction and an integral power of 2.  */
#define frexp(Val1, Val2) __TGMATH_BINARY_FIRST_REAL_ONLY (Val1, Val2, frexp)

/* X times (two to the EXP power).  */
#define ldexp(Val1, Val2) __TGMATH_BINARY_FIRST_REAL_ONLY (Val1, Val2, ldexp)

/* Natural logarithm of X.  */
#define log(Val) __TGMATH_UNARY_REAL_IMAG (Val, log, clog)

/* Base-ten logarithm of X.  */
#ifdef __USE_GNU
# define log10(Val) __TGMATH_UNARY_REAL_IMAG (Val, log10, __clog10)
#else
# define log10(Val) __TGMATH_UNARY_REAL_ONLY (Val, log10)
#endif

/* Return exp(X) - 1.  */
#define expm1(Val) __TGMATH_UNARY_REAL_ONLY (Val, expm1)

/* Return log(1 + X).  */
#define log1p(Val) __TGMATH_UNARY_REAL_ONLY (Val, log1p)

/* Return the base 2 signed integral exponent of X.  */
#define logb(Val) __TGMATH_UNARY_REAL_ONLY (Val, logb)

/* Compute base-2 exponential of X.  */
#define exp2(Val) __TGMATH_UNARY_REAL_ONLY (Val, exp2)

/* Compute base-2 logarithm of X.  */
#define log2(Val) __TGMATH_UNARY_REAL_ONLY (Val, log2)


/* Power functions.  */

/* Return X to the Y power.  */
#define pow(Val1, Val2) __TGMATH_BINARY_REAL_IMAG (Val1, Val2, pow, cpow)

/* Return the square root of X.  */
#define sqrt(Val) __TGMATH_UNARY_REAL_IMAG (Val, sqrt, csqrt)

/* Return `sqrt(X*X + Y*Y)'.  */
#define hypot(Val1, Val2) __TGMATH_BINARY_REAL_ONLY (Val1, Val2, hypot)

/* Return the cube root of X.  */
#define cbrt(Val) __TGMATH_UNARY_REAL_ONLY (Val, cbrt)


/* Nearest integer, absolute value, and remainder functions.  */

/* Smallest integral value not less than X.  */
#define ceil(Val) __TGMATH_UNARY_REAL_ONLY (Val, ceil)

/* Absolute value of X.  */
#define fabs(Val) __TGMATH_UNARY_REAL_IMAG_RET_REAL (Val, fabs, cabs)

/* Largest integer not greater than X.  */
#define floor(Val) __TGMATH_UNARY_REAL_ONLY (Val, floor)

/* Floating-point modulo remainder of X/Y.  */
#define fmod(Val1, Val2) __TGMATH_BINARY_REAL_ONLY (Val1, Val2, fmod)

/* Round X to integral valuein floating-point format using current
   rounding direction, but do not raise inexact exception.  */
#define nearbyint(Val) __TGMATH_UNARY_REAL_ONLY (Val, nearbyint)

/* Round X to nearest integral value, rounding halfway cases away from
   zero.  */
#define round(Val) __TGMATH_UNARY_REAL_ONLY (Val, round)

/* Round X to the integral value in floating-point format nearest but
   not larger in magnitude.  */
#define trunc(Val) __TGMATH_UNARY_REAL_ONLY (Val, trunc)

/* Compute remainder of X and Y and put in *QUO a value with sign of x/y
   and magnitude congruent `mod 2^n' to the magnitude of the integral
   quotient x/y, with n >= 3.  */
#define remquo(Val1, Val2, Val3) \
     __TGMATH_TERNARY_FIRST_SECOND_REAL_ONLY (Val1, Val2, Val3, remquo)

/* Round X to nearest integral value according to current rounding
   direction.  */
#define lrint(Val) __TGMATH_UNARY_REAL_RET_ONLY (Val, long int, lrint)
#define llrint(Val) __TGMATH_UNARY_REAL_RET_ONLY (Val, long long int, llrint)

/* Round X to nearest integral value, rounding halfway cases away from
   zero.  */
#define lround(Val) __TGMATH_UNARY_REAL_RET_ONLY (Val, long int, lround)
#define llround(Val) __TGMATH_UNARY_REAL_RET_ONLY (Val, long long int, llround)


/* Return X with its signed changed to Y's.  */
#define copysign(Val1, Val2) __TGMATH_BINARY_REAL_ONLY (Val1, Val2, copysign)

/* Error and gamma functions.  */
#define erf(Val) __TGMATH_UNARY_REAL_ONLY (Val, erf)
#define erfc(Val) __TGMATH_UNARY_REAL_ONLY (Val, erfc)
#define tgamma(Val) __TGMATH_UNARY_REAL_ONLY (Val, tgamma)
#define lgamma(Val) __TGMATH_UNARY_REAL_ONLY (Val, lgamma)


/* Return the integer nearest X in the direction of the
   prevailing rounding mode.  */
#define rint(Val) __TGMATH_UNARY_REAL_ONLY (Val, rint)

/* Return X + epsilon if X < Y, X - epsilon if X > Y.  */
#define nextafter(Val1, Val2) __TGMATH_BINARY_REAL_ONLY (Val1, Val2, nextafter)
#define nexttoward(Val1, Val2) \
     __TGMATH_BINARY_FIRST_REAL_ONLY (Val1, Val2, nexttoward)

/* Return the remainder of integer divison X / Y with infinite precision.  */
#define remainder(Val1, Val2) __TGMATH_BINARY_REAL_ONLY (Val1, Val2, remainder)

/* Return X times (2 to the Nth power).  */
#if defined __USE_MISC || defined __USE_XOPEN_EXTENDED
# define scalb(Val1, Val2) __TGMATH_BINARY_REAL_ONLY (Val1, Val2, scalb)
#endif

/* Return X times (2 to the Nth power).  */
#define scalbn(Val1, Val2) __TGMATH_BINARY_FIRST_REAL_ONLY (Val1, Val2, scalbn)

/* Return X times (2 to the Nth power).  */
#define scalbln(Val1, Val2) \
     __TGMATH_BINARY_FIRST_REAL_ONLY (Val1, Val2, scalbln)

/* Return the binary exponent of X, which must be nonzero.  */
#define ilogb(Val) __TGMATH_UNARY_REAL_RET_ONLY (Val, int, ilogb)


/* Return positive difference between X and Y.  */
#define fdim(Val1, Val2) __TGMATH_BINARY_REAL_ONLY (Val1, Val2, fdim)

/* Return maximum numeric value from X and Y.  */
#define fmax(Val1, Val2) __TGMATH_BINARY_REAL_ONLY (Val1, Val2, fmax)

/* Return minimum numeric value from X and Y.  */
#define fmin(Val1, Val2) __TGMATH_BINARY_REAL_ONLY (Val1, Val2, fmin)


/* Multiply-add function computed as a ternary operation.  */
#define fma(Val1, Val2, Val3) \
     __TGMATH_TERNARY_REAL_ONLY (Val1, Val2, Val3, fma)


/* Absolute value, conjugates, and projection.  */

/* Argument value of Z.  */
#define carg(Val) __TGMATH_UNARY_REAL_IMAG_RET_REAL (Val, carg, carg)

/* Complex conjugate of Z.  */
#define conj(Val) __TGMATH_UNARY_IMAG (Val, conj)

/* Projection of Z onto the Riemann sphere.  */
#define cproj(Val) __TGMATH_UNARY_IMAG (Val, cproj)


/* Decomposing complex values.  */

/* Imaginary part of Z.  */
#define cimag(Val) __TGMATH_UNARY_REAL_IMAG_RET_REAL (Val, cimag, cimag)

/* Real part of Z.  */
#define creal(Val) __TGMATH_UNARY_REAL_IMAG_RET_REAL (Val, creal, creal)

#endif /* tgmath.h */
