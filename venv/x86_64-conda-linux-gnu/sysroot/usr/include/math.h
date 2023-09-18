/* Declarations for math functions.
   Copyright (C) 1991-1993, 1995-1999, 2001, 2002, 2004, 2006, 2009
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
 *	ISO C99 Standard: 7.12 Mathematics	<math.h>
 */

#ifndef	_MATH_H
#define	_MATH_H	1

#include <features.h>

__BEGIN_DECLS

/* Get machine-dependent HUGE_VAL value (returned on overflow).
   On all IEEE754 machines, this is +Infinity.  */
#include <bits/huge_val.h>
#ifdef __USE_ISOC99
# include <bits/huge_valf.h>
# include <bits/huge_vall.h>

/* Get machine-dependent INFINITY value.  */
# include <bits/inf.h>

/* Get machine-dependent NAN value (returned for some domain errors).  */
# include <bits/nan.h>
#endif /* __USE_ISOC99 */

/* Get general and ISO C99 specific information.  */
#include <bits/mathdef.h>

/* The file <bits/mathcalls.h> contains the prototypes for all the
   actual math functions.  These macros are used for those prototypes,
   so we can easily declare each function as both `name' and `__name',
   and can declare the float versions `namef' and `__namef'.  */

#define __MATHCALL(function,suffix, args)	\
  __MATHDECL (_Mdouble_,function,suffix, args)
#define __MATHDECL(type, function,suffix, args) \
  __MATHDECL_1(type, function,suffix, args); \
  __MATHDECL_1(type, __CONCAT(__,function),suffix, args)
#define __MATHCALLX(function,suffix, args, attrib)	\
  __MATHDECLX (_Mdouble_,function,suffix, args, attrib)
#define __MATHDECLX(type, function,suffix, args, attrib) \
  __MATHDECL_1(type, function,suffix, args) __attribute__ (attrib); \
  __MATHDECL_1(type, __CONCAT(__,function),suffix, args) __attribute__ (attrib)
#define __MATHDECL_1(type, function,suffix, args) \
  extern type __MATH_PRECNAME(function,suffix) args __THROW

#define _Mdouble_ 		double
#define __MATH_PRECNAME(name,r)	__CONCAT(name,r)
# define _Mdouble_BEGIN_NAMESPACE __BEGIN_NAMESPACE_STD
# define _Mdouble_END_NAMESPACE   __END_NAMESPACE_STD
#include <bits/mathcalls.h>
#undef	_Mdouble_
#undef _Mdouble_BEGIN_NAMESPACE
#undef _Mdouble_END_NAMESPACE
#undef	__MATH_PRECNAME

#if defined __USE_MISC || defined __USE_ISOC99


/* Include the file of declarations again, this time using `float'
   instead of `double' and appending f to each function name.  */

# ifndef _Mfloat_
#  define _Mfloat_		float
# endif
# define _Mdouble_ 		_Mfloat_
# ifdef __STDC__
#  define __MATH_PRECNAME(name,r) name##f##r
# else
#  define __MATH_PRECNAME(name,r) name/**/f/**/r
# endif
# define _Mdouble_BEGIN_NAMESPACE __BEGIN_NAMESPACE_C99
# define _Mdouble_END_NAMESPACE   __END_NAMESPACE_C99
# include <bits/mathcalls.h>
# undef	_Mdouble_
# undef _Mdouble_BEGIN_NAMESPACE
# undef _Mdouble_END_NAMESPACE
# undef	__MATH_PRECNAME

# if (__STDC__ - 0 || __GNUC__ - 0) \
     && (!defined __NO_LONG_DOUBLE_MATH || defined __LDBL_COMPAT)
#  ifdef __LDBL_COMPAT

#   ifdef __USE_ISOC99
extern float __nldbl_nexttowardf (float __x, long double __y)
				  __THROW __attribute__ ((__const__));
#    ifdef __REDIRECT_NTH
extern float __REDIRECT_NTH (nexttowardf, (float __x, long double __y),
			     __nldbl_nexttowardf)
     __attribute__ ((__const__));
extern double __REDIRECT_NTH (nexttoward, (double __x, long double __y),
			      nextafter) __attribute__ ((__const__));
extern long double __REDIRECT_NTH (nexttowardl,
				   (long double __x, long double __y),
				   nextafter) __attribute__ ((__const__));
#    endif
#   endif

/* Include the file of declarations again, this time using `long double'
   instead of `double' and appending l to each function name.  */

#   undef __MATHDECL_1
#   define __MATHDECL_2(type, function,suffix, args, alias) \
  extern type __REDIRECT_NTH(__MATH_PRECNAME(function,suffix), \
			     args, alias)
#   define __MATHDECL_1(type, function,suffix, args) \
  __MATHDECL_2(type, function,suffix, args, __CONCAT(function,suffix))
#  endif

#  ifndef _Mlong_double_
#   define _Mlong_double_	long double
#  endif
#  define _Mdouble_ 		_Mlong_double_
#  ifdef __STDC__
#   define __MATH_PRECNAME(name,r) name##l##r
#  else
#   define __MATH_PRECNAME(name,r) name/**/l/**/r
#  endif
#  define _Mdouble_BEGIN_NAMESPACE __BEGIN_NAMESPACE_C99
#  define _Mdouble_END_NAMESPACE   __END_NAMESPACE_C99
#  include <bits/mathcalls.h>
#  undef _Mdouble_
# undef _Mdouble_BEGIN_NAMESPACE
# undef _Mdouble_END_NAMESPACE
#  undef __MATH_PRECNAME

# endif /* __STDC__ || __GNUC__ */

#endif	/* Use misc or ISO C99.  */
#undef	__MATHDECL_1
#undef	__MATHDECL
#undef	__MATHCALL


#if defined __USE_MISC || defined __USE_XOPEN
/* This variable is used by `gamma' and `lgamma'.  */
extern int signgam;
#endif


/* ISO C99 defines some generic macros which work on any data type.  */
#ifdef __USE_ISOC99

/* Get the architecture specific values describing the floating-point
   evaluation.  The following symbols will get defined:

    float_t	floating-point type at least as wide as `float' used
		to evaluate `float' expressions
    double_t	floating-point type at least as wide as `double' used
		to evaluate `double' expressions

    FLT_EVAL_METHOD
		Defined to
		  0	if `float_t' is `float' and `double_t' is `double'
		  1	if `float_t' and `double_t' are `double'
		  2	if `float_t' and `double_t' are `long double'
		  else	`float_t' and `double_t' are unspecified

    INFINITY	representation of the infinity value of type `float'

    FP_FAST_FMA
    FP_FAST_FMAF
    FP_FAST_FMAL
		If defined it indicates that the `fma' function
		generally executes about as fast as a multiply and an add.
		This macro is defined only iff the `fma' function is
		implemented directly with a hardware multiply-add instructions.

    FP_ILOGB0	Expands to a value returned by `ilogb (0.0)'.
    FP_ILOGBNAN	Expands to a value returned by `ilogb (NAN)'.

    DECIMAL_DIG	Number of decimal digits supported by conversion between
		decimal and all internal floating-point formats.

*/

/* All floating-point numbers can be put in one of these categories.  */
enum
  {
    FP_NAN,
# define FP_NAN FP_NAN
    FP_INFINITE,
# define FP_INFINITE FP_INFINITE
    FP_ZERO,
# define FP_ZERO FP_ZERO
    FP_SUBNORMAL,
# define FP_SUBNORMAL FP_SUBNORMAL
    FP_NORMAL
# define FP_NORMAL FP_NORMAL
  };

/* Return number of classification appropriate for X.  */
# ifdef __NO_LONG_DOUBLE_MATH
#  define fpclassify(x) \
     (sizeof (x) == sizeof (float) ? __fpclassifyf (x) : __fpclassify (x))
# else
#  define fpclassify(x) \
     (sizeof (x) == sizeof (float)					      \
      ? __fpclassifyf (x)						      \
      : sizeof (x) == sizeof (double)					      \
      ? __fpclassify (x) : __fpclassifyl (x))
# endif

/* Return nonzero value if sign of X is negative.  */
# ifdef __NO_LONG_DOUBLE_MATH
#  define signbit(x) \
     (sizeof (x) == sizeof (float) ? __signbitf (x) : __signbit (x))
# else
#  define signbit(x) \
     (sizeof (x) == sizeof (float)					      \
      ? __signbitf (x)							      \
      : sizeof (x) == sizeof (double)					      \
      ? __signbit (x) : __signbitl (x))
# endif

/* Return nonzero value if X is not +-Inf or NaN.  */
# ifdef __NO_LONG_DOUBLE_MATH
#  define isfinite(x) \
     (sizeof (x) == sizeof (float) ? __finitef (x) : __finite (x))
# else
#  define isfinite(x) \
     (sizeof (x) == sizeof (float)					      \
      ? __finitef (x)							      \
      : sizeof (x) == sizeof (double)					      \
      ? __finite (x) : __finitel (x))
# endif

/* Return nonzero value if X is neither zero, subnormal, Inf, nor NaN.  */
# define isnormal(x) (fpclassify (x) == FP_NORMAL)

/* Return nonzero value if X is a NaN.  We could use `fpclassify' but
   we already have this functions `__isnan' and it is faster.  */
# ifdef __NO_LONG_DOUBLE_MATH
#  define isnan(x) \
     (sizeof (x) == sizeof (float) ? __isnanf (x) : __isnan (x))
# else
#  define isnan(x) \
     (sizeof (x) == sizeof (float)					      \
      ? __isnanf (x)							      \
      : sizeof (x) == sizeof (double)					      \
      ? __isnan (x) : __isnanl (x))
# endif

/* Return nonzero value is X is positive or negative infinity.  */
# ifdef __NO_LONG_DOUBLE_MATH
#  define isinf(x) \
     (sizeof (x) == sizeof (float) ? __isinff (x) : __isinf (x))
# else
#  define isinf(x) \
     (sizeof (x) == sizeof (float)					      \
      ? __isinff (x)							      \
      : sizeof (x) == sizeof (double)					      \
      ? __isinf (x) : __isinfl (x))
# endif

/* Bitmasks for the math_errhandling macro.  */
# define MATH_ERRNO	1	/* errno set by math functions.  */
# define MATH_ERREXCEPT	2	/* Exceptions raised by math functions.  */

/* By default all functions support both errno and exception handling.
   In gcc's fast math mode and if inline functions are defined this
   might not be true.  */
# ifndef __FAST_MATH__
#  define math_errhandling	(MATH_ERRNO | MATH_ERREXCEPT)
# endif

#endif /* Use ISO C99.  */

#ifdef	__USE_MISC
/* Support for various different standard error handling behaviors.  */
typedef enum
{
  _IEEE_ = -1,	/* According to IEEE 754/IEEE 854.  */
  _SVID_,	/* According to System V, release 4.  */
  _XOPEN_,	/* Nowadays also Unix98.  */
  _POSIX_,
  _ISOC_	/* Actually this is ISO C99.  */
} _LIB_VERSION_TYPE;

/* This variable can be changed at run-time to any of the values above to
   affect floating point error handling behavior (it may also be necessary
   to change the hardware FPU exception settings).  */
extern _LIB_VERSION_TYPE _LIB_VERSION;
#endif


#ifdef __USE_SVID
/* In SVID error handling, `matherr' is called with this description
   of the exceptional condition.

   We have a problem when using C++ since `exception' is a reserved
   name in C++.  */
# ifdef __cplusplus
struct __exception
# else
struct exception
# endif
  {
    int type;
    char *name;
    double arg1;
    double arg2;
    double retval;
  };

# ifdef __cplusplus
extern int matherr (struct __exception *__exc) throw ();
# else
extern int matherr (struct exception *__exc);
# endif

# define X_TLOSS	1.41484755040568800000e+16

/* Types of exceptions in the `type' field.  */
# define DOMAIN		1
# define SING		2
# define OVERFLOW	3
# define UNDERFLOW	4
# define TLOSS		5
# define PLOSS		6

/* SVID mode specifies returning this large value instead of infinity.  */
# define HUGE		3.40282347e+38F

#else	/* !SVID */

# ifdef __USE_XOPEN
/* X/Open wants another strange constant.  */
#  define MAXFLOAT	3.40282347e+38F
# endif

#endif	/* SVID */


/* Some useful constants.  */
#if defined __USE_BSD || defined __USE_XOPEN
# define M_E		2.7182818284590452354	/* e */
# define M_LOG2E	1.4426950408889634074	/* log_2 e */
# define M_LOG10E	0.43429448190325182765	/* log_10 e */
# define M_LN2		0.69314718055994530942	/* log_e 2 */
# define M_LN10		2.30258509299404568402	/* log_e 10 */
# define M_PI		3.14159265358979323846	/* pi */
# define M_PI_2		1.57079632679489661923	/* pi/2 */
# define M_PI_4		0.78539816339744830962	/* pi/4 */
# define M_1_PI		0.31830988618379067154	/* 1/pi */
# define M_2_PI		0.63661977236758134308	/* 2/pi */
# define M_2_SQRTPI	1.12837916709551257390	/* 2/sqrt(pi) */
# define M_SQRT2	1.41421356237309504880	/* sqrt(2) */
# define M_SQRT1_2	0.70710678118654752440	/* 1/sqrt(2) */
#endif

/* The above constants are not adequate for computation using `long double's.
   Therefore we provide as an extension constants with similar names as a
   GNU extension.  Provide enough digits for the 128-bit IEEE quad.  */
#ifdef __USE_GNU
# define M_El		2.7182818284590452353602874713526625L  /* e */
# define M_LOG2El	1.4426950408889634073599246810018921L  /* log_2 e */
# define M_LOG10El	0.4342944819032518276511289189166051L  /* log_10 e */
# define M_LN2l		0.6931471805599453094172321214581766L  /* log_e 2 */
# define M_LN10l	2.3025850929940456840179914546843642L  /* log_e 10 */
# define M_PIl		3.1415926535897932384626433832795029L  /* pi */
# define M_PI_2l	1.5707963267948966192313216916397514L  /* pi/2 */
# define M_PI_4l	0.7853981633974483096156608458198757L  /* pi/4 */
# define M_1_PIl	0.3183098861837906715377675267450287L  /* 1/pi */
# define M_2_PIl	0.6366197723675813430755350534900574L  /* 2/pi */
# define M_2_SQRTPIl	1.1283791670955125738961589031215452L  /* 2/sqrt(pi) */
# define M_SQRT2l	1.4142135623730950488016887242096981L  /* sqrt(2) */
# define M_SQRT1_2l	0.7071067811865475244008443621048490L  /* 1/sqrt(2) */
#endif


/* When compiling in strict ISO C compatible mode we must not use the
   inline functions since they, among other things, do not set the
   `errno' variable correctly.  */
#if defined __STRICT_ANSI__ && !defined __NO_MATH_INLINES
# define __NO_MATH_INLINES	1
#endif

#if defined __USE_ISOC99 && __GNUC_PREREQ(2,97)
/* ISO C99 defines some macros to compare number while taking care for
   unordered numbers.  Many FPUs provide special instructions to support
   these operations.  Generic support in GCC for these as builtins went
   in before 3.0.0, but not all cpus added their patterns.  We define
   versions that use the builtins here, and <bits/mathinline.h> will
   undef/redefine as appropriate for the specific GCC version in use.  */
# define isgreater(x, y)	__builtin_isgreater(x, y)
# define isgreaterequal(x, y)	__builtin_isgreaterequal(x, y)
# define isless(x, y)		__builtin_isless(x, y)
# define islessequal(x, y)	__builtin_islessequal(x, y)
# define islessgreater(x, y)	__builtin_islessgreater(x, y)
# define isunordered(u, v)	__builtin_isunordered(u, v)
#endif

/* Get machine-dependent inline versions (if there are any).  */
#ifdef __USE_EXTERN_INLINES
# include <bits/mathinline.h>
#endif

#ifdef __USE_ISOC99
/* If we've still got undefined comparison macros, provide defaults.  */

/* Return nonzero value if X is greater than Y.  */
# ifndef isgreater
#  define isgreater(x, y) \
  (__extension__							      \
   ({ __typeof__(x) __x = (x); __typeof__(y) __y = (y);			      \
      !isunordered (__x, __y) && __x > __y; }))
# endif

/* Return nonzero value if X is greater than or equal to Y.  */
# ifndef isgreaterequal
#  define isgreaterequal(x, y) \
  (__extension__							      \
   ({ __typeof__(x) __x = (x); __typeof__(y) __y = (y);			      \
      !isunordered (__x, __y) && __x >= __y; }))
# endif

/* Return nonzero value if X is less than Y.  */
# ifndef isless
#  define isless(x, y) \
  (__extension__							      \
   ({ __typeof__(x) __x = (x); __typeof__(y) __y = (y);			      \
      !isunordered (__x, __y) && __x < __y; }))
# endif

/* Return nonzero value if X is less than or equal to Y.  */
# ifndef islessequal
#  define islessequal(x, y) \
  (__extension__							      \
   ({ __typeof__(x) __x = (x); __typeof__(y) __y = (y);			      \
      !isunordered (__x, __y) && __x <= __y; }))
# endif

/* Return nonzero value if either X is less than Y or Y is less than X.  */
# ifndef islessgreater
#  define islessgreater(x, y) \
  (__extension__							      \
   ({ __typeof__(x) __x = (x); __typeof__(y) __y = (y);			      \
      !isunordered (__x, __y) && (__x < __y || __y < __x); }))
# endif

/* Return nonzero value if arguments are unordered.  */
# ifndef isunordered
#  define isunordered(u, v) \
  (__extension__							      \
   ({ __typeof__(u) __u = (u); __typeof__(v) __v = (v);			      \
      fpclassify (__u) == FP_NAN || fpclassify (__v) == FP_NAN; }))
# endif

#endif

__END_DECLS


#endif /* math.h  */
