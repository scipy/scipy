/* Copyright (C) 2001, 2004 Free Software Foundation, Inc.
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

#if !defined _MATH_H && !defined _COMPLEX_H
# error "Never use <bits/mathdef.h> directly; include <math.h> instead"
#endif

#if defined __USE_ISOC99 && defined _MATH_H && !defined _MATH_H_MATHDEF
# define _MATH_H_MATHDEF	1

# include <bits/wordsize.h>

# if __WORDSIZE == 64 || (defined __FLT_EVAL_METHOD__ && __FLT_EVAL_METHOD__ == 0)
/* The x86-64 architecture computes values with the precission of the
   used type.  Similarly for -m32 -mfpmath=sse.  */
typedef float float_t;		/* `float' expressions are evaluated as `float'.  */
typedef double double_t;	/* `double' expressions are evaluated
				   as `double'.  */
# else
/* The ix87 FPUs evaluate all values in the 80 bit floating-point format
   which is also available for the user as `long double'.  Therefore we
   define:  */
typedef long double float_t;	/* `float' expressions are evaluated as
				   `long double'.  */
typedef long double double_t;	/* `double' expressions are evaluated as
				   `long double'.  */
# endif

/* The values returned by `ilogb' for 0 and NaN respectively.  */
# define FP_ILOGB0	(-2147483647 - 1)
# define FP_ILOGBNAN	(-2147483647 - 1)

#endif	/* ISO C99 */
