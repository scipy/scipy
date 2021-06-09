/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_fp_source.h                                                  *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares macros for floating point arithmetic.                    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifndef UNUR_FP_SOURCE_H_SEEN
#define UNUR_FP_SOURCE_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Comparisons                                                               */

int _unur_FP_cmp( double x1, double x2, double eps);
/* compare two floats:                                                       */
/*   x1 eq x2 iff |x1-x2| <= min(|x1|,|x2) * eps                             */
/* return:                                                                   */
/*   -1 if x1 < x2                                                           */
/*    0 if x1 eq x2                                                          */
/*   +1 if x1 > x2                                                           */

/* macros for different levels of accuracy                                   */
#define _unur_FP_cmp_same(a,b) (_unur_FP_cmp((a),(b),DBL_EPSILON))
#define _unur_FP_cmp_equal(a,b) (_unur_FP_cmp((a),(b),UNUR_EPSILON))
#define _unur_FP_cmp_approx(a,b) (_unur_FP_cmp((a),(b),UNUR_SQRT_DBL_EPSILON))

/* a == b (except precision bit) */
#define _unur_FP_same(a,b) (_unur_FP_cmp((a),(b),DBL_EPSILON)==0)

/* a == b */
#define _unur_FP_equal(a,b) (_unur_FP_cmp((a),(b),UNUR_EPSILON)==0)

/* a is approximately equal to b */
#define _unur_FP_approx(a,b) (_unur_FP_cmp((a),(b),UNUR_SQRT_DBL_EPSILON)==0)

/* a < b */
#define _unur_FP_less(a,b) ((_unur_FP_cmp((a),(b),UNUR_EPSILON)<0) ? TRUE : FALSE)

/* a > b */
#define _unur_FP_greater(a,b) ((_unur_FP_cmp((a),(b),UNUR_EPSILON)>0) ? TRUE : FALSE)

/*---------------------------------------------------------------------------*/
/* Comparing floating point with == or != is unsafe.                         */
/* However, we assume that comparing with 0.0 and powers of 2.0 is safe.     */
/* Thus we use the followig functions to mark these "safe" comparisons in    */
/* the code and thus we can use the GCC to detect all other comparisons.     */
/* For the latter _unur_FP_cmp_same() must be used.                          */

/* defined as macros */
#define _unur_iszero(x)     ((x)==0.0)
#define _unur_isone(x)      ((x)==1.0)
#define _unur_isfsame(x,y)  ((x)==(y))

/* defined as functions: only used for debugging.                            */
/* (it switches off GCC warnings)                                            */
#ifndef _unur_iszero
int _unur_iszero (const double x);
#endif

#ifndef _unur_isone
int _unur_isone (const double x);
#endif

#ifndef _unur_isfsame
int _unur_isfsame (const double x, const double y);
#endif

/*---------------------------------------------------------------------------*/
/* Infinity and NaN (Not a Number)                                           */

/* wrapper / replacement for corresponding C99 functions */
int _unur_isfinite (const double x);
int _unur_isnan (const double x);
int _unur_isinf (const double x);

/*---------------------------------------------------------------------------*/
/* Other checks for infinity                                                 */

/* +oo */
#define _unur_FP_is_infinity(a)  ((a) >= INFINITY)

/* -oo */
#define _unur_FP_is_minus_infinity(a)  ((a) <= -INFINITY)

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_FP_SOURCE_H_SEEN */
/*---------------------------------------------------------------------------*/






