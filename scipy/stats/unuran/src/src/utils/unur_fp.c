/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_fp.c                                                         *
 *                                                                           *
 *   miscelleanous routines for floating point arithmetic                    *
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

#include <unur_source.h>

/*---------------------------------------------------------------------------*/

int
_unur_FP_cmp( double x1, double x2, double eps)
     /*----------------------------------------------------------------------*/
     /* Compare two floats:                                                  */
     /*   x1 eq x2 iff |x1-x2| <= min(|x1|,|x2) * eps                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x1  ... double                                                     */
     /*   x2  ... double                                                     */
     /*   eps ... maximal relative deviation                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   -1 if x1 < x2                                                      */
     /*    0 if x1 eq x2                                                     */
     /*   +1 if x1 > x2                                                      */
     /*                                                                      */
     /* remark:                                                              */
     /*   This is similar to Knuth's algorithm. However, we use              */
     /*   instead of max(|x1|,|x2).                                          */
     /*   We also have to deal with +-INFINITY correctly.                    */
     /*                                                                      */
     /*   (For an implementation of Knuth's algorithm see                    */
     /*    fcmp 1.2.2 Copyright (c) 1998-2000 Theodore C. Belding            */
     /*    University of Michigan Center for the Study of Complex Systems    */
     /*    Ted.Belding@umich.edu)                                            */
     /*----------------------------------------------------------------------*/
{
  double fx1 = (x1>=0.) ? x1 : -x1;
  double fx2 = (x2>=0.) ? x2 : -x2;
  double delta = eps * _unur_min(fx1,fx2);
  double difference = x1 - x2;

  /* we have to take care about INFINITY */
  if (_unur_isinf(delta)) {
    delta = eps * DBL_MAX;
  }

  /* denormalized numbers (close to zero) may cause problems. */
  /* so we check this special case.                           */
  if (fx1 <= 2.*DBL_MIN && fx2 <= 2.*DBL_MIN)
    return 0;

  if (difference > delta)       /* x1 > x2 */
    return +1;
  else if (difference < -delta) /* x1 < x2 */
    return -1;
  else                          /* -delta <= difference <= delta */
    return 0;                   /* x1 ~=~ x2 */

} /* end of _unur_FP_cmp() */

/*---------------------------------------------------------------------------*/

#ifndef _unur_iszero
int _unur_iszero (const double x)
     /*----------------------------------------------------------------------*/
     /* Check whether x is equal to 0.0                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x ... floating-point number                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   TRUE   if x == 0                                                   */
     /*   FALSE  otherwise                                                   */
     /*----------------------------------------------------------------------*/
{
  return (x==0.);
} /* end of _unur_iszero() */
#endif

/*---------------------------------------------------------------------------*/

#ifndef _unur_isone
int _unur_isone (const double x)
     /*----------------------------------------------------------------------*/
     /* Check whether x is equal to 1.0                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x ... floating-point number                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   TRUE   if x == 1                                                   */
     /*   FALSE  otherwise                                                   */
     /*----------------------------------------------------------------------*/
{
  return (x==1.);
} /* end of _unur_isone() */
#endif

/*---------------------------------------------------------------------------*/

#ifndef _unur_isfsame
int _unur_isfsame (const double x, const double y)
     /*----------------------------------------------------------------------*/
     /* Check whether x == y                                                 */
     /* (Unsafe version that should only be used for powers of 2)            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x ... floating-point number                                        */
     /*   y ... floating-point number                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   TRUE   if x == y                                                   */
     /*   FALSE  otherwise                                                   */
     /*----------------------------------------------------------------------*/
{
  return (x==y);
} /* end of _unur_isfsame() */
#endif

/*---------------------------------------------------------------------------*/

int
_unur_isfinite (const double x)
     /*----------------------------------------------------------------------*/
     /* Check whether x is a finite number.                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x ... floating-point number                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   TRUE   if x is finite                                              */
     /*   FALSE  otherwise (i.e., if x is +/- infinity or NaN)               */
     /*----------------------------------------------------------------------*/
{
#if HAVE_DECL_ISFINITE
  return (isfinite(x) ? TRUE : FALSE);
#elif defined(_MSC_VER) /* Microsoft Visual C++ */
  return (_finite(x) ? TRUE : FALSE);
#elif HAVE_IEEE_COMPARISONS
  if (x < INFINITY && x > -INFINITY)
    return TRUE;
  else
    return FALSE;
#else
# error
# error +--------------------------------------------+
# error ! Sorry, Cannot handle INFINITY correctly! . !
# error ! Please contact <unuran@statmath.wu.ac.at>. !
# error +--------------------------------------------+
# error
#endif
} /* end of _unur_isfinite() */

/*---------------------------------------------------------------------------*/

int
_unur_isnan (const double x)
     /*----------------------------------------------------------------------*/
     /*  Check whether x is a NaN (not a number) value.                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x ... floating-point number                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   TRUE   if x is NaN                                                 */
     /*   FALSE  otherwise                                                   */
     /*----------------------------------------------------------------------*/
{
#if HAVE_DECL_ISNAN
  return (isnan(x) ? TRUE : FALSE);
#elif defined(_MSC_VER) /* Microsoft Visual C++ */
  return (_isnan(x) ? TRUE : FALSE);
#elif HAVE_IEEE_COMPARISONS
  return ((x!=x) ? TRUE : FALSE);
#else
# error
# error +--------------------------------------------+
# error ! Sorry, Cannot handle NaN correctly! ...... !
# error ! Please contact <unuran@statmath.wu.ac.at>. !
# error +--------------------------------------------+
# error
#endif
} /* end of _unur_isnan() */

/*---------------------------------------------------------------------------*/

int
_unur_isinf (const double x)
     /*----------------------------------------------------------------------*/
     /*  Check whether x is infinity.                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x ... floating-point number                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   -1  ... if x is -infinity                                          */
     /*    1  ... if x is +infinity                                          */
     /*    0  ... otherwise                                                  */
     /*----------------------------------------------------------------------*/
{
#if HAVE_DECL_ISINF
  return isinf(x);
#elif defined(_MSC_VER) /* Microsoft Visual C++ */
  int fpc = _fpclass(x);

  if (fpc == _FPCLASS_PINF)
    return +1;
  else if (fpc == _FPCLASS_NINF)
    return -1;
  else 
    return 0;
#elif HAVE_IEEE_COMPARISONS
  if (x>=INFINITY)
    return 1;
  else if (x<=-INFINITY)
    return -1;
  else
    return 0;
#else
# error
# error +--------------------------------------------+
# error ! Sorry, Cannot handle INFINITY correctly! . !
# error ! Please contact <unuran@statmath.wu.ac.at>. !
# error +--------------------------------------------+
# error
#endif
} /* end of _unur_isinf() */

/*---------------------------------------------------------------------------*/
