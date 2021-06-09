/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: umath.c                                                           *
 *                                                                           *
 *   miscelleanous mathematical routines                                     *
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

/* If the macro INFINITY is not already defined we store infinity in the     */
/* global variable INFINITY.                                                 */
#ifndef INFINITY
const double INFINITY = 1.0 / 0.0;
#endif

/*---------------------------------------------------------------------------*/

#define ARCMEAN_HARMONIC   (1.e3)   /* use harmonic mean when abs larger than this value */
#define ARCMEAN_ARITHMETIC (1.e-6)  /* use harmonic mean when abs larger than this value */

double
_unur_arcmean( double x0, double x1 )
     /*----------------------------------------------------------------------*/
     /* compute "arctan mean" of two numbers.                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x0, x1 ... two numbers                                             */
     /*                                                                      */
     /* return:                                                              */
     /*   mean                                                               */
     /*                                                                      */
     /* comment:                                                             */
     /*   "arctan mean" = tan(0.5*(arctan(x0)+arctan(x1)))                   */
     /*                                                                      */
     /*   a combination of arithmetical mean (for x0 and x1 close to 0)      */
     /*   and the harmonic mean (for |x0| and |x1| large).                   */
     /*----------------------------------------------------------------------*/
{
  double a0,a1;
  double r;

  /* we need x0 < x1 */
  if (x0>x1) {double tmp = x0; x0=x1; x1=tmp;}

  if (x1 < -ARCMEAN_HARMONIC || x0 > ARCMEAN_HARMONIC)
    /* use harmonic mean */
    return (2./(1./x0 + 1./x1));

  a0 = (x0<=-INFINITY) ? -M_PI/2. : atan(x0);
  a1 = (x1>= INFINITY) ?  M_PI/2. : atan(x1);

  if (fabs(a0-a1) < ARCMEAN_ARITHMETIC)
    /* use arithmetic mean */
    r = 0.5*x0 + 0.5*x1;

  else
    /* use "arc mean" */
    r = tan((a0 + a1)/2.);

  return r;
} /* end of _unur_arcmean() */

/*---------------------------------------------------------------------------*/




