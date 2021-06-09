/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: umath_source.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares macros, constants, structures, function prototypes, etc. *
 *         for using mathematics in UNU.RAN.                                 *
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
#ifndef UMATH_SOURCE_H_SEEN
#define UMATH_SOURCE_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Prototypes for various functions used in library                          */  

/* Compute a "mean" defined as combibation of arithmetic and harmonic mean   */
double _unur_arcmean( double x0, double x1 );

/*---------------------------------------------------------------------------*/
/* Macros                                                                    */

#define _unur_min(x,y)   (((x)<(y)) ? (x) : (y))
#define _unur_max(x,y)   (((x)>(y)) ? (x) : (y))

/*---------------------------------------------------------------------------*/
#endif  /* UMATH_SOURCE_H_SEEN */
/*---------------------------------------------------------------------------*/






