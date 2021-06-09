/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: srou_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method SROU                               *
 *         (Simple universal generator, Ratio-Of-Uniforms method)            *
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
/* Information for constructing the generator                                */

struct unur_srou_par { 
  double  r;                 /* parameter for power transformation           */
  double  Fmode;             /* cdf at mode                                  */
  double  um;                /* square root of pdf at mode                   */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_srou_gen { 
  double  um;                /* height of rectangle: square root of f(mode)  */
  double  vl, vr;            /* left and right boundary of rectangle         */
  double  xl, xr;            /* ratios vl/um and vr/um                       */
  double  Fmode;             /* cdf at mode                                  */

  /* parameters for generalized SROU                                         */
  double  r;                 /* parameter for power transformation           */
  double  p;                 /* construction point for bounding curve        */
  double  a, b;              /* parameters for bounding curve                */
  double  log_ab;            /* parameter for bounding curve                 */
};

/*---------------------------------------------------------------------------*/
