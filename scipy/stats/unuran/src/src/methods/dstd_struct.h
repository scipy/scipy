/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: dstd_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method DSTD                               *
 *         (wrapper for special generators for                               *
 *         Discrete STanDard distributions)                                  *
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

struct unur_dstd_par { 
  int dummy;              /* no special parameters                           */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_dstd_gen { 
  double *gen_param;      /* parameters for the generator                    */
  int     n_gen_param;    /* number of parameters for the generator          */

  int    *gen_iparam;     /* integer parameters for generator                */
  int     n_gen_iparam;   /* number of integer parameters for the generator  */

  double  Umin;           /* cdf at left boundary of domain                  */
  double  Umax;           /* cdf at right boundary of domain                 */
  int  is_inversion;      /* indicate whether method is inversion method     */     
  const char *sample_routine_name; /* name of sampling routine               */
};

/*---------------------------------------------------------------------------*/
