/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: vempk_struct.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method VEMPK                              *
 *         ((Vector) EMPirical distribution with Kernel smoothing)           *
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

struct unur_vempk_par {
  /* the observed sample is stored in the distribution object */
  double  smoothing;   /* determines how "smooth" the estimated density will be */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_vempk_gen {
  double *observ;      /* pointer to the array of the observations           */
  int     n_observ;    /* number of observations                             */
  int     dim;         /* dimension of distribution                          */

  UNUR_GEN *kerngen;   /* random variate generator for kernel                */

  double  smoothing;   /* determines how "smooth" the estimated density will be */

  double  hopt;        /* for bandwidth selection                            */
  double  hact;        /* actually used value for bandwith                   */
  double  corfac;      /* correction for variance correction                 */
  double *xbar;        /* mean vector of sample, for variance correction     */
};

/*---------------------------------------------------------------------------*/

