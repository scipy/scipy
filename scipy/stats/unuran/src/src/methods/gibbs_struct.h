/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: gibbs_struct.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method GIBBS                              *
 *         (Markov Chain - GIBBS sampler)                                    *
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

struct unur_gibbs_par { 
  int thinning;             /* thinning factor for generated chain           */
  int burnin;               /* length of burn-in for chain                   */
  double  c_T;              /* parameter c for transformation T_c            */
  const double *x0;         /* starting point of chain                       */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_gibbs_gen {
  int dim;                  /* dimension of distribution                     */
  int thinning;             /* thinning factor for generated chain           */
  double  c_T;              /* parameter c for transformation T_c            */

  double *state;            /* state of chain / current point                */

  struct unur_distr *distr_condi; /* conditional distribution                */

  int    coord;             /* current coordinate used for GIBBS chain       */
  double *direction;        /* working array for random direction            */

  int burnin;               /* length of burn-in for chain                   */
  double *x0;               /* starting point of chain                       */
};

/*---------------------------------------------------------------------------*/

