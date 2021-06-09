/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: empk_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method EMPK                               *
 *         (EMPirical distribution with Kernel smoothing)                    *
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

struct unur_empk_par {
  /* the observed sample is stored in the distribution object */

  const UNUR_GEN *kerngen;  /* random variate generator for kernel
			     (given by user)                                */
  UNUR_GEN *kernel;    /* random variate generator for kernel
			  (provided by UNURAN)                              */

  double  alpha;       /* alpha is used to compute the optimal bandwidth from
			  the point of view of minimizing the mean integrated
			  square error (MISE).
			  alfa depends on the type of kernel being used.     */

  double  beta;        /* beta is used to compute the optimal bandwidth from
			  the point of view of minimizing the mean integrated
			  square error (MISE).
			  beta depends on the (unknown) distribution of the
			  sampled data points. Thus its value contains some
			  guess on this distribution.                        */

  double  smoothing;   /* determines how "smooth" the estimated density will be */

  double  kernvar;     /* variance of used kernel, only used if varcor == 1  */

};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_empk_gen {
  double *observ;      /* pointer to the array of the observations           */
  int     n_observ;    /* number of observations                             */

  UNUR_GEN *kerngen;   /* random variate generator for kernel                */

  double  smoothing;   /* determines how "smooth" the estimated density will be */
  double  kernvar;     /* variance of used kernel, only used if varcor == 1  */

  double  bwidth;      /* bandwidth for kernel density estimation            */
  double  bwidth_opt;  /* optimal bandwith for kernel density estimation     */
  double  mean_observ; /* mean of observed data                              */
  double  stddev_observ; /* standard deviation of oberved data               */
  double  sconst;      /* constant used for variance corrected version 
			  of kernel method                                   */
  double  alpha;       /* parameter for kernel width                         */
  double  beta;        /* parameter for kernel width                         */
};

/*---------------------------------------------------------------------------*/
