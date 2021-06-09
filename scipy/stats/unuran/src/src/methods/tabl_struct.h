/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: tabl_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method TABL                               *
 *         (Ahren's TABLe method: piecewise constant hat)                    *
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

struct unur_tabl_par { 

  const double *slopes; /* slopes <a_i,b_i>, i.e.\ f(a_i) >= f(b_i)          */
  int     n_slopes;     /* number of slopes                                  */
  double  bleft;        /* left border of the domain                         */
  double  bright;       /* right border of the domain                        */

  int     max_ivs;      /* maximum number of intervals                       */
  double  max_ratio;    /* limit for ratio r_n = A(squeeze) / A(hat)         */

  const double *cpoints; /* pointer to array of points for constr. slopes    */
  int     n_cpoints;    /* number of points for constructing slopes          */
  int     n_stp;        /* number of points for hat for start                */
  double  area_fract;   /* parameter for equal area rule                     */
  double  darsfactor;   /* factor for (derandomized) ARS                     */

  double  guide_factor; /* relative size of guide table                      */
};

/*---------------------------------------------------------------------------*/
/* storing information about generator                                       */

struct unur_tabl_interval {

  double  xmax;         /* maximum of pdf in interval                        */
  double  fmax;         /* maximal value of pdf in interval                  */
  double  xmin;         /* minimum of pdf in interval                        */
  double  fmin;         /* minimal value of pdf in interval                  */

  double  Ahat;         /* area of total bar (below hat)                     */
  double  Asqueeze;     /* area of bar below squeeze                         */
  double  Acum;         /* cumulated area of bars                            */

  struct unur_tabl_interval *next;  /* pointer to next element in list       */

#ifdef UNUR_COOKIES
  unsigned cookie;      /* magic cookie                                      */
#endif
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_tabl_gen { 
  double  Atotal;               /* total area below hat                      */
  double  Asqueeze;             /* area of squeeze polygon                   */

  double  bleft;                /* left boundary of domain                   */
  double  bright;               /* right boundary of domain                  */

  struct unur_tabl_interval **guide; /* pointer to guide table               */
  int     guide_size;           /* size of guide table                       */
  double  guide_factor;         /* relative size of guide table              */

  double  Umin, Umax;           /* bounds for iid random variable in respect to
				   the given (truncated) domain of the distr.*/

  struct unur_tabl_interval *iv;     /* pointer to linked list of intervals  */
  int     n_ivs;                /* number of intervals                       */
  int     max_ivs;              /* maximum number of intervals               */
  double  max_ratio;            /* limit for ratio r_n = A(squeeze) / A(hat) */

  double  darsfactor;           /* factor for (derandomized) ARS             */
#ifdef UNUR_ENABLE_INFO
  int     max_ivs_info;         /* maximum number of intervals (as given)    */
#endif
};

/*---------------------------------------------------------------------------*/
