/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: tdr_struct.h                                                      *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method TDR                                *
 *         (Transformed Density Rejection)                                   *
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

struct unur_tdr_par { 

  double  guide_factor;         /* relative size of guide table              */

  const double *starting_cpoints; /* pointer to array of starting points     */
  int     n_starting_cpoints;   /* number of construction points at start    */

  const double *percentiles;    /* percentiles of hat for c. points of new hat  */
  int n_percentiles;            /* number of percentiles                     */
  int retry_ncpoints;           /* number of cpoints for second trial of reinit */

  int     max_ivs;              /* maximum number of intervals               */
  double  max_ratio;            /* bound for ratio r_n = Atotal / Asqueeze   */
  double  bound_for_adding;     /* lower bound for relative area             */

  double  c_T;                  /* parameter c for transformation T_c        */           
  double  darsfactor;           /* factor for derandomized ARS               */
  int     darsrule;             /* rule for finding splitting points in DARS */
};

/*---------------------------------------------------------------------------*/
/* store data for segments                                                   */

struct unur_tdr_interval {

  double  x;                    /* (left) construction point (cp)            */
  double  fx;                   /* value of PDF at cp                        */ 
  double  Tfx;                  /* value of transformed PDF at cp            */ 
  double  dTfx;                 /* derivative of transformed PDF at cp       */
  double  sq;                   /* slope of transformed squeeze in interval  */
  double  ip;                   /* intersection point between two tangents   */
  double  fip;                  /* value of PDF at ip (for PS and IA only)   */
                                 
  double  Acum;                 /* cumulated area of intervals               */
  double  Ahat;                 /* area below hat                            */
  double  Ahatr;                /* area below hat on right side              */
  double  Asqueeze;             /* area squeeze                              */

  struct unur_tdr_interval *next; /* pointer to next interval in list        */
  struct unur_tdr_interval *prev; /* pointer to previous interval in list    
				     (for PS and IA only)                    */

#ifdef UNUR_COOKIES
  unsigned cookie;              /* magic cookie                              */
#endif
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_tdr_gen { 

  double  Atotal;               /* area below hat                            */
  double  Asqueeze;             /* area below squeeze                        */

  double  c_T;                  /* parameter c for transformation T_c        */           
  double  Umin, Umax;           /* bounds for iid random variable in respect to
				   the given (truncated) domain of the distr.*/

  struct unur_tdr_interval *iv; /* pointer to linked list of intervals       */
  int     n_ivs;                /* number of intervals                       */
  int     max_ivs;              /* maximum number of intervals               */
  double  max_ratio;            /* bound for ratio r_n = Atotal / Asqueeze   */
  double  bound_for_adding;     /* lower bound for relative area             */

  struct unur_tdr_interval **guide; /* pointer to guide table                */
  int     guide_size;           /* size of guide table                       */
  double  guide_factor;         /* relative size of guide table              */

  double  center;               /* approximate location of mode              */

  double *starting_cpoints;     /* pointer to array of starting points       */
  int     n_starting_cpoints;   /* number of construction points at start    */

  double *percentiles;          /* percentiles of hat for c. points of new hat  */
  int n_percentiles;            /* number of percentiles                        */
  int retry_ncpoints;           /* number of cpoints for second trial of reinit */

  double  darsfactor;           /* factor for derandomized ARS               */
  int     darsrule;             /* rule for finding splitting points in DARS */
#ifdef UNUR_ENABLE_INFO
  int     max_ivs_info;         /* maximum number of intervals (as given)    */
#endif
};

/*---------------------------------------------------------------------------*/
