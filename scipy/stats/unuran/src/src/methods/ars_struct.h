/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: ars_struct.h                                                      *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method ARS                                *
 *         (Adaptive Rejection Sampling - Gilks & Wild)                      *
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

struct unur_ars_par { 

  const double *starting_cpoints; /* pointer to array of starting points     */
  int n_starting_cpoints;         /* number of construction points at start  */
  const double *percentiles; /* percentiles of hat for c. points of new hat  */
  int n_percentiles;         /* number of percentiles                        */
  int retry_ncpoints;        /* number of cpoints for second trial of reinit */

  int max_ivs;               /* maximum number of intervals                  */
  int max_iter;              /* maximum number of iterations                 */
};

/*---------------------------------------------------------------------------*/
/* store data for segments                                                   */

struct unur_ars_interval {

  double  x;              /* (left hand side) construction point (cp)        */
  double  logfx;          /* value of logPDF at cp                           */
  double  dlogfx;         /* derivative of logPDF at cp                      */
  double  sq;             /* slope of transformed squeeze in interval        */

  double  Acum;           /* cumulated area of intervals                     */
  double  logAhat;        /* log of area below hat                           */
  double  Ahatr_fract;    /* fraction of area below hat on r.h.s.            */

  struct unur_ars_interval *next; /* pointer to next interval in list        */

#ifdef DEBUG_STORE_IP 
  double  ip;             /* intersection point between two tangents         */
#endif
#ifdef UNUR_COOKIES
  unsigned cookie;        /* magic cookie                                    */
#endif
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_ars_gen { 

  double  Atotal;               /* area below hat                            */
  double  logAmax;              /* log of maximum area in intervals          */

  struct unur_ars_interval *iv; /* pointer to linked list of intervals       */
  int     n_ivs;                /* number of intervals                       */
  int     max_ivs;              /* maximum number of intervals               */
  int     max_iter;             /* maximum number of iterations              */

  double *starting_cpoints;     /* pointer to array of starting points       */
  int     n_starting_cpoints;   /* number of construction points at start    */

  double *percentiles;       /* percentiles of hat for c. points of new hat  */
  int n_percentiles;         /* number of percentiles                        */
  int retry_ncpoints;        /* number of cpoints for second trial of reinit */
};

/*---------------------------------------------------------------------------*/
