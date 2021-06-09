/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: hist_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method HIST                               *
 *         (HISTogram of empirical distribution)                             *
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

struct unur_hist_par {
  /* the histogram is stored in the distribution object */
  int dummy;
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_hist_gen {
  int     n_hist;       /* number of bins in histogram                       */
  double *prob;         /* probabilities for bins                            */
  double *bins;         /* location of bins (when different width)           */
  double  hmin, hmax;   /* lower and upper bound for histogram               */
  double  hwidth;       /* width of bins (when equal width)                  */
  double  sum;          /* sum of all probabilities = cumpv[len-1]           */
  double *cumpv;        /* pointer to the vector of cumulated probabilities  */
  int    *guide_table;  /* pointer to guide table                            */
};

/*---------------------------------------------------------------------------*/


