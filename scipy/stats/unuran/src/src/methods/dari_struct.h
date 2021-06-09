/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: dari_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method DARI                               *
 *         ((Discrete) Alias-Urn)                                            *
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

struct unur_dari_par { 
  int     squeeze;       /* should the squeeze be used  
                            0.. no squeeze,  1..squeeze                      */
  int     size;          /* size of table for speeding up generation         */
  double  c_factor;      /* constant for choosing the design points          */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_dari_gen { 
  double  vt;            /* total volume below hat                           */
  double  vc;            /* volume below center part                         */
  double  vcr;           /* volume center and right together                 */

  double  xsq[2];        /* value necesary for the squeeze computation       */
  double  y[2];          /* value of the transformed density in points of contact */
  double  ys[2];         /* the slope of the transformed hat                 */
  double  ac[2];         /* left and right starting point of the uniform hat 
                            in the center                                    */

  double  pm;            /* mode probability                                 */
  double  Hat[2];        /* point where the hat starts for the left and
                            the right tail                                   */
  double  c_factor;      /* constant for choosing the design points          */

  int     m;             /* mode                                             */
  int     x[2];          /* points of contact left and right of the mode     */
  int     s[2];          /* first and last integer of the center part        */
  int     n[2];          /* contains the first and the last i 
                            for which values are stored in table             */
  int     size;          /* size of the auxiliary tables                     */
  int     squeeze;       /* use squeeze yes/no                               */

  double *hp;            /* pointer to double array of length size           */
  char   *hb;            /* pointer to boolean array of length size          */
};

/*---------------------------------------------------------------------------*/


