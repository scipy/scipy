/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: deprecated_vmt_struct.h                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method VMT                                *
 *         (Vector Matrix Transformation)                                    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  THIS METHOD AND THE CORRESPONDING ROUTINES SHOULD NOT BE USED ANY MORE!  *
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

struct unur_vmt_par { 
  int dummy;
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_vmt_gen { 
  struct unur_gen **marginalgen_list;   /* list of generators for marginal distributions */
  double *cholesky;         /* cholesky factor of covariance matrix          */
  int    dim;               /* dimension of distribution                     */
};

/*---------------------------------------------------------------------------*/

