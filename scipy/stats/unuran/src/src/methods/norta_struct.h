/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: norta_struct.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method NORTA                              *
 *         (NORmal To Anything)                                              *
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

struct unur_norta_par { 
  int dummy;
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_norta_gen { 
  int    dim;                          /* dimension of distribution          */
  double *copula;                      /* pointer to intermediate copula     */
  struct unur_distr *normaldistr;      /* standard normal distribution       */
  struct unur_gen **marginalgen_list;  /* list of generators for marginal distributions */

  /* Remark: We use gen->gen_aux to store the pointer to the                 */
  /*         multinormal generator.                                          */
  /*         It is accessed via the macro 'MNORMAL'.                         */
};

/*---------------------------------------------------------------------------*/

