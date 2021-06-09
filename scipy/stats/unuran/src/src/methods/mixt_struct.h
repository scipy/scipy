/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: mixt_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for meta method MIXT                          *
 *         (MIXTure of distributions)                                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2010 Wolfgang Hoermann and Josef Leydold                  *
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

struct unur_mixt_par { 
  int n_comp;                   /* number of components                      */
  const double *prob;           /* probabilities (weights) for components    */
  struct unur_gen **comp;       /* array of pointers to components           */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_mixt_gen { 
  int is_inversion;             /* whether inversion is used                 */

  /* components are stored in slot 'gen_aux_list'                            */
  /* probabilities are stored in slot 'gen_aux' as generator with method DGT */
};

/*---------------------------------------------------------------------------*/
