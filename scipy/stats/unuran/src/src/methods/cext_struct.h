/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: cext_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method CEXT                               *
 *         (wrapper for Continuous EXTernal generators)                      *
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

struct unur_cext_par {   /* data for external generator */
  int    (*init)  (UNUR_GEN *gen);   /* pointer to initialization routine    */
  double (*sample)(UNUR_GEN *gen);   /* pointer to sampling routine          */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_cext_gen { 
  int    (*init)  (UNUR_GEN *gen);   /* pointer to initialization routine    */
  double (*sample)(UNUR_GEN *gen);   /* pointer to sampling routine          */

  void *param;               /* parameters for the generator                 */
  size_t size_param;         /* size of parameter object                     */
};

/*---------------------------------------------------------------------------*/
