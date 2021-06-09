/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: lobatto_struct.h                                                  *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for Gauss-Lobatto integration                 *
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
/* 'Lobatto object':                                                         */
/*   store integrand, boundaries and integrals of subintervals computed      */
/*   during adaptive Gauss-Lobatto integration.                              */

struct unur_lobatto_nodes {
  double x;   /* right boundary of subinterval */
  double u;   /* integral of PDF over subinterval */
}; 

struct unur_lobatto_table {
  struct unur_lobatto_nodes *values; /* boundaries and integral values       */
  int n_values;              /* number of stored integral values (nodes)     */
  int cur_iv;                /* position of first entry whose x value is
				larger than left boundary of current interval*/
  int size;                  /* size of table                                */
  
  UNUR_LOBATTO_FUNCT *funct; /* pointer to integrand                         */
  struct unur_gen *gen;      /* pointer to generator object                  */
  double tol;                /* tolerated ABSOLUTE integration error         */
  UNUR_LOBATTO_ERROR *uerror; /* function for estimating error               */
  double bleft;              /* left boundary of computational domain        */
  double bright;             /* right boundary of computational domain       */
  double integral;           /* integral over whole domain                   */
}; 


/*---------------------------------------------------------------------------*/
