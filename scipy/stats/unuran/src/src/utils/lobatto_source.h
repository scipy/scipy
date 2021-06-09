/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      lobatto_source.c                                             *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         macros and function prototypes for Gauss-Lobatto integration      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2009 Wolfgang Hoermann and Josef Leydold                  *
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
/* typedefs                                                                  */

/* Integrand */
typedef double UNUR_LOBATTO_FUNCT(double x, struct unur_gen *gen);

/* Error estimate */
typedef double UNUR_LOBATTO_ERROR(struct unur_gen *gen, double delta, double x);

struct unur_lobatto_table;

/*---------------------------------------------------------------------------*/
/* Function prototypes                                                       */

double _unur_lobatto_adaptive (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen,
			       double x, double h, double tol, UNUR_LOBATTO_ERROR uerror);
/*---------------------------------------------------------------------------*/
/* numerical integration of 'funct' over the interval (x,x+h) using          */
/* adaptive Gauss-Lobatto integration with 5 points for each recursion.      */
/*---------------------------------------------------------------------------*/

struct unur_lobatto_table *
_unur_lobatto_init (UNUR_LOBATTO_FUNCT funct, struct unur_gen *gen,
		    double left, double center, double right,
		    double tol, UNUR_LOBATTO_ERROR uerror, int size);
/*---------------------------------------------------------------------------*/
/* create object for Lobatto integral.                                       */
/*---------------------------------------------------------------------------*/

int _unur_lobatto_find_linear (struct unur_lobatto_table *Itable, double x);
/*---------------------------------------------------------------------------*/
/* find first subinterval where left boundary is not less than x.            */ 
/*---------------------------------------------------------------------------*/

double _unur_lobatto_eval_diff (struct unur_lobatto_table *Itable, double x, double h, double *fx);
/*---------------------------------------------------------------------------*/
/* evaluate integration object over the interval (x,x+h).                    */
/*---------------------------------------------------------------------------*/

double _unur_lobatto_eval_CDF (struct unur_lobatto_table *Itable, double x);
/*---------------------------------------------------------------------------*/
/* evaluate integration object over the interval (-INFINITY, x).             */
/* it is important, that the integration object 'Itable' already exists.     */
/*---------------------------------------------------------------------------*/

double _unur_lobatto_integral (struct unur_lobatto_table *Itable );
/*---------------------------------------------------------------------------*/
/* get value of integral from Lobatto object.                                */
/*---------------------------------------------------------------------------*/

void _unur_lobatto_free (struct unur_lobatto_table **Itable);
/*---------------------------------------------------------------------------*/
/* destroy Lobatto object and set pointer to NULL.                           */
/*---------------------------------------------------------------------------*/

void _unur_lobatto_debug_table (struct unur_lobatto_table *Itable,
				const struct unur_gen *gen, int print_Itable );
/*---------------------------------------------------------------------------*/
/* print size and entries of table of integral values.                       */
/*---------------------------------------------------------------------------*/

int _unur_lobatto_size_table (struct unur_lobatto_table *Itable);
/*---------------------------------------------------------------------------*/
/* size of table of integral values.                                         */
/*---------------------------------------------------------------------------*/

