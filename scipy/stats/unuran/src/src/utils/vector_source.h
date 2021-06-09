/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: vector_source.h                                                   *
 *                                                                           *
 *   Routines for computations with vectors (arrays).                        *
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
				                                                                                    
/*--------------------------------------------------------------------------*/

double *_unur_vector_new(int dim);
/* allocate memory for new vector and initialize to 0 */

void _unur_vector_free(double *v);
/* free allocated memory used by vector */

double _unur_vector_norm(int dim, double *v);
/* calculation of euclidean (L2) norm of vector */

double _unur_vector_scalar_product(int dim, double *v1, double *v2);
/* calculation of scalar product */

void _unur_vector_normalize(int dim, double *v);
/* normalize a vector to have unit norm */

/*---------------------------------------------------------------------------*/
