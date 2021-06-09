/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: vector.c                                                          *
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

#include <unur_source.h>

/*---------------------------------------------------------------------------*/

/* allocate memory for new vector */
double *
_unur_vector_new(int dim)
{
  int i;
  double *v;

  v = _unur_xmalloc(dim*sizeof(double));

  /* setting all coordinates to 0 */
  for (i=0; i<dim; i++) v[i] = 0.;

  return v;
} /* end of _unur_vector_new() */

/*--------------------------------------------------------------------------*/

/* free allocated memory used by vector structure */
void 
_unur_vector_free(double *v)
{
  if (v) free(v);
} /* end of _unur_vector_free() */

/*--------------------------------------------------------------------------*/

/* calculation of euclidean (L2) norm of vector */
/* avoid overflow by using                      */
/*   ||v|| = v_m * sqrt(sum (v_i/v_m)^2)        */
/* where v_m = max |v_i|                        */
double 
_unur_vector_norm(int dim, double *v)
{
  int i;
  double vsum;
  double vmax;
  double p;

  /* checking if v is NULL */
  if (v==NULL) return 0.; 

  /* determining the largest element (absolute values) */
  vmax = 0.;
  for (i=0; i<dim; i++) {
    if (vmax < fabs(v[i])) vmax = fabs(v[i]); 
  }
  
  /* case: null vector */
  if (vmax <= 0) return 0.;

  /* it's nummerically more stable to calculate the norm this way */
  vsum = 0.;
  for (i=0; i<dim; i++) {
    p=v[i]/vmax;
    vsum += p*p;
  }

  /* return L2 norm */
  return vmax * sqrt(vsum);

} /* end of _unur_vector_norm() */

/*--------------------------------------------------------------------------*/

/* normalize a vector to have unit norm */
void 
_unur_vector_normalize(int dim, double *v)
{
  int i;
  double norm;
  
  /* checking if v is NULL */
  if (v==NULL) return; 
  
  norm = _unur_vector_norm(dim, v);
  
  for (i=0; i<dim; i++)  v[i] /= norm;

} /* end of _unur_vector_normalize() */

/*--------------------------------------------------------------------------*/

/* calculation of scalar product */
double 
_unur_vector_scalar_product(int dim, double *v1, double *v2)
{
  int i;
  double scalar_product;
  
  /* checking if v1 or v2 are NULL */
  if (v1==NULL || v2==NULL) return 0.; 
  
  scalar_product = 0.;
  for (i=0; i<dim; i++) {
    scalar_product += v1[i]*v2[i];
  }

  return scalar_product;
} /* end of _unur_vector_scalar_product() */

/*--------------------------------------------------------------------------*/

