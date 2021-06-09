/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      moments.c                                                    *
 *                                                                           *
 *   compute central moments of samples                                      *
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
 *****************************************************************************
 *                                                                           *
 *   REFERENCES:                                                             *
 *   [1] Spicer C.C. (1972): Algorithm AS 52: Calculation of Power Sums of   *
 *       Deviations about the mean, Applied Statistics 21(2), pp. 226-227.   *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <methods/unur_methods_source.h>
#include <methods/x_gen_source.h>
#include "unuran_tests.h"

/*---------------------------------------------------------------------------*/
static char test_name[] = "Moments";
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

int
unur_test_moments( UNUR_GEN *gen, double *moments, int n_moments, int samplesize,
		   int verbosity, FILE *out )
     /*----------------------------------------------------------------------*/
     /*  compute central moments of samples.                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   moments    ... array for storing moments                           */
     /*                  for multivariate distributions the moments are      */
     /*                  stored consecutively for each dimension.            */
     /*                  arraylength = (n_moments+1) * dim                   */ 
     /*   n_moments  ... number of moments to be calculated (at most 4)      */
     /*   samplesize ... sample size                                         */
     /*   verbosity  ... verbosity level, 0 = no output, 1 = output          */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
#define idx(d,n) ((d)*(n_moments+1)+(n))

  double an, an1, dx, dx2;
  int n, mom, d, dim;
  double *x;

  /* check parameter */
  _unur_check_NULL(test_name, gen, UNUR_ERR_NULL);

  /* type of distribution */
  if (! ( ((gen->method & UNUR_MASK_TYPE) == UNUR_METH_DISCR) ||
	  ((gen->method & UNUR_MASK_TYPE) == UNUR_METH_CONT)  ||
	  ((gen->method & UNUR_MASK_TYPE) == UNUR_METH_VEC) )) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"dont know how to compute moments for distribution");
    return UNUR_ERR_GENERIC;
  }
  /* array for storing moments */
  CHECK_NULL(moments,0);
  if (n_moments <= 0 || n_moments > 4) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"number of moments < 1 or > 4");
    return UNUR_ERR_GENERIC;
  }

  /* sample size >= 10 */
  if (samplesize < 10) 
    samplesize = 10;

  /* number of dimensions can only be > 1 for multivariate case  */  
  dim = 1;
  if ((gen->method & UNUR_MASK_TYPE) == UNUR_METH_VEC) dim = gen->distr->dim;
  
  /* allocating memory for sampling "vector" */
  x = _unur_xmalloc(dim * sizeof(double));
    
  /* clear array of moments */
  for (d=0; d<dim; d++) {
    moments[idx(d,0)] = 1.; /* dummy field */
    for (mom = 1; mom <= n_moments; mom++ )
      moments[idx(d,mom)] = 0.;
  }
   
  /* sampling */
  /* compute moments: we use a recurrence relation by Spicer [1]. */
  
  for (n=1; n<=samplesize; n++) {

    /* which type of distribution */
    switch (gen->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      x[0] = _unur_sample_discr(gen); break;
    case UNUR_METH_CONT:
      x[0] = _unur_sample_cont(gen); break;
    case UNUR_METH_VEC:
      _unur_sample_vec(gen, x); break;
    }

    for (d=0; d<dim; d++) {
      an = (double)n;
      an1 = an-1.;
      dx = (x[d] - moments[idx(d,1)]) / an;
      dx2 = dx * dx;
   
      switch (n_moments) {
      case 4:
        moments[idx(d,4)] -= dx * (4.*moments[idx(d,3)] - dx * (6.*moments[idx(d,2)] + an1*(1. + an1*an1*an1)*dx2));
      case 3:
        moments[idx(d,3)] -= dx * (3.*moments[idx(d,2)] - an*an1*(an-2.)*dx2);
      case 2:
        moments[idx(d,2)] += an * an1 * dx2;
      case 1:
        moments[idx(d,1)] += dx;
      }
    }
  }

  /* compute moments */
  /* moments[1] is already the first moment. no division necessary. */
  for (d=0; d<dim; d++) {
    for (mom = 2; mom <= n_moments; mom++ )
      moments[idx(d,mom)] /= samplesize;
  
    /* now print results */
    if (verbosity) {
      if (dim==1) fprintf(out,"\nCentral MOMENTS:\n");
      else        fprintf(out,"\nCentral MOMENTS for dimension #%d:\n", d);
      for (mom = 1; mom <= n_moments; mom++ )
        fprintf(out,"\t[%d] =\t%g\n",mom,moments[idx(d,mom)]);
      fprintf(out,"\n");
    }
  }

  /* release memory */
  free(x);
  
  return UNUR_SUCCESS;

#undef idx  
} /* end of unur_test_moments() */

/*---------------------------------------------------------------------------*/

