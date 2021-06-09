/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      correlation.c                                                *
 *                                                                           *
 *   compute correlation coefficient of two samples                          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   REFERENCES:                                                             *
 *   [1] Spicer C.C. (1972): Algorithm AS 52: Calculation of Power Sums of   *
 *       Deviations about the mean, Applied Statistics 21(2), pp. 226-227.   *
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

#include <unur_source.h>
#include <distr/cont.h>
#include <utils/matrix_source.h>
#include <methods/unur_methods_source.h>
#include <methods/x_gen_source.h>
#include "unuran_tests.h"


/* #include <limits.h> */
/* #include <distr/cvec.h> */
/* #include <distr/discr.h> */
/* #include <distr/distr_source.h> */
/* #include <distributions/unur_distributions.h> */
/* #include <specfunct/unur_specfunct_source.h> */
/* #include <utils/matrix_source.h> */

/*---------------------------------------------------------------------------*/
static char test_name[] = "Correlation";
/*---------------------------------------------------------------------------*/

#define CORR_DEFAULT_SAMPLESIZE 10000
/* default sample size used when the given size is <= 0 */

#define CORR_MAX_SAMPLESIZE 10000000
/* maximal sample size to prevent extremely long run times */

/*---------------------------------------------------------------------------*/

double
unur_test_correlation( UNUR_GEN *genx, UNUR_GEN *geny, int samplesize, int verbosity, FILE *out )
     /*----------------------------------------------------------------------*/
     /*  compute correlation coefficient of two samples.                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   genx       ... pointer to generator object                         */
     /*   geny       ... pointer to another generator object                 */
     /*   samplesize ... sample size                                         */
     /*   verbosity  ... verbosity level, 0 = no output, 1 = output          */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*     correlation coefficient                                          */
     /*                                                                      */
     /* error:                                                               */
     /*   -2. ... missing data                                               */
     /*   -3. ... other errors                                               */
     /*----------------------------------------------------------------------*/
{
  double x  =0., y =0.;  /* contains random numbers           */ 
  double mx =0., my=0.;  /* mean values (analog to s1 in [1]  */
  double dx =0., dy=0.;  /* analog to delta in [1]            */
  double sx =0., sy=0.;  /* analog to s2 in [1]               */
  double sxy=0.;         /* analog to s2 in [1]               */
  double factor;
  int n;

  /* check parameter */
  _unur_check_NULL(test_name,genx,-3.);
  _unur_check_NULL(test_name,geny,-3.);

  /* type of distribution */
  if (! ( ((genx->method & UNUR_MASK_TYPE) == UNUR_METH_DISCR) ||
	  ((genx->method & UNUR_MASK_TYPE) == UNUR_METH_CONT) )) {
    _unur_error(test_name,UNUR_ERR_GENERIC,
         "dont know how to compute correlation coefficient for distribution");
    return -2.;
  }
  if (! ( ((geny->method & UNUR_MASK_TYPE) == UNUR_METH_DISCR) ||
	  ((geny->method & UNUR_MASK_TYPE) == UNUR_METH_CONT) )) {
    _unur_error(test_name,UNUR_ERR_GENERIC,
         "dont know how to compute correlation coefficient for distribution");
    return -2.;
  }

  /* samplesize */
  if( samplesize <= 0 ) samplesize = CORR_DEFAULT_SAMPLESIZE;
  samplesize = _unur_min( samplesize, CORR_MAX_SAMPLESIZE );

  /* sampling */  
  for (n=1; n<=samplesize; n++) {

    /* which type of distribution */
    switch (genx->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      x = _unur_sample_discr(genx); break;
    case UNUR_METH_CONT:
      x = _unur_sample_cont(genx); break;
    }
    switch (geny->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      y = _unur_sample_discr(geny); break;
    case UNUR_METH_CONT:
      y = _unur_sample_cont(geny); break;
    }

    factor = (double) ( n*(n-1) );

    dx = (x - mx) / n;
    dy = (y - my) / n;
    mx += dx;
    my += dy;

    sx  += factor * dx*dx;
    sy  += factor * dy*dy;
    sxy += factor * dx*dy;
    /* fprintf(out,"sx:%f sy:%f sxy:%f, c:%f\n",sx,sy,sxy, sxy/(sqrt(sx*sy)) ); */

  }

  /* now print results */
  if (verbosity) {
    fprintf(out,"\nCorrelation coefficient: %g\n\n", sxy/(sqrt(sx*sy)) );
  }

  return ( sxy/(sqrt(sx*sy)) );

} /* end of unur_test_correlation() */

/*---------------------------------------------------------------------------*/

int
unur_test_cvec_rankcorr( double *rc, struct unur_gen *gen, int samplesize, int verbose, FILE *out )
     /*----------------------------------------------------------------------*/
     /*  compute rank correlations between components of random vector.      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   rc         ... pointer to dimxdim array to result                  */
     /*   gen        ... pointer to multivariate generator object            */
     /*   samplesize ... sample size                                         */
     /*   verbose    ... verbosity level, 0 = no output, 1 = output          */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
#define DISTR   gen->distr->data.cvec
#define idx(a,b) ((a)*dim+(b))

  double *mx;   /* mean values (analog to s1 in [1]  */
  double *dx;   /* analog to delta in [1]            */
  double factor;

  int dim;         /* dimension of multivariate distribution */

  UNUR_DISTR **marginals;  /* pointer to marginal distributions */
  UNUR_FUNCT_CONT **marginal_cdf;  /* pointer to CDFs of marginal distributions */

  double *X;             /* sampling vector */
  double *U;             /* X transformed to uniform */

  int i, j, n;     /* auxiliary variables */

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  /* we do not check magic cookies here */

  if (verbose >= 1)
    fprintf(out,"\nRank correlations of random vector:\n");

  /* samplesize */
  if( samplesize <= 0 ) samplesize = CORR_DEFAULT_SAMPLESIZE;
  samplesize = _unur_min( samplesize, CORR_MAX_SAMPLESIZE );

  /* dimension of distribution */
  dim = gen->distr->dim;
  if (dim < 1) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"distribution dimension < 1 ?");
    return UNUR_ERR_GENERIC;
  }

  /* type of distribution */
  if ( (gen->method & UNUR_MASK_TYPE) != UNUR_METH_VEC ) {
    _unur_error(test_name,UNUR_ERR_GENERIC,
		"rank correlation coefficients cannot be computed");
    return UNUR_ERR_GENERIC;
  }

  /* we need all marginal distributions */
  if (DISTR.marginals==NULL) {
    _unur_error(gen->distr->name,UNUR_ERR_DISTR_REQUIRED,"marginal distributions");
    return UNUR_ERR_DISTR_REQUIRED; }
  marginals = _unur_xmalloc(dim * sizeof(UNUR_DISTR *));
  marginal_cdf = _unur_xmalloc(dim * sizeof(UNUR_FUNCT_CONT *));
  for (i=0; i<dim; i++) {
    marginals[i] = DISTR.marginals[i];
    marginal_cdf[i] = unur_distr_cont_get_cdf(DISTR.marginals[i]);
    if (marginals[i]==NULL || marginal_cdf[i]==NULL) {
      _unur_error(gen->distr->name,UNUR_ERR_DISTR_REQUIRED,"CDF of continuous marginal");
      free (marginals);  free (marginal_cdf);
      return UNUR_ERR_DISTR_REQUIRED; }
  }

  /* allocate working space memory */
  X   = _unur_xmalloc( dim * sizeof(double));
  U   = _unur_xmalloc( dim * sizeof(double));
  mx  = _unur_xmalloc( dim * sizeof(double));
  dx  = _unur_xmalloc( dim * sizeof(double));

  /* clear working arrays */
  for (i=0; i<dim; i++)
    mx[i] = dx[i] = 0.;
  for (i=0; i<dim*dim; i++) 
    rc[i] = 0.;

  /* now run generator */
  for (n=1; n<=samplesize; n++) {

    factor = ((double)n)*(n-1.);

    /* get random vector X */
    _unur_sample_vec(gen, X);

    for (i=0; i<dim; i++) {
      /* transform to uniform U */
      U[i] = marginal_cdf[i](X[i],marginals[i]);

      dx[i] = (U[i] - mx[i]) / n;
      mx[i] += dx[i];
    }

    for (i=0; i<dim; i++)
      for (j=i; j<dim; j++)
	rc[idx(i,j)] += factor * dx[i] * dx[j];
  }

  /* compute rank correlation */
  for (i=0; i<dim; i++) {
    for (j=0; j<i; j++)
      rc[idx(i,j)] = rc[idx(j,i)];
    for (j=i+1; j<dim; j++)
      rc[idx(i,j)] /= sqrt(rc[idx(i,i)] * rc[idx(j,j)]);
    rc[idx(i,i)] = 1.;
  }

  /* print result */
  if (verbose >= 1)
    _unur_matrix_print_matrix( dim, rc, "rank correlation =", 
			       out, "", "\t");

  /* free memory */
  if (X)    free(X);
  if (U)    free(U);
  if (mx)   free(mx);
  if (dx)   free(dx);
  if (marginals)  free (marginals);
  if (marginal_cdf)  free (marginal_cdf);

  /* return result of test */
  return UNUR_SUCCESS;

#undef DISTR
#undef idx
} /* end of unur_test_cvec_rankcorr() */

/*---------------------------------------------------------------------------*/
