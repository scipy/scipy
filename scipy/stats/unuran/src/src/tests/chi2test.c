/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      chi2test.c                                                   *
 *                                                                           *
 *   Chi^2 Tests (to validate implementation)                                *
 *                                                                           *
 *   WARNING!                                                                *
 *   A succesfull chi^2 test shows that the implementation is (probably)     *
 *   correct. It DOES NOT mean that the generator is of good quality!        *
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

#include <limits.h>
#include <unur_source.h>
#include <methods/unur_methods_source.h>
#include <methods/x_gen_source.h>
#include <distr/cont.h>
#include <distr/cvec.h>
#include <distr/discr.h>
#include <distr/distr_source.h>
#include <distributions/unur_distributions.h>
#include <specfunct/unur_specfunct_source.h>
#include <utils/matrix_source.h>
#include "unuran_tests.h"

/*---------------------------------------------------------------------------*/
/* constants */

#define CHI2_SAMPLEFAC  40
/* if samplesize<=0 use samplesize = CHI2_SAMPLEFAC * intervals^dim */

#define CHI2_CLASSMIN_DEFAULT  20
/* default number of observations in class */

#define CHI2_INTERVALS_DEFAULT 50
/* default number of intervals for chi^2 test if given number is <= 0 */

#define CHI2_DEFAULT_SAMPLESIZE 10000
/* default sample size used when the given size is <= 0 */

#define CHI2_MAX_SAMPLESIZE 1000000
/* maximal sample size to prevent extremely long run times */

#define CHI2_MAX_TOTALINTERVALS 1000000
/* maximal product of intervals used in chi2vec test */

/*---------------------------------------------------------------------------*/
static char test_name[] = "Chi^2-Test";
/*---------------------------------------------------------------------------*/

static double _unur_test_chi2_discr( struct unur_gen *gen, int samplesize, int classmin,
				     int verbose, FILE *out );

static double _unur_test_chi2_cont( struct unur_gen *gen, int n_intervals, int samplesize, int classmin,
				    int verbose, FILE *out );

static double _unur_test_chi2_cemp( struct unur_gen *gen, int n_intervals, int samplesize, int classmin,
				    int verbose, FILE *out );

static double _unur_test_chi2_vec( struct unur_gen *gen, int n_intervals, int samplesize, int classmin,
				   int verbose, FILE *out );

static double _unur_test_chi2_cvemp( struct unur_gen *gen, int n_intervals, int samplesize, int classmin,
				     int verbose, FILE *out );

static double _unur_test_chi2test( double *prob, int *observed, int len, int classmin,
				   int verbose, FILE *out );

/*---------------------------------------------------------------------------*/

double
unur_test_chi2( struct unur_gen *gen,
    int n_intervals,
    int samplesize,
    int classmin,
    int verbose,
    FILE *out )
     /*----------------------------------------------------------------------*/
     /* Chi^2 test for univariate distributions                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator                                */
     /*   n_intervals .. number if intervals                                 */
     /*                  case probability vector: use its length             */
     /*                  otherwise and if <= 0 use default                   */
     /*   samplesize ... samplesize for test                                 */
     /*                  (if <= 0, (#classes)^2 is used as default.)         */
     /*   classmin   ... minimum number of expected occurrences for each class */
     /*                  (if <= 0, a default value is used.)                 */
     /*   verbose    ... verbosity level                                     */
     /*                  0 = no output on out                                */
     /*                  1 = print summary                                   */
     /*                  2 = print classes and summary                       */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   p-value of test statistics under H_0                               */
     /*                                                                      */
     /* error:                                                               */
     /*   -2. ... missing data                                               */
     /*   -1. ... other errors                                               */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL(test_name,gen,-1.);

  if (verbose >= 1)
    fprintf(out,"\nGOODNESS-OF-FIT TESTS:\n");

  switch (gen->method & UNUR_MASK_TYPE) {

  case UNUR_METH_DISCR:
    return _unur_test_chi2_discr(gen, samplesize, classmin, verbose, out);

  case UNUR_METH_CONT:
    return _unur_test_chi2_cont(gen, n_intervals, samplesize, classmin, verbose, out);

  case UNUR_METH_CEMP:
    return _unur_test_chi2_cemp(gen, n_intervals, samplesize, classmin, verbose, out);

  case UNUR_METH_VEC:
    return _unur_test_chi2_vec(gen, n_intervals, samplesize, classmin, verbose, out);

  case UNUR_METH_CVEMP:
    return _unur_test_chi2_cvemp(gen, n_intervals, samplesize, classmin, verbose, out);

  default:
    _unur_error(test_name,UNUR_ERR_GENERIC,"Not implemented for such distributions!");
    return -1.;
  }
} /* end of unur_test_chi2() */

/*---------------------------------------------------------------------------*/

double
_unur_test_chi2_discr( struct unur_gen *gen,
           int samplesize,
           int classmin,
           int verbose,
           FILE *out)
     /*----------------------------------------------------------------------*/
     /* Chi^2 test for univariate discrete distributions                     */
     /* with given probability vector.                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   samplesize ... samplesize for test                                 */
     /*                  (if <= 0, (#classes)^2 is used as default.)         */
     /*   classmin   ... minimum number of expected occurrences for each class */
     /*                  (if <= 0, a default value is used.)                 */
     /*   verbose    ... verbosity level                                     */
     /*                  0 = no output                                       */
     /*                  1 = print summary                                   */
     /*                  2 = print classes and summary                       */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   p-value of test statistics under H_0                               */
     /*                                                                      */
     /* error:                                                               */
     /*   -2. ... missing data                                               */
     /*   -1. ... other errors                                               */
     /*----------------------------------------------------------------------*/
{
#define DISTR   gen->distr->data.discr
  double *pv;           /* pointer to probability vectors */
  int n_pv;             /* length of probability vector   */

  int *observed;        /* vector for observed occurrences */
  double pval;          /* p-value */
  int had_PV;           /* whether we had a PV for test or not */
  int i,j;

  /* check arguments */
  CHECK_NULL(gen,-1.);

  /* probability vector */
  if (DISTR.pv == NULL) {
    had_PV = FALSE;
    /* no PV given --> try to compute PV */
    if (!unur_distr_discr_make_pv( gen->distr )) {
      /* not successful */
      return -2.;
    }
  }
  else
    had_PV = TRUE;
  /* pointer to PV */
  pv = DISTR.pv;
  n_pv = DISTR.n_pv;

  /* allocate memory for observations */
  observed = _unur_xmalloc( n_pv * sizeof(int));

  /* clear array */
  for( i=0; i<n_pv; i++ )
    observed[i] = 0;

  /* samplesize */
  if( samplesize <= 0 ) {
    samplesize = (INT_MAX/n_pv > n_pv) ? n_pv * n_pv : 1000000;
    samplesize = _unur_max(samplesize,1000000);
  }

  samplesize = _unur_min( samplesize, CHI2_MAX_SAMPLESIZE );

  /* now run generator */
  for( i=0; i<samplesize; i++ ) {
    /* sample */
    j = _unur_sample_discr(gen);
    if (verbose >= 3) fprintf(out,"i = %d\n",j);
    /* shift vector */
    j -= DISTR.domain[0];
    /* check range of random variates !! */
    if (j >= 0 && j < n_pv) 
      ++observed[j];
    /* else: ignore number --> chop off tail */
  }

  if (verbose >= 1) {
    fprintf(out,"\nChi^2-Test for discrete distribution with given probability vector:");
    fprintf(out,"\n  length     = %d\n",n_pv);
  }

  /* and now make chi^2 test */
  pval = _unur_test_chi2test(pv, observed, n_pv, classmin, verbose, out);

  /* free memory */
  free(observed);
  if (!had_PV) {
    free (DISTR.pv);
    DISTR.pv = NULL;
    DISTR.n_pv = 0;
  }

  /* return result of test */
  return pval;

#undef DISTR
} /* end of _unur_test_chi2_discr() */

/*---------------------------------------------------------------------------*/

double
_unur_test_chi2_cont( struct unur_gen *gen,
          int n_intervals,
          int samplesize,
          int classmin,
          int verbose,
          FILE *out )
     /*----------------------------------------------------------------------*/
     /* Chi^2 test for univariate continuous distributions.                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   n_intervals .. number of intervals in which (0,1) is partitioned   */
     /*   samplesize ... samplesize for test                                 */
     /*                  (if <= 0, intervals^2 is used as default.)          */
     /*   classmin   ... minimum number of expected occurrences for each class */
     /*                  (if <= 0, a default value is used.)                 */
     /*   verbose    ... verbosity level                                     */
     /*                  0 = no output                                       */
     /*                  1 = print summary                                   */
     /*                  2 = print classes and summary                       */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   p-value of test statistics under H_0                               */
     /*                                                                      */
     /* error:                                                               */
     /*   -2. ... missing data                                               */
     /*   -1. ... other errors                                               */
     /*----------------------------------------------------------------------*/
{
#define DISTR   gen->distr->data.cont

  double F, Fl, Fr, Fdelta;  /* value of CDF (at left and right boundary point) */
  UNUR_FUNCT_CONT *cdf;      /* pointer to CDF */
  int *observed;             /* vector for observed occurrences */
  double pval;               /* p-value */
  int i,j;

  /* check arguments */
  CHECK_NULL(gen,-1.);
  /* we do not check magic cookies here */

  /* CDF required */
  cdf = DISTR.cdf;
  if (DISTR.cdf == NULL) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"CDF required for continuous random variates!");
    return -2.;
  }

  /* check given number of intervals */
  if (n_intervals <= 2)
    n_intervals = CHI2_INTERVALS_DEFAULT;

  /* allocate memory for observations */
  observed = _unur_xmalloc( n_intervals * sizeof(int));

  /* clear array */
  for( i=0; i<n_intervals; i++ )
    observed[i] = 0;

  /* samplesize */
  if( samplesize <= 0 )
    samplesize = (INT_MAX/n_intervals > n_intervals) ? n_intervals*n_intervals : INT_MAX;

  samplesize = _unur_min( samplesize, CHI2_MAX_SAMPLESIZE );

  /* compute Fl and Fr */
  if (gen->distr->set & UNUR_DISTR_SET_TRUNCATED) {
    Fl = (DISTR.trunc[0] <= -INFINITY) ? 0. : cdf(DISTR.trunc[0], gen->distr);
    Fr = (DISTR.trunc[1] >=  INFINITY) ? 1. : cdf(DISTR.trunc[1], gen->distr);
  }
  else {
    Fl = (DISTR.domain[0] <= -INFINITY) ? 0. : cdf(DISTR.domain[0], gen->distr);
    Fr = (DISTR.domain[1] >=  INFINITY) ? 1. : cdf(DISTR.domain[1], gen->distr);
  }
  Fdelta = Fr - Fl;

  /* Fr - Fl <= 0. is a fatal error */
  if (Fdelta <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_GENERIC,"Fdelta <= 0.");
    free (observed);
    return -1.;
  }

  /* now run generator */
  for( i=0; i<samplesize; i++ ) {
    if (verbose >= 3) {
      double x = _unur_sample_cont(gen);
      F = cdf( x, gen->distr );
      fprintf(out,"x = %g\n",x);
    }
    else {
      F = cdf( _unur_sample_cont(gen), gen->distr );
    }
    F = (F-Fl)/Fdelta;
    j = (int)(n_intervals * F);
    if (j > n_intervals) {
      _unur_warning(test_name,UNUR_ERR_GENERIC,"F(x) > Fmax (out of domain).");
      j = n_intervals-1;
    }
    if (j >= n_intervals)    /* cdf() might return 1. */
      j = n_intervals-1;
    if (j < 0 ) {           /* there is something wrong with the boundaries */
      _unur_warning(test_name,UNUR_ERR_GENERIC,"F(x) < 0 (out of domain).");
      j = 0;
    }
    ++observed[j];
  }

  if (verbose >= 1) {
    fprintf(out,"\nChi^2-Test for continuous distribution:");
    fprintf(out,"\n  intervals  = %d\n",n_intervals);
  }

  /* and now make chi^2 test */
  pval = _unur_test_chi2test(NULL, observed, n_intervals, classmin, verbose, out );

  /* free memory */
  free(observed);

  /* return result of test */
  return pval;

#undef DISTR
} /* end of _unur_test_chi2_cont() */

/*---------------------------------------------------------------------------*/

double
_unur_test_chi2_cemp( struct unur_gen *gen,
          int n_intervals,
          int samplesize,
          int classmin,
          int verbose,
          FILE *out )
     /*----------------------------------------------------------------------*/
     /* Chi^2 test for continuous empirical distributions.                   */
     /* Tests for standard normal distribution only!                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   n_intervals .. number of intervals in which (0,1) is partitioned   */
     /*   samplesize ... samplesize for test                                 */
     /*                  (if <= 0, intervals^2 is used as default.)          */
     /*   classmin   ... minimum number of expected occurrences for each class */
     /*                  (if <= 0, a default value is used.)                 */
     /*   verbose    ... verbosity level                                     */
     /*                  0 = no output                                       */
     /*                  1 = print summary                                   */
     /*                  2 = print classes and summary                       */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   p-value of test statistics under H_0                               */
     /*                                                                      */
     /* error:                                                               */
     /*   -2. ... missing data                                               */
     /*   -1. ... other errors                                               */
     /*----------------------------------------------------------------------*/
{
  UNUR_DISTR *distr_normal;  /* object for standard normal distribution */
  UNUR_FUNCT_CONT *cdf;      /* pointer to CDF */
  double F;                  /* value of CDF */
  int *observed;             /* vector for observed occurrences */
  double pval;               /* p-value */
  int i,j;

  /* check arguments */
  CHECK_NULL(gen,-1.);
  /* we do not check magic cookies here */

  /* CDF for standard normal distribution */
  distr_normal = unur_distr_normal( NULL, 0 );
  cdf = distr_normal->data.cont.cdf;

  /* check given number of intervals */
  if (n_intervals <= 2)
    n_intervals = CHI2_INTERVALS_DEFAULT;

  /* allocate memory for observations */
  observed = _unur_xmalloc( n_intervals * sizeof(int));

  /* clear array */
  for( i=0; i<n_intervals; i++ )
    observed[i] = 0;

  /* samplesize */
  if( samplesize <= 0 )
    samplesize = (INT_MAX/n_intervals > n_intervals) ? n_intervals*n_intervals : INT_MAX;
  samplesize = _unur_min( samplesize, CHI2_MAX_SAMPLESIZE );

  /* now run generator */
  for( i=0; i<samplesize; i++ ) {
    F = cdf( _unur_sample_cont(gen), distr_normal );
    j = (int)(n_intervals * F);
    if (j > n_intervals) {
      _unur_warning(test_name,UNUR_ERR_GENERIC,"F(x) > Fmax (out of domain).");
      j = n_intervals-1;
    }
    if (j >= n_intervals)    /* cdf() might return 1. */
      j = n_intervals-1;
    if (j < 0 ) {           /* there is something wrong with the boundaries */
      _unur_warning(test_name,UNUR_ERR_GENERIC,"F(x) < 0 (out of domain).");
      j = 0;
    }
    ++observed[j];
  }

  if (verbose >= 1) {
    fprintf(out,"\nChi^2-Test for continuous empirical distribution:");
    fprintf(out,"\n(Assumes standard normal distribution!)");
    fprintf(out,"\n  intervals  = %d\n",n_intervals);
  }

  /* and now make chi^2 test */
  pval = _unur_test_chi2test(NULL, observed, n_intervals, classmin, verbose, out );

  /* free memory */
  _unur_distr_free(distr_normal);
  free(observed);

  /* return result of test */
  return pval;

} /* end of _unur_test_chi2_cemp() */

/*---------------------------------------------------------------------------*/

double
_unur_test_chi2_vec ( struct unur_gen *gen,
          int n_intervals,
          int samplesize,
          int classmin,
          int verbose,
          FILE *out )
     /*----------------------------------------------------------------------*/
     /* Chi^2 test for multivariate (NORMAL) continuous distributions.       */
     /* It runs a chi^2 test on all marginal distributions and on the        */
     /* multivariate distributions.                                          */
     /* For the marginal distributions the number of intervals is set to     */
     /* n_intervals; for the global test we use n_intervals for the first    */
     /* dimension, n_intervals/2 for the second, n_intervals/3 for the third */
     /* dimension as long as the total number of cells does not exceed       */
     /* CHI2_MAX_TOTALINTERVALS.                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   n_intervals... number of intervals in which (0,1) is partitioned   */
     /*                  (each dimension is partitioned in n_intervals)      */
     /*   samplesize ... samplesize for test                                 */
     /*                  (if <= 0, CHI2_SAMPLEFAC * intervals^dim is used)   */
     /*   classmin   ... minimum number of expected occurrences for each class */
     /*                  (if <= 0, a default value is used.)                 */
     /*   verbose    ... verbosity level                                     */
     /*                  0 = no output                                       */
     /*                  1 = print summary                                   */
     /*                  2 = print classes and summary                       */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   minimal p-value * number_of_tests of all test statistics under H_0 */
     /*                                                                      */
     /* error:                                                               */
     /*   -2. ... missing data                                               */
     /*   -1. ... other errors                                               */
     /*----------------------------------------------------------------------*/
{
#define DISTR     gen->distr->data.cvec
#define idx(i,j)  ((i)*dim+(j))

  int dim;                       /* dimension of multivariate distribution */
  double *Fl = NULL;             /* value of CDF at left and right boundary point */
  double *Fr = NULL;
  double *Fdelta = NULL;

  UNUR_DISTR **marginals = NULL; /* pointer to marginal distributions */
  UNUR_FUNCT_CONT **marginal_cdf = NULL;  /* pointer to CDFs of marginal distributions */

  const double *mean = NULL;     /* pointer to mean vector */
  const double *L = NULL;        /* pointer to Cholesky factor */
  double *Linv = NULL;           /* pointer to inverse Cholesky factor */
  double Linv_det;               /* determinant of Linv */

  double *X = NULL;              /* sampling vector */
  double *U = NULL;              /* X transformed to uniform */

  int *bm = NULL;                /* array for counting bins for marginals */
  double pval, pval_min;         /* p-value */

  int i, j, k;                   /* auxiliary variables */

  /* check arguments */
  CHECK_NULL(gen,-1.);
  /* we do not check magic cookies here */

  /* we need all marginal distributions */
  if (DISTR.marginals==NULL) {
    _unur_error(gen->distr->name,UNUR_ERR_DISTR_REQUIRED,"marginals");
    return -2.; 
  }

  /* when the domain has been changed by a unur_distr_cvec_set_domain_...() */
  /* call and the covariance matrix is not the identity matrix then         */
  /* the test may fail anyway.                                              */
  if ((gen->distr->set & UNUR_DISTR_SET_DOMAINBOUNDED) &&
      !(gen->distr->set & UNUR_DISTR_SET_COVAR_IDENT) )
    _unur_warning(test_name,UNUR_ERR_GENERIC,"correlated and domain truncated --> test might fail");

  /* check given number of intervals */
  if (n_intervals <= 2)
    n_intervals = CHI2_INTERVALS_DEFAULT;

  /* samplesize */
  if( samplesize <= 0 ) samplesize = CHI2_DEFAULT_SAMPLESIZE;
  samplesize = _unur_min( samplesize, CHI2_MAX_SAMPLESIZE );

  /* dimension of distribution */
  dim = gen->distr->dim;  
  if (dim < 1) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"distribution dimension < 1 ?");
    return -1.;
  }

  /* we need mean vector and covariance matrix */
  mean = unur_distr_cvec_get_mean(gen->distr);
  L = unur_distr_cvec_get_cholesky(gen->distr);
  /* remark: mean or L might be NULL */

  /* get marginal distributions and their CDFs */
  marginals = _unur_xmalloc(dim * sizeof(UNUR_DISTR *));
  marginal_cdf = _unur_xmalloc(dim * sizeof(UNUR_FUNCT_CONT *));
  for (i=0; i<dim; i++) {
    marginals[i] = DISTR.marginals[i];
    marginal_cdf[i] = unur_distr_cont_get_cdf(DISTR.marginals[i]);
    if (marginals[i]==NULL || marginal_cdf[i]==NULL) {
      _unur_error(gen->distr->name,UNUR_ERR_DISTR_REQUIRED,"CDF of continuous standardized marginal");
      pval_min = -2.; goto free_memory;
    }
  }

  /* compute Fl and Fr */
  Fl  = _unur_xmalloc(dim * sizeof(double));
  Fr  = _unur_xmalloc(dim * sizeof(double));
  Fdelta = _unur_xmalloc(dim * sizeof(double));
  if (gen->distr->set & UNUR_DISTR_SET_DOMAINBOUNDED) {
    for (i=0; i<dim; i++) {
      Fl[i] = ( (!_unur_isfinite(DISTR.domainrect[2*i])) ? 0. 
		: marginal_cdf[i](DISTR.domainrect[2*i],marginals[i]) );
      Fr[i] = ( (!_unur_isfinite(DISTR.domainrect[2*i+1])) ? 1. 
		: marginal_cdf[i](DISTR.domainrect[2*i+1],marginals[i]) );
      Fdelta[i] = Fr[i] - Fl[i];
      /* Fr - Fl <= 0. is a fatal error */
      if (Fdelta[i] <= 0.) {
	_unur_error(gen->genid,UNUR_ERR_GENERIC,"Fdelta <= 0.");
	pval_min = -1.; goto free_memory;
      }
    }
  }
  else {
    for (i=0; i<dim; i++) {
      Fl[i] = 0.;
      Fr[i] = 1.;
      Fdelta[i] = 1.;
    }
  }

  /* allocate working space memory */
  X   = _unur_xmalloc( dim * sizeof(double));
  U   = _unur_xmalloc( dim * sizeof(double));
  bm  = _unur_xmalloc( dim * n_intervals * sizeof(int)); /* bins for marginal tests */
  memset(bm , 0, dim * n_intervals  * sizeof(int));

  /* calculate inverse Cholesky factor */
  if (L != NULL) {
    Linv  = _unur_xmalloc( dim * dim * sizeof(double));
    if (_unur_matrix_invert_matrix (dim, L, Linv, &Linv_det) != UNUR_SUCCESS) {
      _unur_error(test_name,UNUR_ERR_DISTR_DATA,"cannot compute inverse of Cholesky factor");
      pval_min = -2.; goto free_memory;
    }
  }

  /* now run generator */
  for( i=0; i<samplesize; i++ ) {

    /* get random vector X */
    _unur_sample_vec(gen, X);

    /* standardize vector: Z = L^{-1} (X - mean) */
    /* and transform to uniform U                */
    if (mean) { for (j=0; j<dim; j++)  X[j] -= mean[j]; }
    for (j=0; j<dim; j++) {
      double Z=0;
      if (Linv)
	for (k=0; k<=j; k++)  Z += Linv[idx(j,k)] * X[k]; 
      else
	Z = X[j];
      U[j] = (marginal_cdf[j](Z,marginals[j]) - Fl[j]) / Fdelta[j];
    }

    /* increase bins for tests for marginal distributions */
    for (j=0; j<dim; j++) {
      int iv;
      iv = (int)( n_intervals * U[j] );
      if (iv>=n_intervals) iv = n_intervals-1;
      if (iv < 0) iv = 0;
      bm[j*n_intervals + iv] += 1;
    }
  }

  /* ----------------------------------------------------------------------------*/

  /* make chi^2 test (marginal) */
  pval_min = 1.;
  for (j=0; j<dim; j++) {
    if (verbose >= 1) {
      fprintf(out,"\nChi^2-Test for marginal distribution [%d]\n",j);
    }

    pval = _unur_test_chi2test(NULL, bm+(j*n_intervals), n_intervals, classmin, verbose, out );
    pval_min = _unur_min(pval_min,pval);
  }
  
  if (verbose >= 1) {
    fprintf(out,"\nSummary:\n");
    fprintf(out,"  Minimal p-value * number_of_tests = %g:\n\n",pval_min * dim);
  }

  /* ----------------------------------------------------------------------------*/

free_memory:
  /* free memory */
  if (X)    free(X);
  if (U)    free(U);
  if (bm)   free(bm);
  if (Linv) free(Linv);
  if (marginals)  free (marginals);
  if (marginal_cdf)  free (marginal_cdf);
  if (Fl)   free(Fl);
  if (Fr)   free(Fr);
  if (Fdelta)   free(Fdelta);

  /* return result of test */
  return pval_min * dim;

#undef idx
#undef DISTR
} /* end of _unur_test_chi2_vec() */

/*---------------------------------------------------------------------------*/

double
_unur_test_chi2_cvemp ( struct unur_gen *gen,
          int n_intervals,
          int samplesize,
          int classmin,
          int verbose,
          FILE *out )
     /*----------------------------------------------------------------------*/
     /* Chi^2 test for continuous multivariate empirical distributions       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   n_intervals... number of intervals in which (0,1) is partitioned   */
     /*                  (each dimension is partitioned in n_intervals)      */
     /*   samplesize ... samplesize for test                                 */
     /*                  (if <= 0, CHI2_SAMPLEFAC * intervals^dim is used)   */
     /*   classmin   ... minimum number of expected occurrences for each class */
     /*                  (if <= 0, a default value is used.)                 */
     /*   verbose    ... verbosity level                                     */
     /*                  0 = no output                                       */
     /*                  1 = print summary                                   */
     /*                  2 = print classes and summary                       */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   minimal p-value * number_of_tests of all test statistics under H_0 */
     /*                                                                      */
     /* error:                                                               */
     /*   -2. ... missing data                                               */
     /*   -1. ... other errors                                               */
     /*----------------------------------------------------------------------*/
{
#define DISTR     gen->distr->data.cvec
#define idx(i,j)  ((i)*dim+(j))

  int dim;                       /* dimension of multivariate distribution */
  UNUR_DISTR *distr_normal;      /* object for standard normal distribution */
  UNUR_FUNCT_CONT *cdf;          /* pointer to CDF */
  double *X = NULL;              /* sampling vector */
  double *U = NULL;              /* X transformed to uniform */

  int *bm = NULL;                /* array for counting bins for marginals */
  double pval, pval_min;         /* p-value */

  int i, j;                      /* auxiliary variables */

  /* check arguments */
  CHECK_NULL(gen,-1.);
  /* we do not check magic cookies here */

  /* check given number of intervals */
  if (n_intervals <= 2)
    n_intervals = CHI2_INTERVALS_DEFAULT;

  /* samplesize */
  if( samplesize <= 0 ) samplesize = CHI2_DEFAULT_SAMPLESIZE;
  samplesize = _unur_min( samplesize, CHI2_MAX_SAMPLESIZE );

  /* dimension of distribution */
  dim = gen->distr->dim;  
  if (dim < 1) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"distribution dimension < 1 ?");
    return -1.;
  }

  /* we need the CDF for the standard normal distribution */
  distr_normal = unur_distr_normal( NULL, 0 );
  cdf = distr_normal->data.cont.cdf;

  /* allocate working space memory */
  X   = _unur_xmalloc( dim * sizeof(double));
  U   = _unur_xmalloc( dim * sizeof(double));
  bm  = _unur_xmalloc( dim * n_intervals * sizeof(int)); /* bins for marginal tests */
  memset(bm , 0, dim * n_intervals  * sizeof(int));

  /* now run generator */
  for( i=0; i<samplesize; i++ ) {

    /* get random vector X */
    _unur_sample_vec(gen, X);

    /* transform to uniform U */
    for (j=0; j<dim; j++)
      U[j] = cdf(X[j], distr_normal);

    /* increase bins for tests for marginal distributions */
    for (j=0; j<dim; j++) {
      int iv;
      iv = (int)( n_intervals * U[j] );
      if (iv>=n_intervals) iv = n_intervals-1;
      if (iv < 0) iv = 0;
      bm[j*n_intervals + iv] += 1;
    }
  }

  /* make chi^2 test (marginal) */
  pval_min = 1.;
  for (j=0; j<dim; j++) {
    if (verbose >= 1) {
      fprintf(out,"\nChi^2-Test for marginal distribution [%d]\n",j);
    }

    pval = _unur_test_chi2test(NULL, bm+(j*n_intervals), n_intervals, classmin, verbose, out );
    pval_min = _unur_min(pval_min,pval);
  }
  
  if (verbose >= 1) {
    fprintf(out,"\nSummary:\n");
    fprintf(out,"  Minimal p-value * number_of_tests = %g:\n\n",pval_min * dim);
  }

  /* free memory */
  _unur_distr_free(distr_normal);
  free(X);
  free(U);
  free(bm);

  /* return result of test */
  return pval_min * dim;

#undef idx
#undef DISTR
} /* end of _unur_test_chi2_cvemp() */

/*---------------------------------------------------------------------------*/

double
_unur_test_chi2test( double *prob,
         int *observed,
         int len,
         int classmin,
         int verbose,
         FILE *out )
     /*----------------------------------------------------------------------*/
     /* Chi^2 test for discrete distributions.                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   prob     ... probability vector (need not sum to 1)                */
     /*                  NULL indicates a uniform distribution on (0,len-1)  */
     /*   observed ... vector containing the observed number of occurrences  */
     /*   len      ... length of (both) vectors                              */
     /*   classmin ... minimum number of expected occurrences for each class */
     /*                  (if <= 0, a default value is used.)                 */
     /*   verbose  ... verbosity level                                       */
     /*                  0 = no output                                       */
     /*                  1 = print summary                                   */
     /*                  2 = print classes and summary                       */
     /*   out      ... output stream                                         */
     /*                                                                      */
     /* return:                                                              */
     /*   p-value of test statistics under H_0                               */
     /*                                                                      */
     /* error:                                                               */
     /*   -2. ... missing data                                               */
     /*   -1. ... other errors                                               */
     /*                                                                      */
     /* comment:                                                             */
     /*   joines classes if expected number of occurrences in some classes   */
     /*   is too small, i.e., less than "classmin".                          */
     /*----------------------------------------------------------------------*/
{
  double chi2 = 0.;     /* chi2-value */
  double df;            /* degrees of freedom for chi^2 distribution */
  double pval;          /* p-value */
  double clexpd = 0.;   /* expected occurrences in class */
  int clobsd = 0;       /* observed occurrences in class */
  int classes = 0;      /* (total) number of classes     */
  double probsum = 0.;  /* sum of all "probabilities" (need not be 1, for convenience) */
  int samplesize = 0;
  double factor;        /* factor for calculating expected number of occurrences */
  int i;

  UNUR_DISTR *distr_chisquare = NULL; /* distribution object to chi^2 distribution */

  /* check arguments */
  CHECK_NULL(observed,-1.);

  /* minimum number of occurrences in a class */
  classmin = (classmin > 0) ? classmin : CHI2_CLASSMIN_DEFAULT;

  /* compute sample size */
  for( samplesize=0, i=0; i<len; i++) 
    samplesize += observed[i];

  /* sum of probabilities (if not uniform distribution) */
  if (prob != NULL) {
    for( probsum=0., i=0; i<len; i++ )
      probsum += prob[i];
    factor = samplesize/probsum;
  }
  else   /* uniform distribution on (0,len-1) --> prob[i] = 1/len */
    factor = ((double)samplesize)/len;

  /* compute chi^2 value */
  clexpd = 0.;
  clobsd = 0;
  classes = 0;
  for( chi2=0., i=0; i<len; i++ ) {
    clexpd += (prob) ? prob[i]*factor : factor;  /* expected occurrences in this class */
    clobsd += observed[i];                       /* observed occurrences in this class */

    if (clexpd >= classmin || i == len-1) {
      /* number of expected occurrences is large enough or end of array */
      if (clobsd <= 0 && clexpd <= 0.) break;
      chi2 += (clobsd-clexpd)*(clobsd-clexpd)/clexpd;
      if (verbose >= 2)
	fprintf(out,"Class #%d:\tobserved %d\texpected %.2f\n",classes,clobsd,clexpd);
      clexpd = 0.;
      clobsd = 0;
      classes++;
    }
  }

  /* there must be at least two classes */
  if (classes < 2) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"too few classes!");
    if (verbose >= 1)
      fprintf(out,"\nCannot run chi^2-Test: too few classes\n");
    return -1.;
  }

  /* evaluate test statistics */
  /* chisquare distribution with df degrees of freedom */
  df = (double)(classes-1);
  distr_chisquare = unur_distr_chisquare( &df, 1 );
  /* p-value */
  if (distr_chisquare->data.cont.cdf) {
    pval = 1. - _unur_cont_CDF( chi2, distr_chisquare );
  }
  else {
    _unur_error(test_name,UNUR_ERR_GENERIC,"CDF for CHI^2 distribution required");
    pval = -2.;
  }
  _unur_distr_free(distr_chisquare);

  /* print result (if requested) */
  if (verbose >= 1 && pval >= 0.) {
    fprintf(out,"\nResult of chi^2-Test:\n  samplesize = %d\n",samplesize);
    fprintf(out,"  classes    = %d\t (minimum per class = %d)\n", classes, classmin);
    fprintf(out,"  chi2-value = %g\n  p-value    = %g\n\n", chi2, pval);
  }

  /* return result of test */
  return pval;

} /* end of _unur_test_chi2test() */

/*---------------------------------------------------------------------------*/
