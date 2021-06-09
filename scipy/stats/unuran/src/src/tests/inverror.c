/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file: inverror.c                                                        *
 *                                                                           *
 *   Estimate u and x-errors for inversion methods                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2008-2010 Wolfgang Hoermann and Josef Leydold             *
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
#include <methods/unur_methods_source.h>
#include <distr/distr_source.h>
#include <methods/x_gen.h>
#include <methods/x_gen_source.h>
#include <methods/cstd.h>
#include <methods/cstd_struct.h>
#include <methods/dgt.h>
#include <methods/dstd.h>
#include <methods/dstd_struct.h>
#include <methods/hinv.h>
#include <methods/mixt.h>
#include <methods/mixt_struct.h>
#include <methods/ninv.h>
#include <methods/pinv.h>
#include "unuran_tests.h"

/*---------------------------------------------------------------------------*/

static char test_name[] = "InvError";

static double uerror_cont ( const UNUR_GEN *gen, 
			    double *max_error, double *MAE, double threshold,
			    int samplesize, int randomized, int testtails, 
			    int verbosity, FILE *out );

static double uerror_discr( const UNUR_GEN *gen, 
			    double *max_error, double *MAE, double threshold,
			    int samplesize, int randomized, int testtails, 
			    int verbosity, FILE *out );

static double qrng (int i, int samplesize);

/*---------------------------------------------------------------------------*/

double
unur_test_u_error( const UNUR_GEN *gen, 
		   double *max_error, double *MAE, double threshold,
		   int samplesize, int randomized, int testtails, 
		   int verbosity, FILE *out )
     /*----------------------------------------------------------------------*/
     /* Estimate maximal u-error and mean absolute error (MAE) by means of   */
     /* (Quasi-) Monte-Carlo simulation.                                     */
     /* In addition a penalty value is computed and returned.                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   max_error  ... pointer for storing maximal u-error                 */
     /*   MEA        ... pointer for storing MA u-error                      */
     /*   threshold  ... maximum allowed error                               */
     /*   samplesize ... sample size for Monte Carlo simulation              */
     /*   randomized ... use pseudo-random (TRUE) or quasi-random (FALSE)    */
     /*   testtails  ... when TRUE then run a special test for tails         */
     /*   verbosity  ... verbosity level, 0 = no output, 1 = output          */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   penalty score                                                      */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL(test_name,gen,-1.);
  if (verbosity) { _unur_check_NULL(test_name,out,-1.); }
  if (samplesize < 1000) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"samplesize too small --> increased to 1000");
    samplesize = 1000;
  }

  /* type of generator object */
  if ( (gen->method == UNUR_METH_HINV) ||
       (gen->method == UNUR_METH_NINV) ||
       (gen->method == UNUR_METH_PINV) ||
       (gen->method == UNUR_METH_CSTD && ((struct unur_cstd_gen*)gen->datap)->is_inversion) ||
       (gen->method == UNUR_METH_MIXT && ((struct unur_mixt_gen*)gen->datap)->is_inversion)
       ) {
    /* continuous distribution */
    return uerror_cont(gen,max_error,MAE,threshold,samplesize,
		       randomized,testtails,verbosity,out);
  }
  if ( (gen->method == UNUR_METH_DGT) ||
       (gen->method == UNUR_METH_DSTD && ((struct unur_dstd_gen*)gen->datap)->is_inversion)
       ) {
    /* continuous distribution */
    return uerror_discr(gen,max_error,MAE,threshold,samplesize,
			randomized,testtails,verbosity,out);
  }

  _unur_error(test_name,UNUR_ERR_GENERIC,"inversion method required");
  return -1.;

} /* end of unur_test_estimate_u_error() */

/*---------------------------------------------------------------------------*/

double
uerror_cont( const UNUR_GEN *gen, 
	     double *max_error, double *MAE, double threshold,
	     int samplesize, int randomized, int testtails, 
	     int verbosity, FILE *out )
     /*----------------------------------------------------------------------*/
     /* Estimate u-error for continuous distributions.                       */
     /*                                                                      */
     /* return:                                                              */
     /*   penalty score                                                      */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
#define DISTR   gen->distr->data.cont
  
  double CDFmin, CDFmax;     /* minimum and maximum of CDF in given domain */

  double (*quantile)(const UNUR_GEN *, double);  /* pointer to quantile function */

  double U;                  /* uniform and non-uniform (Q)RN */
  double X;                  /* non-uniform RV */
  double cdfX;               /* CDF at X */

  double uerror, umax, usum; /* last, maximum, sum of U-error(s) */
  double penalty = 0;        /* score for error */

  int j;                     /* aux variable */

  /* get pointer to function that approximates quantiles */
  switch (gen->method) {
  case UNUR_METH_HINV:
    quantile = unur_hinv_eval_approxinvcdf;
    break;

  case UNUR_METH_NINV:
    quantile = unur_ninv_eval_approxinvcdf;
    break;

  case UNUR_METH_PINV:
    quantile = unur_pinv_eval_approxinvcdf;
    break;

  case UNUR_METH_CSTD:
    if (! (((struct unur_cstd_gen*)gen->datap)->is_inversion))
      return -1.;
    quantile = unur_cstd_eval_invcdf;
    break;

  case UNUR_METH_MIXT:
    if (! (((struct unur_mixt_gen*)gen->datap)->is_inversion))
      return -1.;
    quantile = unur_cstd_eval_invcdf;
    break;

  default:
    _unur_error(test_name,UNUR_ERR_GENERIC,"inversion method required");
    return -1.;
  }

  /* CDF required */
  if (DISTR.cdf == NULL) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"CDF required");
    return -2.;
  }

  /* range of CDF */
  CDFmin = (DISTR.trunc[0] > -INFINITY) ? _unur_cont_CDF((DISTR.trunc[0]),(gen->distr)) : 0.;
  CDFmax = (DISTR.trunc[1] < INFINITY)  ? _unur_cont_CDF((DISTR.trunc[1]),(gen->distr)) : 1.;

  /* initialize variables */
  umax = 0.;
  usum = 0.;

  /* get sample */
  for(j=0;j<samplesize;j++) {

    /* uniform random number */
    if (randomized)
      U = _unur_call_urng(gen->urng);
    else
      U = (testtails) ? qrng(j,samplesize) : (j+0.5) / ((double) samplesize);

    /* compute inverse CDF */
    X = quantile(gen,U);

    /* compute CDF at X */
    cdfX = _unur_cont_CDF(X,gen->distr);

    /* compute U-error:  U-error = | U - CDF(X) |
       However, we have first to rescale CDF(X):
          CDFrescaled = (CDF(X) - CDFmin) / (CDFmax - CDFmin)
       as we have to deal with truncated distribution.
       But then the computed U error might be too large; thus 
       we have to correct it by the factor (CDFmax - CDFmin)
       to unavoidably high U-errors.
     */
    uerror = fabs( U*(CDFmax - CDFmin) - (cdfX-CDFmin));

    /* update error estimates */
    usum += uerror;
    if (uerror > umax) {
      umax = uerror;
    }

    /* update penalty */
    if (_unur_FP_less(threshold,uerror)) {
      penalty += 1. + 10.*(uerror - threshold) / threshold;
      if (verbosity)
	fprintf(out,"\tmax u-error exceeded at %g: %g (>%g)\n",
		X,uerror,threshold);
    }
  }

  /* save data */
  *max_error = umax;
  *MAE = usum/samplesize;

  /* return penalty */
  return penalty/samplesize;

#undef DISTR
} /* end of uerror_cont() */

/*---------------------------------------------------------------------------*/

double
uerror_discr( const UNUR_GEN *gen, 
		double *max_error, double *MAE, double threshold,
		int samplesize, int randomized,
		int testtails ATTRIBUTE__UNUSED, 
		int verbosity, FILE *out )
     /*----------------------------------------------------------------------*/
     /* Estimate u-error for discrete distributions.                         */
     /*                                                                      */
     /* return:                                                              */
     /*   penalty score                                                      */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
#define DISTR   gen->distr->data.discr
  
  double CDFmin, CDFmax;     /* minimum and maximum of CDF in given domain */

  int (*quantile)(const UNUR_GEN *, double);  /* pointer to quantile function */

  double U;                  /* uniform and non-uniform (Q)RN */
  int K;                     /* non-uniform RV */

  double cdfK;               /* CDF at K */
  double uerror, umax, usum; /* uerror, maximum u-error and sum of all uerrors */
  double penalty = 0;        /* score for error */

  int j;                     /* aux variable */

  /* get pointer to function that approximates quantiles */
  switch (gen->method) {
  case UNUR_METH_DGT:
    quantile = unur_dgt_eval_invcdf;
    break;

  case UNUR_METH_DSTD:
    if (! (((struct unur_dstd_gen*)gen->datap)->is_inversion))
      return -1.;
    quantile = unur_dstd_eval_invcdf;
    break;

  default:
    _unur_error(test_name,UNUR_ERR_GENERIC,"inversion method required");
    return -1.;
  }

  /* CDF required */
  if (DISTR.cdf == NULL) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"CDF required");
    return -2.;
  }

  /* range of CDF */
  CDFmin = (DISTR.trunc[0] <= INT_MIN) ? 0. : _unur_discr_CDF((DISTR.trunc[0]),(gen->distr));
  CDFmax = _unur_discr_CDF((DISTR.trunc[1]),(gen->distr));

  /* initialize variables */
  umax = 0.;
  usum = 0.;

  /* get sample */
  for(j=0;j<samplesize;j++) {

    /* uniform random number */
    if (randomized)
      U = _unur_call_urng(gen->urng);
    else
      U = (j+0.5) / ((double) samplesize);

    /* compute inverse CDF */
    K = (int) quantile(gen,U);

    /* estimate uerror */
    uerror = 0.;
    cdfK = _unur_discr_CDF(K,gen->distr);
    if (cdfK < U) {
      uerror = U - cdfK;
    }
    else {
      cdfK = _unur_discr_CDF(K-1,gen->distr);
      uerror = cdfK - U;
      uerror = _unur_max(0.,uerror);
    }
    /* Remark: we could save a lot of time if we store CDF values in a table. */

    /* Remark: we do not consider truncated distributions here !! */

    /* update error estimates */
    usum += uerror;
    if (uerror > umax)
      umax = uerror;


    /* update error estimates */
    if (uerror > umax) {
      umax = uerror;
    }

    /* update penalty */
    if (_unur_FP_less(threshold,uerror)) {
      penalty += 1. + 10.*(uerror - threshold) / threshold;
      if (verbosity)
	fprintf(out,"\tmax u-error exceeded at U=%g: %g (>%g)\n",
		U,uerror,threshold);
    }
  }

  /* save data */
  *max_error = umax;
  *MAE = usum/samplesize;

  /* return penalty */
  return penalty/samplesize;

#undef DISTR
} /* end of uerror_discr() */

/*---------------------------------------------------------------------------*/

double 
qrng (int i, int samplesize)
     /*----------------------------------------------------------------------*/
     /* Quasi-random sequence for testing inversion error.                   */
     /* 90% of the sample size is used for equidistributed points in (0,1),  */
     /* 5% for left tail in (0,1e-5) and                                     */
     /* 5% for the right tail in (1-1e-5,1).                                 */ 
     /*                                                                      */
     /* parameters:                                                          */
     /*   i          ... index of pointe in sequence                         */
     /*   samplesize ... sample size for Monte Carlo simulation              */
     /*                                                                      */
     /* return:                                                              */
     /*   point in sequence                                                  */
     /*----------------------------------------------------------------------*/
{
  double U;
  int tail;

  /* theshold values for tail parts */
  tail = (int) (0.05 * samplesize);

  /* i must be in 0, ..., samplesize-1 */
  i = i % samplesize;

  if (i < tail) {
    /* lower tail */
    U = (i+0.5) / (1.e5 * tail); 
  }
  else if (i >= samplesize-tail) {
    /* upper tail */
    i -= samplesize-tail;
    U = 1. - (i+0.5) / (1.e5 * tail);
  }
  else {
    /* central part */
    i -= tail;
    U = (i+0.5) / (samplesize - 2.*tail);
  }

  /* return point */
  return U;

}  /* end of qrng() */

/*---------------------------------------------------------------------------*/
