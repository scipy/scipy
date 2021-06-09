/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      corder.c                                                     *
 *                                                                           *
 *   manipulate univariate continuous order statistics                       *
 *                                                                           *
 *   return:                                                                 *
 *     UNUR_SUCCESS ... on success                                           *
 *     error code   ... on error                                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2010  Wolfgang Hoermann and Josef Leydold            *
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
#include <specfunct/unur_specfunct_source.h>
#include <distributions/unur_stddistr.h>
#include "distr.h"
#include "cont.h"
#include "corder.h"
#include "distr_source.h"

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "order statistics";

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cont    /* underlying (base) distribution          */
#define OS    os->data.cont       /* order statistics                        */
#define LOGNORMCONSTANT (os->data.cont.norm_constant)

/*---------------------------------------------------------------------------*/

/* function prototypes                                                       */
static double _unur_pdf_corder( double x, const struct unur_distr *os );
static double _unur_dpdf_corder( double x, const struct unur_distr *os );

static double _unur_cdf_corder( double x, const struct unur_distr *os );
static int _unur_upd_area_corder( struct unur_distr *os );

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** univariate continuous order statistics                                  **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_corder_new( const struct unur_distr *distr, int n, int k )
     /*----------------------------------------------------------------------*/
     /* Create an object for order statistics of for a                       */
     /* sample size n and rank k.                                            */
     /* `distr' must be a pointer to a univariate continuous distribution.   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to  univariate continuous distribution.          */
     /*   n     ... sample size                                              */
     /*   k     ... rank                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to distribution object                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_distr *os;

  /* check arguments */
  _unur_check_NULL( distr_name,distr,NULL );

  /* check distribution */
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);

  if (distr->id == UNUR_DISTR_CORDER) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,"No order statistics of order statistics allowed");
    return NULL; 
  }

  /* check parameters n and k */
  if (n < 2 || k < 1 || k > n) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,"n < 2 or k < 1 or k > n");
    return NULL;
  }

  /* get distribution object for generic continuous univariate distribution */
  os = unur_distr_cont_new();
  if (!os) return NULL;

  /* set id to distribution of order statistics */
  os->id = UNUR_DISTR_CORDER;

  /* name of distribution */
  os->name = distr_name;

  /* this is a derived distribution */
  /* clone base distribution ... */
  os->base = _unur_distr_cont_clone( distr );
  if (!os->base) { free(os); return NULL; }

  /* set parameters for order statistics */
  OS.n_params = 2;                 /* two parameters: n and k                */
  OS.params[0] = (double) n;
  OS.params[1] = (double) k;

  /* copy data */
  OS.area = DISTR.area;               /* area below PDF (same as for distr)  */
  OS.trunc[0] = OS.domain[0] = DISTR.domain[0];  /* left boundary of domain  */
  OS.trunc[1] = OS.domain[1] = DISTR.domain[1];  /* right boundary of domain */
  
  /* pointer to PDF, its derivative, and CDF */
  if (DISTR.cdf) {
    OS.cdf = _unur_cdf_corder;        /* pointer to CDF    */
    if (DISTR.pdf) {
      OS.pdf = _unur_pdf_corder;      /* pointer to PDF    */
      if (DISTR.dpdf)
	OS.dpdf = _unur_dpdf_corder;  /* derivative of PDF */
    }
  }

  /* there is no necessity for a function that computes the area below PDF   */
  OS.upd_area  = _unur_upd_area_corder;

  /* parameters set */
  os->set = distr->set & ~UNUR_DISTR_SET_MODE; /* mode not derived from distr */

  /* log of normalization constant */
  if (_unur_upd_area_corder(os)==UNUR_SUCCESS)
    os->set |= UNUR_DISTR_SET_PDFAREA;

  /* return pointer to object */
  return os;

} /* end of unur_distr_corder_new() */

/*---------------------------------------------------------------------------*/

const struct unur_distr *
unur_distr_corder_get_distribution( const struct unur_distr *os )
     /*----------------------------------------------------------------------*/
     /* get pointer to distribution object for underlying distribution       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   os ... pointer to order statistics                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to underlying distribution                                 */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( distr_name, os, NULL );
  _unur_check_distr_object( os, CONT, NULL );

  /* check distribution */
  if (os->id != UNUR_DISTR_CORDER) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_INVALID,"");
    return NULL; 
  }

  return os->base;
} /* end of unur_distr_corder_get_distribution() */

/*---------------------------------------------------------------------------*/

int
unur_distr_corder_set_rank( struct unur_distr *os, int n, int k )
     /*----------------------------------------------------------------------*/
     /* change sample size n and rank k of order statistics.                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   os    ... pointer to order statistics                              */
     /*   n     ... sample size                                              */
     /*   k     ... rank                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( distr_name, os, UNUR_ERR_NULL );
  _unur_check_distr_object( os, CONT, UNUR_ERR_DISTR_INVALID );

  /* check distribution */
  if (os->id != UNUR_DISTR_CORDER) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return UNUR_ERR_DISTR_INVALID; }
  COOKIE_CHECK(os,CK_DISTR_CONT,UNUR_ERR_COOKIE);

  /* check parameters n and k */
  if (n < 2 || k < 1 || k > n) {
    _unur_error(distr_name,UNUR_ERR_DISTR_SET,"n < 2 or k < 1 or k > n");
    return UNUR_ERR_DISTR_SET;
  }

  /* changelog */
  os->set &= ~UNUR_DISTR_SET_MODE; /* mode unknown */

  /* copy parameters */
  OS.params[0] = (double) n;
  OS.params[1] = (double) k;

  /* log of normalization constant */
  _unur_upd_area_corder(os);

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_corder_set_rank() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_corder_get_rank( const struct unur_distr *os, int *n, int *k )
     /*----------------------------------------------------------------------*/
     /* get sample size n and rank k of order statistics.                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   os    ... pointer to order statistics                              */
     /*   n     ... sample size                                              */
     /*   k     ... rank                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( distr_name, os, UNUR_ERR_NULL );
  _unur_check_distr_object( os, CONT, UNUR_ERR_DISTR_INVALID );

  /* check distribution */
  if (os->id != UNUR_DISTR_CORDER) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return UNUR_ERR_DISTR_INVALID; }
  COOKIE_CHECK(os,CK_DISTR_CONT,UNUR_ERR_COOKIE);

  /* copy parameters */
  *n = (int)(OS.params[0] + 0.5);
  *k = (int)(OS.params[1] + 0.5);

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_corder_get_rank() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** PDF, its derivative and CDF of order statistics                         **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

double
_unur_pdf_corder( double x, const struct unur_distr *os )
     /* 
	PDF(x) = b(F(x)) * f(x)

	b(.) ... PDF of beta(k,n-k+1) distribution
	f(.) ... PDF of underlying distribution
	F(.) ... CDF of underlying distribution
     */
{ 
  double Fx;    /* CDF of underlying distribution at x */
  double fx;    /* PDF of underlying distribution at x */
  double p,q;   /* shape parameters for beta distribution */

  /* check arguments */
  _unur_check_NULL( NULL, os, INFINITY );
  _unur_check_distr_object( os, CONT, INFINITY );
  CHECK_NULL( os->base, INFINITY );
  _unur_check_distr_object( os->base, CONT, INFINITY );
  CHECK_NULL( os->base->data.cont.cdf, INFINITY );
  CHECK_NULL( os->base->data.cont.pdf, INFINITY );

  Fx = (*(os->base->data.cont.cdf)) (x, os->base);
  fx = (*(os->base->data.cont.pdf)) (x, os->base);

  p = OS.params[1];                       /* k     */
  q = OS.params[0] - OS.params[1] + 1.;   /* n-k+1 */

  /* PDF(x) = b(F(x)) * f(x) */
  if (fx <= 0. || Fx <= 0. || Fx >= 1.)
    return 0.;
  else
    return exp(log(fx) + (p-1.)*log(Fx) + (q-1.)*log(1.-Fx) - LOGNORMCONSTANT);

} /* end of _unur_pdf_corder() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_corder( double x, const struct unur_distr *os )
     /* 
	PDF'(x) = b'(F(x)) * f(x)^2 + b(F(x)) * f'(x)

	b(.) ... PDF of beta(k,n-k+1) distribution
	f(.) ... PDF of underlying distribution
	F(.) ... CDF of underlying distribution
     */
{
  double Fx;    /* CDF of underlying distribution at x */
  double fx;    /* PDF of underlying distribution at x */
  double dfx;   /* derivative of PDF of underlying distribution at x */
  double p,q;   /* shape parameters for beta distribution */
  double dpdf;  /* derivative of PDF of order statistics */
  double lFx, lFy;

  /* check arguments */
  _unur_check_NULL( NULL, os, INFINITY );
  _unur_check_distr_object( os, CONT, INFINITY );
  CHECK_NULL( os->base, INFINITY );
  _unur_check_distr_object( os->base, CONT, INFINITY );
  CHECK_NULL( os->base->data.cont.cdf, INFINITY );
  CHECK_NULL( os->base->data.cont.pdf, INFINITY );
  CHECK_NULL( os->base->data.cont.dpdf, INFINITY );

  Fx = (*(os->base->data.cont.cdf)) (x, os->base);
  fx = (*(os->base->data.cont.pdf)) (x, os->base);
  dfx = (*(os->base->data.cont.dpdf)) (x, os->base);

  p = OS.params[1];                       /* k     */
  q = OS.params[0] - OS.params[1] + 1.;   /* n-k+1 */

  if (fx <= 0. || Fx <= 0. || Fx >= 1.)
    return 0.;

  /* else */

  lFx = log(Fx);
  lFy = log(1.-Fx);

  /* PDF'(x) = b'(F(x)) * f(x)^2 + b(F(x)) * f'(x) */
  dpdf = ( exp(2.*log(fx) + (p-2.)*lFx + (q-2.)*lFy - LOGNORMCONSTANT)
	   * ( (p-1.)*(1.-Fx) - (q-1.)*Fx ));
  dpdf += exp((p-1.)*lFx + (q-1.)*lFy - LOGNORMCONSTANT) * dfx;

  return dpdf;
} /* end of _unur_dpdf_corder() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_corder( double x, const struct unur_distr *os ) 
     /* 
	CDF(x) = B(F(x))

	B(.) ... CDF of beta(k,n-k+1) distribution
	F(.) ... CDF of underlying distribution
     */
{
  double Fx;    /* CDF of underlying distribution at x */
  double p,q;   /* shape parameters for beta distribution */

  /* check arguments */
  _unur_check_NULL( NULL, os, INFINITY );
  _unur_check_distr_object( os, CONT, INFINITY );
  CHECK_NULL( os->base, INFINITY );
  _unur_check_distr_object( os->base, CONT, INFINITY );
  CHECK_NULL( os->base->data.cont.cdf, INFINITY );

  Fx = (*(os->base->data.cont.cdf)) (x, os->base);

  p = OS.params[1];                       /* k     */
  q = OS.params[0] - OS.params[1] + 1.;   /* n-k+1 */

  /* CDF(x) = B(F(x)) */
  return _unur_SF_incomplete_beta(Fx,p,q);

} /* end of _unur_cdf_corder() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_corder( UNUR_DISTR *os )
{

  /* log of normalization constant */
  /* LOGNORMCONSTANT = _unur_SF_ln_gamma(k) + _unur_SF_ln_gamma(n-k+1) - _unur_SF_ln_gamma(n+1); */
  LOGNORMCONSTANT = ( _unur_SF_ln_gamma(OS.params[1]) 
		      + _unur_SF_ln_gamma(OS.params[0] - OS.params[1] + 1.) 
		      - _unur_SF_ln_gamma(OS.params[0] + 1.) );

  /* we assume that the PDF is the derivative of the given CDF ! */

  /* standard distribution --> nothing to do */

  if (!(os->set & UNUR_DISTR_SET_STDDOMAIN)) {
    /* truncated distributions */
    if (OS.cdf == NULL)
      /* no CDF */
      return UNUR_ERR_DISTR_REQUIRED;

    OS.area  = (OS.domain[1] < INFINITY)  ? _unur_cdf_corder(OS.domain[1],os) : 1.;
    OS.area -= (OS.domain[0] > -INFINITY) ? _unur_cdf_corder(OS.domain[0],os) : 0.;
  }

  return (OS.area > 0.) ? UNUR_SUCCESS : UNUR_ERR_DISTR_DATA;

} /* end of _unur_upd_area_corder() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** debug                                                                   **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_distr_corder_debug( const struct unur_distr *os, const char *genid )
     /*----------------------------------------------------------------------*/
     /* write info about distribution into LOG file                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   os    ... pointer to order statistics                              */
     /*   genid ... pointer to generator id                                  */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(os,RETURN_VOID);
  COOKIE_CHECK(os,CK_DISTR_CONT,RETURN_VOID);
  CHECK_NULL(os->base,RETURN_VOID);

  LOG = unur_get_stream();

  /* print data about order statistics */
  fprintf(LOG,"%s: distribution:\n",genid);
  fprintf(LOG,"%s:\ttype = order statistics of continuous univariate distribution\n",genid);
  fprintf(LOG,"%s:\tname = %s\n",genid,os->name);
  fprintf(LOG,"%s:\tsample size\tn = %d\n",genid,(int)(OS.params[0]+0.5));
  fprintf(LOG,"%s:\trank\t\tk = %d\n",genid,(int)(OS.params[1]+0.5));
  fprintf(LOG,"%s:\n",genid);

  if (os->set & UNUR_DISTR_SET_MODE)
    fprintf(LOG,"%s:\tmode = %g\n",genid,OS.mode);
  else
    fprintf(LOG,"%s:\tmode unknown\n",genid);

  fprintf(LOG,"%s:\tdomain = (%g, %g)",genid,OS.domain[0],OS.domain[1]);
  _unur_print_if_default(os,UNUR_DISTR_SET_DOMAIN);

  fprintf(LOG,"\n%s:\tarea below PDF = %g",genid,OS.area);
  _unur_print_if_default(os,UNUR_DISTR_SET_PDFAREA);
  fprintf(LOG,"\n%s:\n",genid);

  /* print data about underlying distribution */
  fprintf(LOG,"%s: Underlying distribution:\n",genid);
  _unur_distr_cont_debug(os->base, genid);

} /* end of _unur_distr_corder_debug() */

/*---------------------------------------------------------------------------*/
#endif    /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
