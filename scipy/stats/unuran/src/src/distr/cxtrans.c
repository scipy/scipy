/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      cxtrans.c                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for manipulating distribution objects of      *
 *         id   CXTRANS  (continuous distribution of transformed RV)         *
 *         type CONT     (continuous univariate distribution)                *
 *                                                                           *
 *   return:                                                                 *
 *     UNUR_SUCCESS ... on success                                           *
 *     error code   ... on error                                             *
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
#include <distributions/unur_stddistr.h>
#include "distr.h"
#include "cxtrans.h"
#include "cont.h"
#include "distr_source.h"

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "transformed RV";

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cont    /* underlying (base) distribution          */
#define CXT   cxt->data.cont      /* distribution of transformed RV          */

#define BD_LEFT     domain[0]     /* left and ...                            */
#define BD_RIGHT    domain[1]     /* ... right boundary of domain            */

#define ALPHA       params[0]     /* parameter for transformation            */
#define MU          params[1]     /* relocation parameter for transformation */
#define SIGMA       params[2]     /* rescaling parameter for transformation  */
#define logPDFPOLE  params[3]     /* logPDF at pole                          */
#define dlogPDFPOLE params[4]     /* derivative of logPDF at pole            */

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */

static int _unur_distr_cxtrans_compute_domain( struct unur_distr *cxt );
/*---------------------------------------------------------------------------*/
/* compute domain for transformed RV                                         */
/*---------------------------------------------------------------------------*/

/* prototypes for CDF, PDF and its derviative for transformed RV             */
static double _unur_cdf_cxtrans( double x, const struct unur_distr *cxt );
static double _unur_pdf_cxtrans( double x, const struct unur_distr *cxt );
static double _unur_logpdf_cxtrans( double x, const struct unur_distr *cxt );
static double _unur_dpdf_cxtrans( double x, const struct unur_distr *cxt );
static double _unur_dlogpdf_cxtrans( double x, const struct unur_distr *cxt );

static double _unur_pdf_at_pole( const struct unur_distr *cxt );
static double _unur_dpdf_at_pole( const struct unur_distr *cxt );

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** transformed univariate continuous random RV                             **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_cxtrans_new( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* Create an object for the distribution of a transformed RV            */
     /* `distr' must be a pointer to a univariate continuous distribution.   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to  univariate continuous distribution.          */
     /*   alpha ... parameter of power transformation                        */
     /*             0.            -> logarithmic transformation              */
     /*             UNUR_INFINITY -> exponential transformation              */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to distribution object                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_distr *cxt;

  /* check arguments */
  _unur_check_NULL( distr_name,distr,NULL );

  /* check distribution */
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);

  /* get distribution object for generic continuous univariate distribution */
  cxt = unur_distr_cont_new();
  if (!cxt) return NULL;

  /* set id to distribution of order statistics */
  cxt->id = UNUR_DISTR_CXTRANS;

  /* name of distribution */
  cxt->name = distr_name;

  /* this is a derived distribution */
  /* clone base distribution ... */
  cxt->base = _unur_distr_cont_clone( distr );
  if (!cxt->base) { free(cxt); return NULL; }

  /* defaults: */
  CXT.n_params = 5;                 /* five parameters                       */
  CXT.ALPHA = 1.;                   /* parameter for transformation (default: identity) */
  CXT.MU = 0.;                      /* relocation parameter for transformation */
  CXT.SIGMA = 1.;                   /* rescaling parameter for transformation  */

  CXT.logPDFPOLE = -UNUR_INFINITY;  /* logPDF at pole                        */
  CXT.dlogPDFPOLE = UNUR_INFINITY;  /* derivative of logPDF at pole          */

  /* copy data */
  CXT.area = DISTR.area;            /* area below PDF (same as for distr)    */
  CXT.BD_LEFT = DISTR.BD_LEFT;      /* boundaries of ...                     */
  CXT.BD_RIGHT = DISTR.BD_RIGHT;    /* ... domain                            */
  CXT.mode = DISTR.mode;            /* mode of distribution                  */

  /* pointer to PDF, its derivative, and CDF */
  if (DISTR.cdf)     CXT.cdf = _unur_cdf_cxtrans;          /* pointer to CDF       */
  if (DISTR.pdf)     CXT.pdf = _unur_pdf_cxtrans;          /* pointer to PDF       */
  if (DISTR.logpdf)  CXT.logpdf = _unur_logpdf_cxtrans;    /* pointer to logPDF    */
  if (DISTR.dpdf)    CXT.dpdf = _unur_dpdf_cxtrans;        /* derivative of PDF    */
  if (DISTR.dlogpdf) CXT.dlogpdf = _unur_dlogpdf_cxtrans;  /* derivative of logPDF */

  /* parameters set */
  cxt->set = distr->set;

  /* return pointer to object */
  return cxt;

} /* end of unur_distr_cxtrans_new() */

/*---------------------------------------------------------------------------*/

const struct unur_distr *
unur_distr_cxtrans_get_distribution( const struct unur_distr *cxt )
     /*----------------------------------------------------------------------*/
     /* get pointer to distribution object for underlying distribution       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   cxt ... pointer to distribution of transformed RV                  */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to underlying distribution                                 */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( distr_name, cxt, NULL );
  _unur_check_distr_object( cxt, CONT, NULL );

  /* check distribution */
  if (cxt->id != UNUR_DISTR_CXTRANS) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,"");
    return NULL;
  }

  return cxt->base;
} /* end of unur_distr_cxtrans_get_distribution() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cxtrans_set_alpha( struct unur_distr *cxt, double alpha )
     /*----------------------------------------------------------------------*/
     /* change parameter of transformation                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   cxt   ... pointer to distribution of transformed RV                */
     /*   alpha ... parameter of power transformation                        */
     /*             0.            -> logarithmic transformation              */
     /*             UNUR_INFINITY -> exponential transformation              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double alpha_bak;

  /* check arguments */
  _unur_check_NULL( distr_name, cxt, UNUR_ERR_NULL );
  _unur_check_distr_object( cxt, CONT, UNUR_ERR_DISTR_INVALID );

  /* check distribution */
  if (cxt->id != UNUR_DISTR_CXTRANS) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return UNUR_ERR_DISTR_INVALID; }
  CHECK_NULL( cxt->base, UNUR_ERR_NULL );

  /* check parameter alpha */
  if (alpha < 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_SET,"alpha < 0");
    return UNUR_ERR_DISTR_SET;
  }

  if (_unur_iszero(alpha) && cxt->base->data.cont.BD_LEFT < 0. ) {
    /* logarithmic transformation */
    _unur_error(distr_name,UNUR_ERR_DISTR_SET,"invalid domain");
    return UNUR_ERR_DISTR_SET;
  }

  /* copy parameters */
  alpha_bak = CXT.ALPHA;  /* store old value of alpha */
  CXT.ALPHA = alpha;

  /* change domain of transformed RV */
  if (_unur_distr_cxtrans_compute_domain(cxt) != UNUR_SUCCESS) {
    /* invalid alpha */
    CXT.ALPHA = alpha_bak;
    return UNUR_ERR_DISTR_SET;
  }

  /* changelog */
  cxt->set &= ~UNUR_DISTR_SET_MODE; /* mode unknown */

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_cxtrans_set_alpha() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cxtrans_get_alpha( const struct unur_distr *cxt )
     /*----------------------------------------------------------------------*/
     /* get parameter of transformation                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   cxt   ... pointer to distribution of transformed RV                */
     /*                                                                      */
     /* return:                                                              */
     /*   parameter alpha ... on success                                     */
     /*   -UNUR_INFINITY  ... on error                                       */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( distr_name, cxt, -UNUR_INFINITY );
  _unur_check_distr_object( cxt, CONT, -UNUR_INFINITY );

  /* check distribution */
  if (cxt->id != UNUR_DISTR_CXTRANS) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return -UNUR_INFINITY; }

  /* o.k. */
  return CXT.ALPHA;
} /* end of unur_distr_cxtrans_get_alpha() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cxtrans_set_rescale( struct unur_distr *cxt, double mu, double sigma )
     /*----------------------------------------------------------------------*/
     /* change relocation and rescaling parameter                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   cxt   ... pointer to distribution of transformed RV                */
     /*   mu    ... relocation parameter                                     */
     /*   sigma ... rescaling parameter                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{

  double mu_bak;
  double sigma_bak;

  /* check arguments */
  _unur_check_NULL( distr_name, cxt, UNUR_ERR_NULL );
  _unur_check_distr_object( cxt, CONT, UNUR_ERR_DISTR_INVALID );

  /* check distribution */
  if (cxt->id != UNUR_DISTR_CXTRANS) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return UNUR_ERR_DISTR_INVALID; }
  CHECK_NULL( cxt->base, UNUR_ERR_NULL );

  /* check parameter sigma */
  if (sigma <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_SET,"sigma <= 0");
    return UNUR_ERR_DISTR_SET;
  }

  /* copy parameters */
  mu_bak = CXT.MU;
  CXT.MU = mu;
  sigma_bak = CXT.SIGMA;
  CXT.SIGMA = sigma;

  /* change domain of transformed RV */
  if (_unur_distr_cxtrans_compute_domain(cxt) != UNUR_SUCCESS) {
    /* invalid alpha */
    CXT.MU = mu_bak;
    CXT.SIGMA = sigma_bak;
    return UNUR_ERR_DISTR_SET;
  }

  /* changelog */
  cxt->set &= ~UNUR_DISTR_SET_MODE; /* mode unknown */

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_cxtrans_set_rescale() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cxtrans_get_mu( const struct unur_distr *cxt )
     /*----------------------------------------------------------------------*/
     /* get relocation parameter                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   cxt   ... pointer to distribution of transformed RV                */
     /*                                                                      */
     /* return:                                                              */
     /*   parameter mu    ... on success                                     */
     /*   -UNUR_INFINITY  ... on error                                       */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( distr_name, cxt, -UNUR_INFINITY );
  _unur_check_distr_object( cxt, CONT, -UNUR_INFINITY );

  /* check distribution */
  if (cxt->id != UNUR_DISTR_CXTRANS) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return -UNUR_INFINITY; }

  /* o.k. */
  return CXT.MU;
} /* end of unur_distr_cxtrans_get_mu() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cxtrans_get_sigma( const struct unur_distr *cxt )
     /*----------------------------------------------------------------------*/
     /* get rescaling parameter                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   cxt   ... pointer to distribution of transformed RV                */
     /*                                                                      */
     /* return:                                                              */
     /*   parameter sigma ... on success                                     */
     /*   -UNUR_INFINITY  ... on error                                       */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( distr_name, cxt, -UNUR_INFINITY );
  _unur_check_distr_object( cxt, CONT, -UNUR_INFINITY );

  /* check distribution */
  if (cxt->id != UNUR_DISTR_CXTRANS) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return -UNUR_INFINITY; }

  /* o.k. */
  return CXT.SIGMA;
} /* end of unur_distr_cxtrans_get_sigma() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cxtrans_set_logpdfpole( struct unur_distr *cxt, double logpdfpole, double dlogpdfpole )
     /*----------------------------------------------------------------------*/
     /* set value for logarithm of PDF and its derivative at poles of the    */
     /* underlying distribution                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   cxt         ... pointer to distribution of transformed RV          */
     /*   logpdfpole  ... value of logPDF at pole                            */
     /*   dlogpdfpole ... value of derviate of logPDF at pole                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( distr_name, cxt, UNUR_ERR_NULL );
  _unur_check_distr_object( cxt, CONT, UNUR_ERR_DISTR_INVALID );

  /* check distribution */
  if (cxt->id != UNUR_DISTR_CXTRANS) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return UNUR_ERR_DISTR_INVALID; }

  /* changelog */
  cxt->set |= UNUR_DISTR_SET_GENERIC; 

  /* copy parameters */
  CXT.logPDFPOLE = logpdfpole;
  CXT.dlogPDFPOLE = dlogpdfpole;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_distr_cxtrans_set_logpdfpole() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cxtrans_set_domain( struct unur_distr *cxt, double left, double right )
     /*----------------------------------------------------------------------*/
     /* set the left and right borders of the domain of the distribution     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   cxt      ... pointer to distribution of transformed RV             */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   the new boundary points may be +/- INFINITY                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, cxt, UNUR_ERR_NULL );
  _unur_check_distr_object( cxt, CONT, UNUR_ERR_DISTR_INVALID );

  /* check distribution */
  if (cxt->id != UNUR_DISTR_CXTRANS) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return UNUR_ERR_DISTR_INVALID; }

  if (_unur_isinf(CXT.ALPHA)==1) {
    /* exponential transformation */
    if (left < CXT.MU) {
      _unur_error(NULL,UNUR_ERR_DISTR_SET,"domain, left < 0");
      return UNUR_ERR_DISTR_SET;
    }
  }

  return unur_distr_cont_set_domain( cxt, left, right );

} /* end of unur_distr_cxtrans_set_domain() */

/*---------------------------------------------------------------------------*/

int _unur_distr_cxtrans_compute_domain( struct unur_distr *cxt )
     /*----------------------------------------------------------------------*/
     /* compute domain for transformed RV                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   cxt   ... pointer to distribution of transformed RV                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double left, right;
  double left_new, right_new;
  double alpha;

  /* check arguments */
  CHECK_NULL( cxt, UNUR_ERR_NULL );
  CHECK_NULL( cxt->base, UNUR_ERR_NULL );

  /* check distribution */
  if (cxt->id != UNUR_DISTR_CXTRANS) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return UNUR_ERR_DISTR_INVALID; }

  /* domain of underlying distribution */
  left  = cxt->base->data.cont.BD_LEFT;
  right = cxt->base->data.cont.BD_RIGHT;

  alpha = CXT.ALPHA;

  /* exponential transformation */
  if (_unur_isinf(alpha)==1) {
    left_new  = (_unur_isfinite(left)) ? exp(left) : 0.;
    right_new = exp(right);
  }

  /* logarithmic transformation */
  else if (_unur_iszero(alpha)) {
    if (left < 0. ) {
      _unur_error(distr_name,UNUR_ERR_DISTR_SET,"invalid domain");
      return UNUR_ERR_DISTR_SET;
    }
    left_new  = (left<=0.) ? -UNUR_INFINITY : log(left);
    right_new = log(right);
  }

  /* power transformation */
  else if (alpha > 0.) {
    left_new  = (left>=0.)  ? pow(left,alpha)  : -pow(-left,alpha);
    right_new = (right>=0.) ? pow(right,alpha) : -pow(-right,alpha);
  }

  /* else: error */
  else {
    _unur_error(distr_name,UNUR_ERR_SHOULD_NOT_HAPPEN,""); 
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }

  /* check for NaNs */
  if (_unur_isnan(left_new) || _unur_isnan(right_new)) {
      _unur_error(distr_name,UNUR_ERR_DISTR_SET,"NaN in now domain boundaries");
      return UNUR_ERR_DISTR_SET;
  }

  /* store new domain */
  /* also set boundary for truncated distributions */

  CXT.trunc[0] = CXT.domain[0] = left_new;     /* left boundary of domain  */
  CXT.trunc[1] = CXT.domain[1] = right_new;    /* right boundary of domain */

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_distr_cxtrans_compute_domain() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** CDF, PDF, and its derivative of transformed RV                          **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

#define CDF(x)      ((*(cxt->base->data.cont.cdf))    ((x), cxt->base))
#define PDF(x)      ((*(cxt->base->data.cont.pdf))    ((x), cxt->base))
#define logPDF(x)   ((*(cxt->base->data.cont.logpdf)) ((x), cxt->base))
#define dPDF(x)     ((*(cxt->base->data.cont.dpdf))   ((x), cxt->base))
#define dlogPDF(x)  ((*(cxt->base->data.cont.dlogpdf))((x), cxt->base))

/* power transformation */
#define POW(x)      ((x>=0.) ? pow(x,1./alpha) : -pow(-x,1./alpha))
#define dPOW(x)     ( pow(fabs(x), 1./alpha-1.) / alpha )          /* alpha != 1. */
#define dlogPOW(x)  ( (1./alpha-1.)*log(fabs(x)) - log(alpha) )
#define ddPOW(x)    ( ((x>=0.)?(1.-alpha):(alpha-1.)) * (_unur_isfsame(alpha,0.5)?1.0:pow(fabs(x),1./alpha-2.)) / (alpha*alpha) )

/* rescale */
#define rescale(x)  (CXT.SIGMA * (x) + CXT.MU)
/* #define drescale(x) (CXT.SIGMA * (x)) */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_cxtrans( double x, const struct unur_distr *cxt )
     /*
	CDF(x) = F(phi^{-1}(x))

	F(.)   ... CDF of underlying distribution
	phi(.) ... transformation
     */
{
  double alpha, s, mu;

  /* check arguments */
  CHECK_NULL( cxt, INFINITY );
  CHECK_NULL( cxt->base, INFINITY );
  CHECK_NULL( cxt->base->data.cont.cdf, INFINITY );

  alpha = CXT.ALPHA;
  s = CXT.SIGMA;
  mu = CXT.MU;

  /* exponential transformation */
  if (_unur_isinf(alpha)==1) {
    return ((x<=0.) ? 0. : CDF(s*log(x)+mu));
  }

  /* logarithmic transformation */
  if (_unur_iszero(alpha)) {
    return CDF(s*exp(x)+mu);
  }

  /* power transformation */
  if (alpha > 0.) {
    return CDF(s*POW(x)+mu);
  }

  /* else: error */
  _unur_error(distr_name,UNUR_ERR_SHOULD_NOT_HAPPEN,""); 
  return INFINITY;

} /* end of _unur_cdf_cxtrans() */

/*---------------------------------------------------------------------------*/

double
_unur_pdf_cxtrans( double x, const struct unur_distr *cxt )
     /*
	PDF(x) = f(phi^{-1}(x))*(phi^{-1})'(x)

	f(.)   ... PDF of underlying distribution
	phi(.) ... transformation
     */
{
  double alpha, s, mu;

  /* check arguments */
  CHECK_NULL( cxt, INFINITY );
  CHECK_NULL( cxt->base, INFINITY );
  CHECK_NULL( cxt->base->data.cont.pdf, INFINITY );

  alpha = CXT.ALPHA;
  s = CXT.SIGMA;
  mu = CXT.MU;

  /* exponential transformation */
  if (_unur_isinf(alpha)==1) {
    if (x<=0.) 
      return -INFINITY;
    else {
      /* PDF(log(x))/x */
      double fx = PDF(s*log(x)+mu);
      return (_unur_isfinite(fx) ? fx * s/x : _unur_pdf_at_pole(cxt));
    }
  }

  /* logarithmic transformation */
  if (_unur_iszero(alpha)) {
    /* PDF(exp(x)) * exp(x) */
    double ex = s * exp(x) + mu;
    if (! _unur_isfinite(ex)) {
      /* the PDF must be zero at +-infinity */
      return 0.;
    }
    else {
      double fx = PDF(ex);
      /* if PDF(ex) is not finite, we assume that it is a pole */
      return (_unur_isfinite(fx) ? fx * s*ex :  _unur_pdf_at_pole(cxt));
    }
  }

  /* identical transformation */
  if (_unur_isone(alpha)) {
    double fx = PDF(s*x+mu);
    return (_unur_isfinite(fx) ? s*fx :  _unur_pdf_at_pole(cxt));
  }

  /* power transformation */
  if (alpha > 0.) {
    double phix = s * POW(x) + mu;
    if (! _unur_isfinite(phix)) {
      /* the PDF must be zero at +-infinity */
      return 0.;
    }
    else {
      double fx = PDF(phix);
      if (_unur_isfinite(fx) && (!_unur_iszero(x) || alpha < 1.)) {
	double fcx =  fx * s * dPOW(x);
	/* if f(phix) is finite but fx*dPOW(x) is not,                 */
	/* we assume that the PDF of the transformed variable is zero. */
	/* (This case is very unlikely to happen, but we should be     */
	/* prepared for round-off error of the FPA.)                   */
	return (_unur_isfinite(fcx) ? fcx : 0.);
      }
      else 
	/* if PDF(phix) is not finite, we assume that it is a pole */
	return  _unur_pdf_at_pole(cxt);
    }
  }

  /* else: error */
  _unur_error(distr_name,UNUR_ERR_SHOULD_NOT_HAPPEN,""); 
  return INFINITY;

} /* end of _unur_pdf_cxtrans() */

/*---------------------------------------------------------------------------*/

double
_unur_logpdf_cxtrans( double x, const struct unur_distr *cxt )
     /* 
	logPDF(x) = logf(phi^{-1}(x)) + log((phi^{-1})'(x))

	logf(.) ... logPDF of underlying distribution
	phi(.)  ... transformation
     */
{
  double alpha, s, logs, mu;

  /* check arguments */
  CHECK_NULL( cxt, INFINITY );
  CHECK_NULL( cxt->base, INFINITY );
  CHECK_NULL( cxt->base->data.cont.logpdf, INFINITY );

  alpha = CXT.ALPHA;
  s = CXT.SIGMA;
  mu = CXT.MU;
  logs = log(CXT.SIGMA);

  /* exponential transformation */
  if (_unur_isinf(alpha)==1) {
    if (x<=0.) 
      return -INFINITY;
    else {
      /* logPDF(log(x))-log(x) */
      double logx = log(x);
      double logfx = logPDF(s*logx+mu);
      return (_unur_isfinite(logfx) ? (logfx - logx + logs) : CXT.logPDFPOLE);
    }
  }

  /* logarithmic transformation */
  if (_unur_iszero(alpha)) {
    /* logPDF(exp(x)) + x */
    double ex = s * exp(x) + mu;
    if (! _unur_isfinite(ex)) {
      /* logPDF must be -infinity (the PDF must be zero) at +-infinity */
      return -INFINITY;
    }
    else {
      double logfx = logPDF(ex);
      /* if logPDF(logx) is not finite, we assume that it is a pole */
      return (_unur_isfinite(logfx) ? (logfx + x + logs) : CXT.logPDFPOLE);
    }
  }

  /* identical transformation */
  if (_unur_isone(alpha)) {
    double logfx = logPDF(s*x+mu);
    return (_unur_isfinite(logfx) ? (logfx + logs) : CXT.logPDFPOLE);
  }

  /* power transformation */
  if (alpha > 0.) {
    double phix = s * POW(x) + mu;
    if (! _unur_isfinite(phix)) {
      /* logPDF must be -infinity (the PDF must be zero) at +-infinity */
      return -INFINITY;
    }
    else {
      double logfx = logPDF(phix);
      if (_unur_isfinite(logfx) && (!_unur_iszero(x) || alpha < 1.)) {
	double logfcx =  logfx + logs + dlogPOW(x);
	/* if logf(phix) is finite but logfx+dlogPOW(x) is not,        */
	/* we assume that the PDF of the transformed variable is zero. */
	/* (This case is very unlikely to happen, but we should be     */
	/* prepared for round-off error of the FPA.)                   */
	return (_unur_isfinite(logfcx) ? logfcx : -INFINITY);
      }
      else 
	/* if PDF(phix) is not finite, we assume that it is a pole */
	return CXT.logPDFPOLE;
    }
  }

  /* else: error */
  _unur_error(distr_name,UNUR_ERR_SHOULD_NOT_HAPPEN,""); 
  return INFINITY;

} /* end of _unur_logpdf_cxtrans() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_cxtrans( double x, const struct unur_distr *cxt )
     /*
	PDF(x) = f'(phi^{-1}(x))*((phi^{-1})'(x))^2
                 + f(phi^{-1}(x))*(phi^{-1})''(x)

	f(.)   ... PDF of underlying distribution
	phi(.) ... transformation
     */
{
  double alpha, s, mu;

  /* check arguments */
  CHECK_NULL( cxt, INFINITY );
  CHECK_NULL( cxt->base, INFINITY );
  CHECK_NULL( cxt->base->data.cont.pdf, INFINITY );
  CHECK_NULL( cxt->base->data.cont.dpdf, INFINITY );

  alpha = CXT.ALPHA;
  s = CXT.SIGMA;
  mu = CXT.MU;

  /* exponential transformation */
  if (_unur_isinf(alpha)==1) {
    if (x<=0.) 
      return 0.;
    else {
      /* dPDF(log(x))/x^2 - PDF(log(x))/x^2 */
      double logx = s*log(x)+mu;
      double fx = PDF(logx);
      double dfx = dPDF(logx);
      return (_unur_isfinite(fx) ? s*(s*dfx - fx)/(x*x) : _unur_dpdf_at_pole(cxt));
    }
  }

  /* logarithmic transformation */
  if (_unur_iszero(alpha)) {
    /* dPDF(exp(x)) * exp(2*x) + PDF(exp(x)) * exp(x) */
    double ex = s*exp(x)+mu;
    if (! _unur_isfinite(ex)) {
      /* dPDF must be zero at +-infinity */
      return 0.;
    }
    else {
      double fx = PDF(ex);
      double dfx = dPDF(ex);
      double dfcx = s * (dfx * s*ex*ex + fx * ex);
      /* if PDF(ex) is not finite, we assume that it is a pole */
      if (! _unur_isfinite(fx) ) return _unur_dpdf_at_pole(cxt);
      /* if derivate of PDF(ex) is not finite (or NaN) we return +/- INFINITY */
      if (! _unur_isfinite(dfcx) ) return (dfx>0 ? INFINITY : -INFINITY);
      /* otherwise return computed value */
      return dfcx;
    }
  }

  /* identical transformation */
  if (_unur_isone(alpha)) {
    double fx = PDF(s*x+mu);
    double dfx = dPDF(s*x+mu);
    return (_unur_isfinite(fx) ? s*dfx : _unur_dpdf_at_pole(cxt));
  }

  if (alpha > 0.) {
    /* power transformation */
    double phix = s*POW(x)+mu;
    if (! _unur_isfinite(phix)) {
      /* dPDF must be zero at +-infinity */
      return 0.;
    }
    else {
      double fx = PDF(phix);
      double dfx = dPDF(phix);
      double dphix = dPOW(x);
      double ddphix = ddPOW(x);
      if (_unur_isfinite(fx) && (!_unur_iszero(x) || alpha <= 0.5)) {
	double dfcx = s*(dfx * s*dphix*dphix + fx * s*ddphix);
	/* if f(phix) is finite but dfcx is not, we assume that dPDF */
	/* of the transformed variable is zero.                      */
	/* (This case is very unlikely to happen, but we should be   */
	/* prepared for round-off error of the FPA.)                 */
	return (_unur_isfinite(dfcx) ? dfcx : 0.);
      }
      else
	/* if PDF(phix) is not finite, we assume that it is a pole */
	return _unur_dpdf_at_pole(cxt);
    }
  }

  /* else: error */
  _unur_error(distr_name,UNUR_ERR_SHOULD_NOT_HAPPEN,""); 
  return INFINITY;

} /* end of _unur_dpdf_cxtrans() */

/*---------------------------------------------------------------------------*/

double
_unur_dlogpdf_cxtrans( double x, const struct unur_distr *cxt )
     /*
	PDF(x) = logf'(phi^{-1}(x)) * (phi^{-1})'(x)
	         + (phi^{-1})''(x) / (phi^{-1})'(x)
	       = logf'(phi^{-1}(x)) * (phi^{-1})'(x)
	         + (log(phi^{-1})'(x))'

	logf(.) ... logPDF of underlying distribution
	phi(.)  ... transformation
     */
{
  double alpha, s, mu;

  /* check arguments */
  CHECK_NULL( cxt, INFINITY );
  CHECK_NULL( cxt->base, INFINITY );
  CHECK_NULL( cxt->base->data.cont.logpdf, INFINITY );
  CHECK_NULL( cxt->base->data.cont.dlogpdf, INFINITY );

  alpha = CXT.ALPHA;
  s = CXT.SIGMA;
  mu = CXT.MU;

  /* exponential transformation */
  if (_unur_isinf(alpha)==1) {
    if (x<=0.) 
      return -INFINITY;
    else {
      /* (dlogPDF(log(x)) - 1) / x */
      double logx = s*log(x)+mu;
      double logfx = logPDF(logx);
      double dlogfx = dlogPDF(logx);
      /* if logPDF(logx) is not finite, we assume that it is a pole */
      return (_unur_isfinite(logfx) ? ((s*dlogfx-1)/x) : CXT.dlogPDFPOLE);
    }
  }

  /* logarithmic transformation */
  if (_unur_iszero(alpha)) {
    /* dlogPDF(exp(x))*exp(x) + 1 */
    double ex = s*exp(x)+mu;
    if (! _unur_isfinite(ex)) {
      /* dlogPDF must be -/+ infinity at +/-infinity */
      return (x>1. ? -INFINITY : INFINITY);
    }
    else {
      double logfx = logPDF(ex);
      double dlogfx = dlogPDF(ex);
      return (_unur_isfinite(logfx) ? s*dlogfx*ex + 1 : CXT.dlogPDFPOLE);
    }
  }

  /* identical transformation */
  if (_unur_isone(alpha)) {
    double logfx = logPDF(x);
    return (_unur_isfinite(logfx) ? s*dlogPDF(x) : CXT.dlogPDFPOLE);
  }

  if (alpha > 0.) {
    /* power transformation */
    double phix = s*POW(x)+mu;
    if (! _unur_isfinite(phix)) {
      /* dlogPDF must be -/+ infinity at +/-infinity */
      return ((x>1. || (x>-1. && x < 0.)) ? -INFINITY : INFINITY);
    }
    else {
      double logfx = logPDF(phix);
      if (_unur_isfinite(logfx) && (!_unur_iszero(x) || alpha <= 1.)) {
	double dlogfcx = ((x>=0.)?1.:-1.) * (dlogPDF(phix) * s*dPOW(x) + (1./alpha-1.)/x);
	if (! _unur_isfinite(dlogfcx)) {
	  /* if logf(phix) is finite but dlogfcx is not, we assume that dPDF */
	  /* of the transformed variable is zero.                            */
	  /* (This case is very unlikely to happen, but we should be         */
	  /* prepared for round-off error of the FPA.)                       */
	  return ((x>1. || (x>-1. && x < 0.)) ? -INFINITY : INFINITY);
	}
	return dlogfcx;
      }
      else
	/* if PDF(phix) is not finite, we assume that it is a pole */
	return CXT.dlogPDFPOLE;
    }
  }

  /* else: error */
  _unur_error(distr_name,UNUR_ERR_SHOULD_NOT_HAPPEN,""); 
  return INFINITY;

} /* end of _unur_dlogpdf_cxtrans() */

/*---------------------------------------------------------------------------*/

double 
_unur_pdf_at_pole( const struct unur_distr *cxt )
     /* return PDF at pole using logPDF at pole */
{
  return exp(CXT.logPDFPOLE);
} /* end of _unur_pdf_at_pole() */

/*---------------------------------------------------------------------------*/
double
_unur_dpdf_at_pole( const struct unur_distr *cxt )
     /* return derivative of PDF at pole using derivative of logPDF at pole */
{
  double fx = _unur_pdf_at_pole(cxt);

  if (! (_unur_isfinite(CXT.logPDFPOLE) && _unur_isfinite(fx)) )
    return UNUR_INFINITY;
  else 
    return (fx * CXT.dlogPDFPOLE);
} /* end of _unur_dpdf_at_pole() */

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
_unur_distr_cxtrans_debug( const struct unur_distr *cxt, const char *genid )
     /*----------------------------------------------------------------------*/
     /* write info about distribution into LOG file                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   cxt   ... pointer to distribution of transformed RV                */
     /*   genid ... pointer to generator id                                  */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(cxt,RETURN_VOID);
  COOKIE_CHECK(cxt,CK_DISTR_CONT,RETURN_VOID);
  CHECK_NULL(cxt->base,RETURN_VOID);

  LOG = unur_get_stream();

  /* print data about distribution */
  fprintf(LOG,"%s: distribution:\n",genid);
  fprintf(LOG,"%s:\ttype = continuous univariate distribution of transformed random variable\n",genid);
  fprintf(LOG,"%s:\tname = %s\n",genid,cxt->name);
  fprintf(LOG,"%s:\talpha = %g\t",genid,CXT.ALPHA);
  if (_unur_isinf(CXT.ALPHA)==1)
    fprintf(LOG,"[ exponential transformation: Y = exp(Z) ]\n"); 
  else if (_unur_iszero(CXT.ALPHA))
    fprintf(LOG,"[ logarithmic transformation: Y = log(Z) ]\n"); 
  else
    fprintf(LOG,"[ power transformation: Y = Z^alpha ]\n"); 
  fprintf(LOG,"%s:\tmu = %g, sigma = %g\t[Z = (X-%g)/%g]\n",genid, CXT.MU, CXT.SIGMA, CXT.MU, CXT.SIGMA);
  fprintf(LOG,"%s:\n",genid);

  fprintf(LOG,"%s:\tvalues used at pole of underlying distribution\n",genid);
  fprintf(LOG,"%s:\t\tlogPDF  = %g\t(PDF  = %g)",genid, CXT.logPDFPOLE, _unur_pdf_at_pole(cxt));
  _unur_print_if_default(cxt,UNUR_DISTR_SET_GENERIC);
  fprintf(LOG,"\n");
  fprintf(LOG,"%s:\t\tdlogPDF = %g\t(dPDF = %g)",genid, CXT.dlogPDFPOLE, _unur_dpdf_at_pole(cxt));
  _unur_print_if_default(cxt,UNUR_DISTR_SET_GENERIC);
  fprintf(LOG,"\n");

  if (cxt->set & UNUR_DISTR_SET_MODE)
    fprintf(LOG,"%s:\tmode = %g\n",genid,CXT.mode);
  else
    fprintf(LOG,"%s:\tmode unknown\n",genid);

  fprintf(LOG,"%s:\tdomain = (%g, %g)",genid,CXT.BD_LEFT,CXT.BD_RIGHT);
  _unur_print_if_default(cxt,UNUR_DISTR_SET_DOMAIN);

  fprintf(LOG,"\n%s:\tarea below PDF = %g",genid,CXT.area);
  _unur_print_if_default(cxt,UNUR_DISTR_SET_PDFAREA);
  fprintf(LOG,"\n%s:\n",genid);

  /* print data about underlying distribution */
  fprintf(LOG,"%s: Underlying distribution:\n",genid);
  _unur_distr_cont_debug(cxt->base, genid);

} /* end of _unur_distr_cxtrans_debug() */

/*---------------------------------------------------------------------------*/
#endif    /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
