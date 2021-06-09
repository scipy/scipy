/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_ig.c                                                       *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *   [2] N.L. Johnson, S. Kotz and N. Balakrishnan                           *
 *       Continuous Univariate Distributions,                                *
 *       Volume 1, 2nd edition                                               *
 *       John Wiley & Sons, Inc., New York, 1994                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  distr: Inverse Gaussian (Wald) distribution [2; ch.15, p.259]            *
 *                                                                           *
 *  pdf:       f(x) = sqrt( lambda/(2*pi*x^3) )                              *
 *                     * exp( -lambda*(x-mu)^2 / (2*mu^2*x) )                *
 *  domain:    0 < x < infinity                                              *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  mu     > 0    ... shape (mean)                                    *
 *     1:  lambda > 0    ... shape                                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2009-2010  Wolfgang Hoermann and Josef Leydold            *
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
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "ig";

/*---------------------------------------------------------------------------*/
/* parameters */
#define mu    params[0]
#define lambda params[1]

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */
static double _unur_pdf_ig( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_ig( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_ig( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_ig( double x, const UNUR_DISTR *distr );
static double _unur_cdf_ig( double x, const UNUR_DISTR *distr );

static int _unur_upd_mode_ig( UNUR_DISTR *distr );
static int _unur_upd_area_ig( UNUR_DISTR *distr );
static int _unur_set_params_ig( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_ig( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;

  if (x<=0.)
    return 0.;

  else
    return
      ( sqrt(lambda/(2*M_PI*x*x*x)) * exp( -lambda*(x-mu)*(x-mu) / (2*mu*mu*x) ) );
} /* end of _unur_pdf_ig() */

/*---------------------------------------------------------------------------*/

double
_unur_logpdf_ig( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;

  if (x<0.)
    return -INFINITY;

  else
    return ( 0.5* log ( lambda/(2*M_PI*x*x*x) )
	     -lambda*(x-mu)*(x-mu) / (2*mu*mu*x) );
} /* end of _unur_logpdf_ig() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_ig( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  double res;

  if (x<=0.)
    return 0.;

  /* else */
  res = -exp( -lambda*(x-mu)*(x-mu) / (2*mu*mu*x) );
  res *= sqrt(lambda / (x*x*x));
  res *= (3*mu*mu*x + lambda*(-mu*mu + x*x));
  res /= 2*mu*mu*sqrt(2*M_PI)*x*x;

  return res;

} /* end of _unur_dpdf_ig() */

/*---------------------------------------------------------------------------*/

double
_unur_dlogpdf_ig( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;

  if (x<=0.)
    return 0.;

  else 
    return 
      (0.5 * lambda*(-1./(mu*mu) + 1./(x*x)) - 3./x);

} /* end of _unur_dlogpdf_ig() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_ig( double x, const UNUR_DISTR *distr ) 
{
  register const double *params = DISTR.params;

#define Phi(x)   (_unur_SF_cdf_normal(x))

  if (x<=0.)
    return 0.;

  return 
    ( Phi(sqrt(lambda/x)*(x/mu-1.)) 
      + exp(2*lambda/mu) * Phi(-sqrt(lambda/x)*(x/mu+1.)) );

#undef Phi
} /* end of _unur_cdf_ig() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_ig( UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;

  DISTR.mode = 
    (-3.*mu*mu + mu*sqrt(4.*lambda*lambda + 9*mu*mu))/(2*lambda);

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_ig() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_ig( UNUR_DISTR *distr )
{
  /* log of normalization constant */
  LOGNORMCONSTANT = 0.;

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }

  /* else */
  DISTR.area = ( _unur_cdf_ig( DISTR.domain[1],distr) 
		 - _unur_cdf_ig( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} /* end of _unur_upd_area_ig() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_ig( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  /* check new parameter for generator */
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter a */
  if (mu <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"mu <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }
  if (lambda <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"lambda <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form */
  DISTR.mu     = mu;
  DISTR.lambda = lambda;

  /* default parameters: none */
  /* copy optional parameters: none */

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0.;             /* left boundary  */
    DISTR.domain[1] = INFINITY;        /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_ig() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Make distribution object                                               **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_ig( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_IG;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = NULL;         /*  _unur_stdgen_ig_init */

  /* functions */
  DISTR.pdf     = _unur_pdf_ig;     /* pointer to PDF                  */
  DISTR.logpdf  = _unur_logpdf_ig;  /* pointer to logPDF               */
  DISTR.dpdf    = _unur_dpdf_ig;    /* pointer to derivative of PDF    */
  DISTR.dlogpdf = _unur_dlogpdf_ig; /* pointer to derivative of logPDF */
  DISTR.cdf     = _unur_cdf_ig;     /* pointer to CDF                  */
  /* DISTR.invcdf  = _unur_invcdf_ig;  pointer to inverse CDF          */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );

  /* set parameters for distribution */
  if (_unur_set_params_ig(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant */
  LOGNORMCONSTANT = 0;

  /* mode and area below p.d.f. */
  _unur_upd_mode_ig(distr);
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_ig;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_ig; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_ig; /* funct for computing area */
                
  /* return pointer to object */
  return distr;

} /* end of unur_distr_ig() */

/*---------------------------------------------------------------------------*/
#undef mu
#undef lambda
#undef DISTR
/*---------------------------------------------------------------------------*/
