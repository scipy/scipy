/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_gamma.c                                                    *
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
 *  distr: Gamma distribution [2; ch.17, p.337]                              *
 *                                                                           *
 *  pdf:       f(x) = ((x-gamma)/beta)^(alpha-1) * exp( -(x-gamma)/beta )    *
 *  domain:    x > gamma                                                     *
 *  constant:  1 / (beta * Gamma(alpha))                                     *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  alpha > 0       ... shape                                         *
 *     1:  beta > 0   (1)  ... scale                                         *
 *     2:  gamma      (0)  ... location                                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:       f(x) = x^(alpha-1) * exp(-x)                                  *
 *  domain:    x > 0                                                         *
 *  constant:  Gamma(alpha)                                                  *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  alpha > 0  ... shape                                              *
 *                                                                           *
 *     1:  beta  = 1                                                         *
 *     2:  gamma = 0                                                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold             *
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

static const char distr_name[] = "gamma";

/*---------------------------------------------------------------------------*/
/* parameters */
#define alpha  params[0]   /* shape */
#define beta   params[1]   /* scale */
#define gamma  params[2]   /* location */

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */
static double _unur_pdf_gamma( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_gamma( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_gamma( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_gamma( double x, const UNUR_DISTR *distr );
static double _unur_cdf_gamma( double x, const UNUR_DISTR *distr );
#ifdef _unur_SF_invcdf_gamma
static double _unur_invcdf_gamma( double x, const UNUR_DISTR *distr );
#endif

static int _unur_upd_mode_gamma( UNUR_DISTR *distr );
static int _unur_upd_area_gamma( UNUR_DISTR *distr );
static double _unur_lognormconstant_gamma(const double *params, int n_params);
static int _unur_set_params_gamma( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_gamma( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;

  if (DISTR.n_params > 1)
    /* standardize */
    x = (x-gamma) / beta;

  /* standard form */

  if (_unur_isone(alpha) && x >= 0.)
    return exp( -x - LOGNORMCONSTANT);

  if (x > 0.)
    /* pow(x,alpha-1.) * exp(-x) */
    return exp( (alpha-1.)*log(x) - x - LOGNORMCONSTANT);

  if (_unur_iszero(x))
    return (alpha>1. ? 0. : INFINITY);

  /* out of domain */
  return 0.;

} /* end of _unur_pdf_gamma() */

/*---------------------------------------------------------------------------*/

double
_unur_logpdf_gamma( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;

  if (DISTR.n_params > 1)
    /* standardize */
    x = (x-gamma) / beta;

  /* standard form */

  if (_unur_isone(alpha) && x >= 0.)
    return ( -x - LOGNORMCONSTANT);

  if (x > 0.)
    return ( (alpha-1.)*log(x) - x - LOGNORMCONSTANT);

  if (_unur_iszero(x))
    return (alpha>1. ? -INFINITY : INFINITY);

  /* out of domain */
  return -INFINITY;

} /* end of _unur_logpdf_gamma() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_gamma( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;

  if (DISTR.n_params > 1)
    /* standardize */
    x = (x-gamma) / beta;

  /* standard form */
  
  if (_unur_isone(alpha) && x>=0)
    return( -exp(-x - LOGNORMCONSTANT) / beta );
  
  if (x > 0.)
    return ( exp( log(x) * (alpha-2.) - x - LOGNORMCONSTANT) *  ((alpha-1.) -x) / beta ); 

  if (_unur_iszero(x) && alpha < 2.)
    return (alpha>1. ? INFINITY : -INFINITY);

  /* out of domain */
  return 0.;

} /* end of _unur_dpdf_gamma() */

/*---------------------------------------------------------------------------*/

double
_unur_dlogpdf_gamma( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;

  if (DISTR.n_params > 1)
    /* standardize */
    x = (x-gamma) / beta;

  /* standard form */

  if (_unur_isone(alpha) && x >= 0.)
    return -1./beta;

  if (x > 0.)
    return ((alpha-1.)/x - 1)/beta;

  if (_unur_iszero(x))
    return (alpha>1. ? INFINITY : -INFINITY);

  /* out of domain */
  return 0.;

} /* end of _unur_dlogpdf_gamma() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_gamma( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;

  if (DISTR.n_params > 1)
    /* standardize */
    x = (x-gamma) / beta;

  /* standard form */

  if (x <= 0.)
    return 0.;

  if (_unur_isinf(x)==1)
    return 1.;

  return _unur_SF_incomplete_gamma(x,alpha);

} /* end of _unur_cdf_gamma() */

/*---------------------------------------------------------------------------*/

#ifdef _unur_SF_invcdf_gamma

double
_unur_invcdf_gamma( double x, const UNUR_DISTR *distr )
{ 
  const double *params = DISTR.params;

  if (DISTR.n_params == 1)
    return _unur_SF_invcdf_gamma(x, alpha, 1.);
  else
    return (gamma + _unur_SF_invcdf_gamma(x, alpha, beta));

} /* end of _unur_invcdf_gamma() */

#endif

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_gamma( UNUR_DISTR *distr )
{
  register double *params = DISTR.params;

  DISTR.mode = (alpha >= 1.) ? (alpha - 1.) : 0.;

  if (DISTR.n_params > 1)
    /* de-standardize */
    DISTR.mode = DISTR.mode * beta + gamma;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  /* for shape parameters less than 1 we also should set */
  /* the center of the distribution to some value > 0.   */
  if (alpha < 1.) {
    double center = alpha * beta + gamma;
    center = _unur_max(center,DISTR.domain[0]);
    center = _unur_min(center,DISTR.domain[1]);
    unur_distr_cont_set_center(distr,center);
  }

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_gamma() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_gamma( UNUR_DISTR *distr )
{
  /* log of normalization constant */
  LOGNORMCONSTANT = _unur_lognormconstant_gamma(DISTR.params,DISTR.n_params);
  
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  
  /* else */
  DISTR.area = ( _unur_cdf_gamma( DISTR.domain[1],distr) 
		 - _unur_cdf_gamma( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;

} /* end of _unur_upd_area_gamma() */

/*---------------------------------------------------------------------------*/

double
_unur_lognormconstant_gamma( const double *params, int n_params )
{
  if (n_params > 1)
    /* non-standard form */
    return ( _unur_SF_ln_gamma(alpha) + log(beta) );

  else
    /* standard form */
    return (_unur_SF_ln_gamma(alpha));

} /* end of _unur_lognormconstant_gamma() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_gamma( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter alpha */
  if (alpha <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"alpha <= 0.");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* check parameter beta */
  if (n_params > 1 && beta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"beta <= 0.");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form */
  DISTR.alpha = alpha;

  /* default parameters */
  DISTR.beta  = 1.;
  DISTR.gamma = 0.;

  /* copy optional parameters */
  switch (n_params) {
  case 3:
    DISTR.gamma = gamma;
  case 2:
    DISTR.beta = beta;
    n_params = 3;           /* number of parameters for non-standard form */
  default:
    break;
  }

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = DISTR.gamma;  /* left boundary  */
    DISTR.domain[1] = INFINITY;     /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_gamma() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Make distribution object                                               **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_gamma( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_GAMMA;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = _unur_stdgen_gamma_init;

  /* functions */
  DISTR.pdf     = _unur_pdf_gamma;     /* pointer to PDF                  */
  DISTR.logpdf  = _unur_logpdf_gamma;  /* pointer to logPDF               */
  DISTR.dpdf    = _unur_dpdf_gamma;    /* pointer to derivative of PDF    */
  DISTR.dlogpdf = _unur_dlogpdf_gamma; /* pointer to derivative of logPDF */
  DISTR.cdf     = _unur_cdf_gamma;     /* pointer to CDF                  */
#ifdef _unur_SF_invcdf_gamma
  DISTR.invcdf  = _unur_invcdf_gamma;  /* pointer to inverse CDF          */
#endif

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA |
		 UNUR_DISTR_SET_MODE );
                
  /* set parameters for distribution */
  if (_unur_set_params_gamma(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant */
  LOGNORMCONSTANT = _unur_lognormconstant_gamma(DISTR.params,DISTR.n_params);

  /* mode and area below p.d.f. */
  _unur_upd_mode_gamma( distr );
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_gamma;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_gamma; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_gamma; /* funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_gamma() */

/*---------------------------------------------------------------------------*/
#undef alpha
#undef beta 
#undef gamma
#undef DISTR
/*---------------------------------------------------------------------------*/
