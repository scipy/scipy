/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_exponential.c                                              *
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
 *  distr: Exponential distribution [2; ch.19, p.494]                        *
 *                                                                           *
 *  pdf:       f(x) = exp( - (x-theta)/sigma )                               *
 *  domain:    x >= theta                                                    *
 *  constant:  1 / sigma                                                     *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  sigma > 0  (1)  ... scale                                         *
 *     1:  theta      (0)  ... location                                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:     f(x) = exp(-x)                                                  *
 *  domain:  x >= 0                                                          *
 *                                                                           *
 *  parameters:                                                              *
 *     none                                                                  *
 *                                                                           *
 *     0:  sigma = 1                                                         *
 *     1:  theta = 0                                                         *
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
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"

/*---------------------------------------------------------------------------*/
static const char distr_name[] = "exponential";

/* parameters */
#define sigma  params[0]
#define theta  params[1]

#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)

/* function prototypes                                                       */
static double _unur_pdf_exponential( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_exponential( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_exponential( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_exponential( double x, const UNUR_DISTR *distr );
static double _unur_cdf_exponential( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_exponential( double u, const UNUR_DISTR *distr );

static int _unur_upd_mode_exponential( UNUR_DISTR *distr );
static int _unur_upd_area_exponential( UNUR_DISTR *distr );
static int _unur_set_params_exponential( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_exponential( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;

  if (DISTR.n_params > 0)
    /* standardize */
    x = (x - theta) / sigma;

  /* standard form */
  return ( (x<0.) ? 0. : exp(-x - LOGNORMCONSTANT) );

} /* end of _unur_pdf_exponential() */

/*---------------------------------------------------------------------------*/
  
double
_unur_logpdf_exponential( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;

  if (DISTR.n_params > 0)
    /* standardize */
    x = (x - theta) / sigma;

  /* standard form */
  return ( (x<0.) ? -INFINITY : (-x - LOGNORMCONSTANT) );

} /* end of _unur_logpdf_exponential() */

/*---------------------------------------------------------------------------*/
  
double
_unur_dpdf_exponential( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;

  if (DISTR.n_params > 0)
    /* standardize */
    x = (x - theta) / sigma;

  /* standard form */
  return ( (x<0.) ? 0. : -exp(-x - 2.*LOGNORMCONSTANT) );

} /* end of _unur_dpdf_exponential() */

/*---------------------------------------------------------------------------*/

double
_unur_dlogpdf_exponential( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;

  if (DISTR.n_params > 0)
    /* standardize */
    x = (x - theta) / sigma;

  /* standard form */
  return ( (x<0.) ? 0. : -1./sigma );

} /* end of _unur_dlogpdf_exponential() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_exponential( double x, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;

  if (DISTR.n_params > 0)
    /* standardize */
    x = (x - theta) / sigma;

  /* standard form */
  
  return ( (x<0.) ? 0. : 1.-exp(-x) );

} /* end of _unur_cdf_exponential() */

/*---------------------------------------------------------------------------*/

double
_unur_invcdf_exponential( double U, const UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;
  double X;

  X = - log( 1. - U );
  return ((DISTR.n_params==0) ? X : theta + sigma * X);
} /* end of _unur_invcdf_exponential() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_exponential( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.theta;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_exponential() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_exponential( UNUR_DISTR *distr )
{
  /* log of normalization constant */
  LOGNORMCONSTANT = log(DISTR.sigma);


  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }

  /* else */
  DISTR.area = ( _unur_cdf_exponential( DISTR.domain[1],distr) 
		 - _unur_cdf_exponential( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
  
} /* end of _unur_upd_area_exponential() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_exponential( UNUR_DISTR *distr, const double *params, int n_params )
{

  /* check number of parameters for distribution */
  if (n_params < 0) n_params = 0;
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  if (n_params > 0)
    CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter sigma */
  if (n_params > 0 && sigma <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"sigma <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form: none */

  /* default parameters */
  DISTR.sigma = 1.;
  DISTR.theta = 0.;

  /* copy optional parameters */
  switch (n_params) {
  case 2:
    DISTR.theta = theta;
  case 1:
    DISTR.sigma = sigma;
    n_params = 2;           /* number of parameters for non-standard form */
  default:
    break;
  }

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
  /* domain */
    DISTR.domain[0] = DISTR.theta;     /* left boundary  */
    DISTR.domain[1] = INFINITY;        /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_exponential() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_exponential( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_EXPONENTIAL;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  DISTR.init = _unur_stdgen_exponential_init;

  /* functions */
  DISTR.pdf     = _unur_pdf_exponential;     /* pointer to PDF                  */
  DISTR.logpdf  = _unur_logpdf_exponential;  /* pointer to logPDF               */
  DISTR.dpdf    = _unur_dpdf_exponential;    /* pointer to derivative of PDF    */
  DISTR.dlogpdf = _unur_dlogpdf_exponential; /* pointer to derivative of logPDF */
  DISTR.cdf     = _unur_cdf_exponential;     /* pointer to CDF                  */
  DISTR.invcdf  = _unur_invcdf_exponential;  /* pointer to inverse CDF          */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );
                
  /* set parameters for distribution */
  if (_unur_set_params_exponential(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant */
  LOGNORMCONSTANT = log(DISTR.sigma);

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.theta;   /* theta */
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_exponential;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_exponential; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_exponential; /* funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_exponential() */

/*---------------------------------------------------------------------------*/
#undef sigma 
#undef theta 
#undef DISTR
/*---------------------------------------------------------------------------*/
