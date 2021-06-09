/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_cauchy.c                                                   *
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
 *  distr: Cauchy distribution [2; ch.16, p.299]                             *
 *                                                                           *
 *  pdf:       f(x) = 1./( 1 + ((x-theta)/lambda)^2 )                        *
 *  domain:    -infinity < x < infinity                                      *
 *  constant:  1 / (pi * lambda)                                             *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  theta       ... location                                          *
 *     1:  lambda > 0  ... scale                                             *
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

static const char distr_name[] = "cauchy";

/* parameters */
#define theta  params[0]
#define lambda params[1]

#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)

/* function prototypes                                                       */
static double _unur_pdf_cauchy( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_cauchy( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_cauchy( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_cauchy( double x, const UNUR_DISTR *distr );
static double _unur_cdf_cauchy( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_cauchy( double u, const UNUR_DISTR *distr );

static int _unur_upd_mode_cauchy( UNUR_DISTR *distr );
static int _unur_upd_area_cauchy( UNUR_DISTR *distr );
static int _unur_set_params_cauchy( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_cauchy(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (DISTR.n_params > 0)
    /* standardize */
    x = (x - theta) / lambda; 

  /* standard form */

  return (1./((1+x*x)*NORMCONSTANT));

} /* end of _unur_pdf_cauchy() */

/*---------------------------------------------------------------------------*/

double
_unur_logpdf_cauchy(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (DISTR.n_params > 0)
    /* standardize */
    x = (x - theta) / lambda; 

  /* standard form */

  return (-log1p(x*x)-log(NORMCONSTANT)); 

} /* end of _unur_logpdf_cauchy() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_cauchy(double x, const UNUR_DISTR *distr)
{
  register const double *params = DISTR.params;

  if (DISTR.n_params > 0)
    /* standardize */
    x = (x - theta) / lambda; 

  /* standard form */

  return ( -2.*x/(lambda*(1.+x*x)*(1.+x*x)*NORMCONSTANT) );

} /* end of _unur_dpdf_cauchy() */

/*---------------------------------------------------------------------------*/

double
_unur_dlogpdf_cauchy(double x, const UNUR_DISTR *distr)
{
  register const double *params = DISTR.params;

  if (DISTR.n_params > 0)
    /* standardize */
    x = (x - theta) / lambda; 

  /* standard form */

  return -2.*x/(lambda*(1.+x*x));

} /* end of _unur_dlogpdf_cauchy() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_cauchy(double x, const UNUR_DISTR *distr)
{
  register const double *params = DISTR.params;
  register double Fx;

  if (DISTR.n_params > 0)
    /* standardize */
    x = (x - theta) / lambda; 

  /* standard form */
  Fx = 0.5 + atan(x)/M_PI;

  /* correct round off error of atan() */
  if (Fx<0.)  Fx = 0.;
  if (Fx>1.)  Fx = 1.;

  return Fx;

} /* end of _unur_cdf_cauchy() */

/*---------------------------------------------------------------------------*/

double
_unur_invcdf_cauchy(double u, const UNUR_DISTR *distr)
{
  register const double *params = DISTR.params;
  double X;

  X = tan( M_PI * (u - 0.5) );
  return ((DISTR.n_params==0) ? X : theta + lambda * X );
} /* end of _unur_invcdf_cauchy() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_cauchy( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.theta; 

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_cauchy() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_cauchy( UNUR_DISTR *distr )
{
  /* normalization constant */
  NORMCONSTANT = M_PI * DISTR.lambda;

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }

  /* else */
  DISTR.area = ( _unur_cdf_cauchy( DISTR.domain[1],distr) 
		 - _unur_cdf_cauchy( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
  
} /* end of _unur_upd_area_cauchy() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_cauchy( UNUR_DISTR *distr, const double *params, int n_params )
{

  /* check number of parameters for distribution */
  if (n_params < 0) n_params = 0;
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  if (n_params > 0)
    CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter lambda */
  if (n_params == 2 && lambda <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"lambda <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form: none */

  /* default parameters */
  DISTR.theta  = 0.;
  DISTR.lambda = 1.;

  /* copy optional parameters */
  switch (n_params) {
  case 2:
    DISTR.lambda = lambda;
  case 1:
    DISTR.theta  = theta;
    n_params = 2;           /* number of parameters for non-standard form */
  default:
    break;
  }

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -INFINITY;   /* left boundary  */
    DISTR.domain[1] = INFINITY;    /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_cauchy() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_cauchy( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_CAUCHY;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  /* DISTR.init = _unur_stdgen_cauchy_init; */

  /* functions */
  DISTR.pdf     = _unur_pdf_cauchy;     /* pointer to PDF                  */
  DISTR.logpdf  = _unur_logpdf_cauchy;  /* pointer to logPDF               */
  DISTR.dpdf    = _unur_dpdf_cauchy;    /* pointer to derivative of PDF    */
  DISTR.dlogpdf = _unur_dlogpdf_cauchy; /* pointer to derivative of logPDF */
  DISTR.cdf     = _unur_cdf_cauchy;     /* pointer to CDF                  */
  DISTR.invcdf  = _unur_invcdf_cauchy;  /* pointer to inverse CDF          */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );
                
  /* set parameters for distribution */
  if (_unur_set_params_cauchy(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* normalization constant */
  NORMCONSTANT = M_PI * DISTR.lambda;

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.theta; 
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_cauchy;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_cauchy; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_cauchy; /* funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_cauchy() */

/*---------------------------------------------------------------------------*/
#undef theta 
#undef lambda
#undef DISTR
/*---------------------------------------------------------------------------*/
