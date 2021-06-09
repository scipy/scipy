/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_logistic.c                                                 *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *   [3] N.L. Johnson, S. Kotz and N. Balakrishnan                           *
 *       Continuous Univariate Distributions,                                *
 *       Volume 2, 2nd edition                                               *
 *       John Wiley & Sons, Inc., New York, 1995                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  distr: Logistic distribution [3; ch.23, p.115]                           *
 *                                                                           *
 *  pdf:       f(x) = exp(-(x-alpha)/beta) * (1 + exp(-(x-alpha)/beta))^(-2) *
 *  cdf:       F(x) = (1 + exp(-(x-alpha)/beta))^(-1)                        *
 *  domain:    infinity < x < infinity                                       *
 *  constant:  1 / beta                                                      *
 *                                                                           *
 *  parameters: 2                                                            *
 *     0:  alpha       (0)  ... location                                     *
 *     1:  beta  > 0   (1)  ... scale                                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:       f(x) = exp(-x) * (1 + exp(-x))^(-2)                           *
 *  cdf:       F(x) = (1 + exp(-x))^(-1)                                     *
 *  domain:    infinity < x < infinity                                       *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters: none                                                         *
 *                                                                           *
 *     0:  alpha = 0                                                         *
 *     1:  beta  = 1                                                         *
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
static const char distr_name[] = "logistic";

/* parameters */
#define alpha  params[0]
#define beta   params[1]

#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)

/* function prototypes                                                       */
static double _unur_pdf_logistic( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_logistic( double x, const UNUR_DISTR *distr );
static double _unur_cdf_logistic( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_logistic( double u, const UNUR_DISTR *distr );

static int _unur_upd_mode_logistic( UNUR_DISTR *distr );
static int _unur_upd_area_logistic( UNUR_DISTR *distr );
static int _unur_set_params_logistic( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_logistic( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  register double ex;

  if (DISTR.n_params > 0)
    /* standardize */
    x = (x - alpha) / beta;

  /* standard form */

  ex = exp( -fabs(x) );

  return (NORMCONSTANT * ex / ((1. + ex) * (1. + ex)));

} /* end of _unur_pdf_logistic() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_logistic( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  register double factor = 1.;
  register double ex;

  if (DISTR.n_params > 0) {
    /* standardize */
    factor = 1. / beta;
    x = (x - alpha) / beta;
  }

  /* standard form */

  ex = exp(-fabs(x));
  if (x<0)
    factor = -factor;

  return (factor * NORMCONSTANT * ex * (ex - 1.) / ((1.+ex)*(1.+ex)*(1.+ex)));

} /* end of unur_dpdf_logistic() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_logistic( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;

  if (DISTR.n_params > 0)
    /* standardize */
    x = (x - alpha) / beta;

  /* standard form */

  return ( 1. / (1. + exp(-x)) );

} /* end of _unur_cdf_logistic() */

/*---------------------------------------------------------------------------*/

double
_unur_invcdf_logistic( double U, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  double X;

  X = -log(1./U - 1.);
  return ((DISTR.n_params==0) ? X : alpha + beta * X );
} /* end of _unur_invcdf_logistic() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_logistic( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.alpha;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_logistic() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_logistic( UNUR_DISTR *distr )
{
  /* normalization constant */
  NORMCONSTANT = 1. / DISTR.beta;

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }

  /* else */
  DISTR.area = ( _unur_cdf_logistic( DISTR.domain[1],distr) 
		 - _unur_cdf_logistic( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
  
} /* end of _unur_upd_area_logistic() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_logistic( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 0) n_params = 0;
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  if (n_params > 0)
    CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter sigma */
  if (n_params > 1 && beta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"beta <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form: none */

  /* default parameters */
  DISTR.alpha = 0.;
  DISTR.beta  = 1.;

  /* copy optional parameters */
  /* copy parameters */
  switch (n_params) {
  case 2:
    DISTR.beta = beta;
  case 1:
    DISTR.alpha = alpha;
    n_params = 2;           /* number of parameters for non-standard form */
  default:
    break;
  }

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -INFINITY;       /* left boundary  */
    DISTR.domain[1] = INFINITY;        /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_logistic() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_logistic( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_LOGISTIC;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  /* DISTR.init = _unur_stdgen_logistic_init; */

  /* functions */
  DISTR.pdf    = _unur_pdf_logistic;    /* pointer to PDF               */
  DISTR.dpdf   = _unur_dpdf_logistic;   /* pointer to derivative of PDF */
  DISTR.cdf    = _unur_cdf_logistic;    /* pointer to CDF               */
  DISTR.invcdf = _unur_invcdf_logistic; /* pointer to inverse CDF       */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
 		 UNUR_DISTR_SET_MODE   |
  		 UNUR_DISTR_SET_PDFAREA );

  /* set parameters for distribution */
  if (_unur_set_params_logistic(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_logistic;

  /* normalization constant */
  NORMCONSTANT = 1. / DISTR.beta;

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.alpha;
  DISTR.area = 1.;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_logistic; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_logistic; /* funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_logistic() */

/*---------------------------------------------------------------------------*/
#undef alpha
#undef beta 
#undef DISTR
/*---------------------------------------------------------------------------*/
