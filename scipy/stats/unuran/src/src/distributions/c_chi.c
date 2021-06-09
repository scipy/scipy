/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_chi.c                                                      *
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
 *  distr: Chi distribution [2; ch.18, p.417]                                *
 *  (Do not confuse with Chi^2 distribution!)                                *
 *                                                                           *
 *  pdf:       f(x) = x^(nu-1) * exp( -x^2/2 )                               *
 *  domain:    0 <= x < infinity                                             *
 *  constant:  1 / (2^((nu/2)-1) * Gamma(nu/2))                              *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  nu > 0   ... shape (degrees of freedom)                           *
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

static const char distr_name[] = "chi";

/*---------------------------------------------------------------------------*/
/* parameters */
#define nu  params[0]

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */
static double _unur_pdf_chi( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_chi( double x, const UNUR_DISTR *distr );
static double _unur_cdf_chi( double x, const UNUR_DISTR *distr );

static int _unur_upd_mode_chi( UNUR_DISTR *distr );
static int _unur_upd_area_chi( UNUR_DISTR *distr );
static int _unur_set_params_chi( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_chi(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (x <= 0.)
    /* out of support */
    return 0.;

  return (exp( log(x)*(nu - 1.) -x*x/2. - LOGNORMCONSTANT ));

} /* end of _unur_pdf_chi() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_chi(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (x <= 0.)
    /* out of support */
    return 0.;

  return ( exp(log(x)*(nu - 2.)-x*x/2. - LOGNORMCONSTANT) * (nu - 1. - x*x) );
} /* end of _unur_dpdf_chi() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_chi(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (x <= 0.)
    /* out of support of p.d.f. */
    return 0.;

  return _unur_SF_incomplete_gamma(x*x/2.,nu/2.);
} /* end of _unur_cdf_chi() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_chi( UNUR_DISTR *distr )
{
  DISTR.mode = (DISTR.nu >= 1.) ? sqrt(DISTR.nu - 1.) : 0.;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_chi() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_chi( UNUR_DISTR *distr )
{
  /* normalization constant */
  LOGNORMCONSTANT = _unur_SF_ln_gamma(DISTR.nu/2.) + M_LN2 * (DISTR.nu/2. - 1.);

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }

  /* else */
  DISTR.area = ( _unur_cdf_chi( DISTR.domain[1],distr) 
		 - _unur_cdf_chi( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
} /* end of _unur_upd_area_chi() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_chi( UNUR_DISTR *distr, const double *params, int n_params )
{

  /* check number of parameters for distribution */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 1) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 1; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter nu */
  if (nu <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"nu <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form */
  DISTR.nu = nu;

  /* copy optional parameters: none */

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0.;          /* left boundary  */
    DISTR.domain[1] = INFINITY;    /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_chi() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_chi( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_CHI;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = _unur_stdgen_chi_init;
   
  /* functions */
  DISTR.pdf  = _unur_pdf_chi;   /* pointer to PDF               */
  DISTR.dpdf = _unur_dpdf_chi;  /* pointer to derivative of PDF */
  DISTR.cdf  = _unur_cdf_chi;   /* pointer to CDF               */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA |
		 UNUR_DISTR_SET_MODE );
                
  /* set parameters for distribution */
  if (_unur_set_params_chi(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant */
  LOGNORMCONSTANT = _unur_SF_ln_gamma(DISTR.nu/2.) + M_LN2 * (DISTR.nu/2. - 1.);

  /* mode and area below p.d.f. */
  DISTR.mode = (DISTR.nu >= 1.) ? sqrt(DISTR.nu - 1.) : 0.;
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_chi;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_chi; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_chi; /* funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_chi() */

/*---------------------------------------------------------------------------*/
#undef nu
#undef DISTR
/*---------------------------------------------------------------------------*/
