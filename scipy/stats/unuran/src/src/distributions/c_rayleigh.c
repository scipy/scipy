/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_rayleigh.c                                                 *
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
 *  distr: Rayleigh distribution [2; ch.18, p.456]                           *
 *                                                                           *
 *  pdf:       f(x) = x * exp( -1/2 * (x/sigma)^2 )                          *
 *  domain:    0 <= x < infinity                                             *
 *  constant:  1 / sigma^2                                                   *
 *                                                                           *
 *  parameters: 1                                                            *
 *     0:  sigma > 0   ... scale                                             *
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
static const char distr_name[] =  "rayleigh";

/* parameters */
#define sigma  params[0]

#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)

/* function prototypes                                                       */
static double _unur_pdf_rayleigh( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_rayleigh( double x, const UNUR_DISTR *distr );
static double _unur_cdf_rayleigh( double x, const UNUR_DISTR *distr );

static int _unur_upd_mode_rayleigh( UNUR_DISTR *distr );
static int _unur_upd_area_rayleigh( UNUR_DISTR *distr );
static int _unur_set_params_rayleigh( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_rayleigh( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  
  if (x<=0.) 
    return 0.;

  /* else */
  return (x * exp(-x*x/(2.*sigma*sigma) - LOGNORMCONSTANT ) ); 

} /* end of _unur_pdf_rayleigh() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_rayleigh( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  register double z;

  z = x*x/(sigma*sigma);
  return ( (x<=0.) ? 0. : exp(-z/2 - LOGNORMCONSTANT) * (1-z) ); 
} /* end of _unur_dpdf_rayleigh() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_rayleigh( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  return ( (x<=0.) ? 0. : 1. - exp(-x*x/(2.*sigma*sigma)) );
} /* end of _unur_cdf_rayleigh() */

/*---------------------------------------------------------------------------*/


int
_unur_upd_mode_rayleigh( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.sigma;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_rayleigh() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_rayleigh( UNUR_DISTR *distr )
{
  /* log of normalization constant */
  LOGNORMCONSTANT =   2. * log(DISTR.sigma);

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }

  /* else */
  DISTR.area = ( _unur_cdf_rayleigh( DISTR.domain[1],distr) 
		 - _unur_cdf_rayleigh( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
  
} /* end of _unur_upd_area_rayleigh() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_rayleigh( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 1) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 1; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter sigma */
  if (sigma <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"sigma <= 0.");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form */
  DISTR.sigma = sigma;

  /* default parameters: none */
  /* copy optional parameters: none */

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0.;              /* left boundary  */
    DISTR.domain[1] = INFINITY;        /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_rayleigh() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_rayleigh( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_RAYLEIGH;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = NULL;    /* _unur_stdgen_rayleigh_init; */
                
  /* functions */
  DISTR.pdf  = _unur_pdf_rayleigh;  /* pointer to PDF               */
  DISTR.dpdf = _unur_dpdf_rayleigh; /* pointer to derivative of PDF */
  DISTR.cdf  = _unur_cdf_rayleigh;  /* pointer to CDF               */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );

  /* set parameters for distribution */
  if (_unur_set_params_rayleigh(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant */
  LOGNORMCONSTANT =   2. * log(DISTR.sigma);

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.sigma;
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_rayleigh;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_rayleigh; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_rayleigh; /* funct for computing area */

  /* return pointer to object */
  return distr; 

} /* end of unur_distr_rayleigh() */

/*---------------------------------------------------------------------------*/
#undef sigma
#undef DISTR
/*---------------------------------------------------------------------------*/
