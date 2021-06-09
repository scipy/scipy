/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_weibull.c                                                  *
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
 *  distr: Weibull distribution [2; ch.21, p.628]                            *
 *                                                                           *
 *  pdf:       f(x) = ((x-zeta)/alpha)^(c-1) * exp(-((x-zeta)/alpha)^c)      *
 *  domain:    zeta < x < infinity                                           *
 *  constant:  c / alpha                                                     *
 *                                                                           *
 *  parameters: 3                                                            *
 *     0:  c     > 0         ... shape                                       *
 *     1:  alpha > 0   (1)   ... scale                                       *
 *     2:  zeta        (0)   ... location                                    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:       f(x) = x^(c-1) * exp(-x^c)                                    *
 *  domain:    0 < x < infinity                                              *
 *  constant:  c                                                             *
 *                                                                           *
 *  parameters: 1                                                            *
 *     0:  c     >0    ... shape                                             *
 *                                                                           *
 *     1:  alpha = 1                                                         *
 *     2:  zeta  = 0                                                         *
 *                                                                           *
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
static const char distr_name[] = "weibull";

/* parameters */
#define c      params[0]
#define alpha  params[1]
#define zeta   params[2]

#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)

/* function prototypes                                                       */
static double _unur_pdf_weibull( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_weibull( double x, const UNUR_DISTR *distr );
static double _unur_cdf_weibull( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_weibull( double u, const UNUR_DISTR *distr );

static int _unur_upd_mode_weibull( UNUR_DISTR *distr );
static int _unur_upd_area_weibull( UNUR_DISTR *distr );
static int _unur_set_params_weibull( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_weibull( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;

  if (DISTR.n_params > 1)
    /* standardize */
    x = (x - zeta) / alpha;

  /* standard form */

  if (x < 0.)
    return 0.;

  if (_unur_iszero(x) && _unur_isone(c))
    return NORMCONSTANT;

  if (_unur_iszero(x) && !_unur_isone(c))
    return 0.;

  /* else */
  return (exp (-pow (x, c) + (c-1.) * log (x)) * NORMCONSTANT);

} /* end of _unur_pdf_weibull() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_weibull( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  register double factor = 1.;
  register double xc;

  if (DISTR.n_params > 1) {
    /* standardize */
    factor = 1. / alpha;
    x = (x - zeta) / alpha;
  }

  /* standard form */

  if (x < 0.)
    return 0.;
  
  if (_unur_iszero(x) && _unur_isone(c))
    return 0.; 

  /* else */
  xc = -pow (x, c);
  return (exp (xc + (c-2.) * log (x)) * (-1. - c * (-xc-1.)) * NORMCONSTANT * factor);

} /* end of unur_dpdf_weibull() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_weibull( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;

  if (DISTR.n_params > 1)
    /* standardize */
    x = (x - zeta) / alpha;

  /* standard form */

  if (x <= 0.)
    return 0.;

  /* else */
  return (1. - exp(-pow (x, c)));

} /* end of _unur_cdf_weibull() */

/*---------------------------------------------------------------------------*/

double
_unur_invcdf_weibull( double U, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  double X;

  X = pow( -log(1.-U), 1./c );
  return ((DISTR.n_params==1) ? X : zeta + alpha * X );
} /* end of _unur_invcdf_weibull() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_weibull( UNUR_DISTR *distr )
{
  DISTR.mode = (DISTR.c<=1.) ? 0. : DISTR.alpha * pow((DISTR.c - 1.)/DISTR.c, 1./DISTR.c) + DISTR.zeta;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_weibull() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_weibull( UNUR_DISTR *distr )
{
  /* normalization constant */
  NORMCONSTANT = DISTR.c / DISTR.alpha;

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }

  /* else */
  DISTR.area = ( _unur_cdf_weibull( DISTR.domain[1],distr) 
		 - _unur_cdf_weibull( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
  
} /* end of _unur_upd_area_weibull() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_weibull( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter c */
  if (c <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"c <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* check parameter alpha */
  if (n_params > 1 && alpha <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"alpha <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form */
  DISTR.c = c;

  /* default parameters */
  DISTR.alpha = 1.;
  DISTR.zeta  = 0.;

  /* copy optional parameters */
  switch (n_params) {
  case 3:
    DISTR.zeta = zeta;
  case 2:
    DISTR.alpha = alpha;
    n_params = 3;           /* number of parameters for non-standard form */
  default:
    break;
  }

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = DISTR.zeta;      /* left boundary  */
    DISTR.domain[1] = INFINITY;        /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_weibull() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_weibull( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_WEIBULL;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  /* DISTR.init = _unur_stdgen_weibull_init; */

  /* functions */
  DISTR.pdf    = _unur_pdf_weibull;    /* pointer to PDF               */
  DISTR.dpdf   = _unur_dpdf_weibull;   /* pointer to derivative of PDF */
  DISTR.cdf    = _unur_cdf_weibull;    /* pointer to CDF               */
  DISTR.invcdf = _unur_invcdf_weibull; /* pointer to inverse CDF       */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
 		 UNUR_DISTR_SET_MODE   |
  		 UNUR_DISTR_SET_PDFAREA );

  /* set parameters for distribution */
  if (_unur_set_params_weibull(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* normalization constant */
  NORMCONSTANT = DISTR.c / DISTR.alpha;

  /* mode and area below p.d.f. */
  DISTR.mode = (DISTR.c<=1.) ? 0. : DISTR.alpha * pow((DISTR.c - 1.)/DISTR.c, 1./DISTR.c) + DISTR.zeta;
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_weibull;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_weibull; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_weibull; /* funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_weibull() */

/*---------------------------------------------------------------------------*/
#undef c    
#undef alpha
#undef zeta 
#undef DISTR
/*---------------------------------------------------------------------------*/
