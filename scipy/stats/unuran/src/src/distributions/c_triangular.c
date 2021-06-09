/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_triangular.c                                               *
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
 *  distr: Triangular distribution [3; ch.26, p.297]                         *
 *                                                                           *
 *  pdf:    f(x) =  2*x / H           for 0 <= x <= H                        *
 *  pdf:    f(x) =  2*(1-x) / (1-H)   for H <= x <= 1                        *
 *  domain:    0 <= x <= 1                                                   *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters:  1                                                           *
 *     0:  H  (0 <= H <= 1)  ... shape                                       *
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
static const char distr_name[] = "triangular";

/* parameters */
#define H   params[0]    /* shape */

#define DISTR distr->data.cont
/* #define NORMCONSTANT (distr->data.cont.norm_constant) */

/* function prototypes                                                       */
static double _unur_pdf_triangular( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_triangular( double x, const UNUR_DISTR *distr );
static double _unur_cdf_triangular( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_triangular( double u, const UNUR_DISTR *distr );

static int _unur_upd_mode_triangular( UNUR_DISTR *distr );
static int _unur_upd_area_triangular( UNUR_DISTR *distr );
static int _unur_set_params_triangular( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_triangular( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;

  if (x <= 0.)
    return 0.;
  if (x <= H)
    return (2.*x/H);
  if (x < 1.)
    return (2.*(1.-x)/(1.-H));
  /* otherwise */
  return 0.;
} /* end of _unur_pdf_triangular() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_triangular( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;

  if (x < 0.)
    return 0.;
  if (x <= H && H > 0.)
    return (2./H);
  if (x <= 1.&& H < 1.)
    return (-2./(1.-H));
  /* otherwise */
  return 0.;
} /* end of unur_dpdf_triangular() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_triangular( double x, const UNUR_DISTR *distr )
{ 
  const double *params = DISTR.params;
  double Fx;

  if (x <= 0.)
    return 0.;
  if (x <= H)
    return (x*x/H);
  if (x < 1.) {
    if ((Fx = ((H + x * (x-2.))/(H-1.))) < 1.)
      return Fx;
  }
  /* otherwise */
  return 1.;
} /* end of _unur_cdf_triangular() */

/*---------------------------------------------------------------------------*/

double
_unur_invcdf_triangular( double U, const UNUR_DISTR *distr )
{ 
  const double *params = DISTR.params;
  double tmp,X;

  if (U<=H) {
    X = sqrt(H*U);
  }
  else {
    tmp = (1.-H)*(1.-U);
    X = (tmp>0.) ? (1.-sqrt(tmp)) : 1.;
  }
  return X;
} /* end of _unur_invcdf_triangular() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_triangular( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.H;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_triangular() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_triangular( UNUR_DISTR *distr )
{
  /* normalization constant: none */

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }

  /* else */
  DISTR.area = ( _unur_cdf_triangular( DISTR.domain[1],distr) 
		 - _unur_cdf_triangular( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
  
} /* end of _unur_upd_area_triangular() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_triangular( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 0) n_params = 0;
  if (n_params > 1) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 1; }
  if (n_params > 0)
    CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter H */
  if (n_params > 0 && (H < 0. || H > 1.)) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"H < 0 || H > 1");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form: none */

  /* default parameters */
  DISTR.H = 0.5;   /* default is symmetric triangular distribution */

  /* copy optional parameters */
  if (n_params == 1)
    DISTR.H = H;

  /* store total number of parameters */
  DISTR.n_params = 1;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0.;        /* left boundary  */
    DISTR.domain[1] = 1.;        /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_triangular() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_triangular( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_TRIANGULAR;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  /* DISTR.init = _unur_stdgen_triangular_init; */

  /* functions */
  DISTR.pdf    = _unur_pdf_triangular;    /* pointer to PDF               */
  DISTR.dpdf   = _unur_dpdf_triangular;   /* pointer to derivative of PDF */
  DISTR.cdf    = _unur_cdf_triangular;    /* pointer to CDF               */
  DISTR.invcdf = _unur_invcdf_triangular; /* pointer to inverse CDF       */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
 		 UNUR_DISTR_SET_MODE   |
  		 UNUR_DISTR_SET_PDFAREA );

  /* set parameters for distribution */
  if (_unur_set_params_triangular(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.H;
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_triangular;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_triangular; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_triangular; /* funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_triangular() */

/*---------------------------------------------------------------------------*/
#undef theta
#undef phi  
#undef DISTR
/*---------------------------------------------------------------------------*/
