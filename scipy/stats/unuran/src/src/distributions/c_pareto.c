/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_pareto.c                                                   *
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
 *  distr: Pareto distribution (of first kind) [2; ch.20, p.574]             *
 *                                                                           *
 *  pdf:       f(x) = x^(-(a+1))                                             *
 *  domain:    x >= k                                                        *
 *  constant:  a * k^a                                                       *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  k > 0   ... location, shape                                       *
 *     1:  a > 0   ... shape                                                 *
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
static const char distr_name[] = "pareto";

/* parameters */
#define k  params[0]
#define a  params[1]

#define DISTR distr->data.cont
/* #define NORMCONSTANT (distr->data.cont.norm_constant) */

/* function prototypes                                                       */
static double _unur_pdf_pareto( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_pareto( double x, const UNUR_DISTR *distr );
static double _unur_cdf_pareto( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_pareto( double u, const UNUR_DISTR *distr );

static int _unur_upd_mode_pareto( UNUR_DISTR *distr );
static int _unur_upd_area_pareto( UNUR_DISTR *distr );
static int _unur_set_params_pareto( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_pareto( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (x<k)
    return 0.;
  /* else */
  return ((a/k) / pow(x/k, a + 1.) );

/*    return ( (x<k) ? 0. : (a/k) / pow(x/k, a + 1.) ); */
} /* end of _unur_pdf_pareto() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_pareto( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  return ( (x<k) ? 0. : a * (-a-1.) / (k * k) * pow(x/k,-a-2.) );
} /* end of _unur_dpdf_pareto() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_pareto( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  return ( (x<k) ? 0. : (1. - pow(k/x,a)) );
} /* end of _unur_cdf_pareto() */

/*---------------------------------------------------------------------------*/

double
_unur_invcdf_pareto( double U, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  double X;

  X = pow(1-U, -1/a);
  X *= k;
  return X;
} /* end of _unur_invcdf_pareto() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_pareto( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.k;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_pareto() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_pareto( UNUR_DISTR *distr )
{
  /* normalization constant: none */

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }

  /* else */
  DISTR.area = ( _unur_cdf_pareto( DISTR.domain[1],distr) 
		 - _unur_cdf_pareto( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
  
} /* end of _unur_upd_area_pareto() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_pareto( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  /* check new parameter for generator */
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter k */
  if (k <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"k <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* check parameter a */
  if (a <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"a <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form */
  DISTR.k = k;
  DISTR.a = a;

  /* default parameters: none */
  /* copy optional parameters: none */

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = DISTR.k;         /* left boundary  */
    DISTR.domain[1] = INFINITY;        /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_pareto() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_pareto( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_PARETO;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  /* DISTR.init = _unur_stdgen_pareto_init; */

  /* functions */
  DISTR.pdf    = _unur_pdf_pareto;    /* pointer to PDF               */
  DISTR.dpdf   = _unur_dpdf_pareto;   /* pointer to derivative of PDF */
  DISTR.cdf    = _unur_cdf_pareto;    /* pointer to CDF               */
  DISTR.invcdf = _unur_invcdf_pareto; /* pointer to inverse CDF       */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
  		 UNUR_DISTR_SET_MODE   |
  		 UNUR_DISTR_SET_PDFAREA );

  /* set parameters for distribution */
  if (_unur_set_params_pareto(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* normalization constant: none */

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.k;
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_pareto;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_pareto; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_pareto; /* funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_pareto() */

/*---------------------------------------------------------------------------*/
#undef k
#undef a
#undef DISTR
/*---------------------------------------------------------------------------*/
