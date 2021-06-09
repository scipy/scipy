/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_lomax.c                                                    *
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
 * distr: Lomax distribution (Pareto distr. of second kind) [2; ch.20, p.575]*
 *                                                                           *
 *  pdf:       f(x) = (x+C)^(-(a+1))                                         *
 *  domain:    x >= 0                                                        *
 *  constant:  a * C^a                                                       *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  a > 0   ... shape                                                 *
 *     1:  C > 0   ... location                                              *
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
static const char distr_name[] = "lomax";

/* parameters */
#define a params[0]
#define C params[1]

#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)

/* function prototypes                                                       */
static double _unur_pdf_lomax( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_lomax( double x, const UNUR_DISTR *distr );
static double _unur_cdf_lomax( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_lomax( double u, const UNUR_DISTR *distr );

static int _unur_upd_mode_lomax( UNUR_DISTR *distr );
static int _unur_upd_area_lomax( UNUR_DISTR *distr );
static int _unur_set_params_lomax( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_lomax( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  if (x<0.) 
    return 0.;
  /* else */
  return ( pow(x+C,-(a+1.)) * NORMCONSTANT );

/*    return ( (x<0.) ? 0. : pow(x+C,-(a+1.)) * NORMCONSTANT ); */
} /* end of _unur_pdf_lomax() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_lomax( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;

  return ( (x<0.) ? 0. : -(a+1.) * pow(x+C,-(a+2.)) * NORMCONSTANT );
} /* end of _unur_dpdf_lomax() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_lomax( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  return ( (x<0.) ? 0. : 1. - pow((C/(x+C)),a) );
} /* end of _unur_cdf_lomax() */

/*---------------------------------------------------------------------------*/

double
_unur_invcdf_lomax( double U, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  double X;

  X = C * ( pow(1-U, -1/a) - 1. );
  return X;
} /* end of _unur_invcdf_lomax() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_lomax( UNUR_DISTR *distr )
{
  DISTR.mode = 0.;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_lomax() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_lomax( UNUR_DISTR *distr )
{
  /* normalization constant */
  NORMCONSTANT = DISTR.a * pow(DISTR.C,DISTR.a);

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }

  /* else */
  DISTR.area = ( _unur_cdf_lomax( DISTR.domain[1],distr) 
		 - _unur_cdf_lomax( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
  
} /* end of _unur_upd_area_lomax() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_lomax( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter a */
  if (a <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"a <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* check parameter C */
  if (n_params > 1 && C <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"C <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form */
  DISTR.a = a;

  /* default parameters */
  DISTR.C = 1.; 

  /* copy optional parameters */
  if (n_params == 2)
    DISTR.C = C;

  /* store total number of parameters */
  DISTR.n_params = 2;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0;               /* left boundary  */
    DISTR.domain[1] = INFINITY;        /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_lomax() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_lomax( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_LOMAX;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  /* DISTR.init = _unur_stdgen_lomax_init; */

  /* functions */
  DISTR.pdf    = _unur_pdf_lomax;    /* pointer to PDF               */
  DISTR.dpdf   = _unur_dpdf_lomax;   /* pointer to derivative of PDF */
  DISTR.cdf    = _unur_cdf_lomax;    /* pointer to CDF               */
  DISTR.invcdf = _unur_invcdf_lomax; /* pointer to inverse CDF       */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
  		 UNUR_DISTR_SET_MODE   | 
  		 UNUR_DISTR_SET_PDFAREA );

  /* set parameters for distribution */
  if (_unur_set_params_lomax(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* normalization constant */
  NORMCONSTANT = DISTR.a * pow(DISTR.C,DISTR.a);

  /* mode and area below p.d.f. */
  DISTR.mode = 0.;
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_lomax;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_lomax; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_lomax; /* funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_lomax() */

/*---------------------------------------------------------------------------*/
#undef a
#undef C
#undef DISTR
/*---------------------------------------------------------------------------*/
