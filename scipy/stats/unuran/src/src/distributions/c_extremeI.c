/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_extremeI.c                                                 *
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
 *  distr: Extreme value type I distribution  [3; ch.22, p.2]                *
 *  (also Gumbel-type distribution)                                          *
 *                                                                           *
 *                                                                           *
 *  Type I (also Gumbel-type distribution)                                   *
 *                                                                           *
 *  cdf:       F(x) = exp( -exp( -(x-zeta)/theta ) )                         *
 *  pdf:       f(x) = exp( -exp( -(x-zeta)/theta ) - (x-zeta)/theta )        *               
 *  domain:    -infinity < x <infinity                                       *
 *  constant:  1 / theta                                                     *
 *                                                                           *
 *  parameters: 2                                                            *
 *     0:  zeta      (0)  ... location                                       *
 *     1:  theta >0  (1)  ... scale                                          *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 * Standard Form:                                                            *
 *                                                                           *
 *  cdf:       F(x) = exp( -exp(-x) )                                        *
 *  pdf:       f(x) = exp( -exp(-x) - x )                                    *
 *  domain:    -infinity < x <infinity                                       *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters: none                                                         *
 *                                                                           *
 *     0:  zeta  = 0                                                         *
 *     1:  theta = 1                                                         *
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

static const char distr_name[] = "extremeI";

/* parameters */
#define zeta   params[0]    /* location */
#define theta  params[1]    /* scale */

#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)

/* function prototypes                                                       */
static double _unur_pdf_extremeI( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_extremeI( double x, const UNUR_DISTR *distr );
static double _unur_cdf_extremeI( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_extremeI( double u, const UNUR_DISTR *distr );

static int _unur_upd_mode_extremeI( UNUR_DISTR *distr );
static int _unur_upd_area_extremeI( UNUR_DISTR *distr );
static int _unur_set_params_extremeI( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_extremeI( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;

  if (DISTR.n_params > 0)
    /* standardize */
    x = (x - zeta) / theta;

  /* standard form */

  return ( exp( -exp(-x) - x - LOGNORMCONSTANT) );

} /* end of _unur_pdf_extremeI() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_extremeI( double x, const UNUR_DISTR *distr )
{ 
  register double factor = 1.;
  register double expx;
  register const double *params = DISTR.params;

  if (DISTR.n_params > 0) {
    /* standardize */
    factor = 1. / (theta * theta);
    x = (x - zeta) / theta;
  }

  /* standard form */

  expx = exp(-x);
  return ( exp( -expx - x ) * (expx - 1.) * factor );

} /* end of unur_dpdf_extremeI() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_extremeI( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;

  if (DISTR.n_params > 0)
    /* standardize */
    x = (x - zeta) / theta;

  /* standard form */

  return exp( -exp( -x) );

} /* end of _unur_cdf_extremeI() */

/*---------------------------------------------------------------------------*/

double
_unur_invcdf_extremeI( double U, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  double X;

  X = -log( -log(U) );
  return ((DISTR.n_params==0) ? X : zeta + theta * X );
} /* end of _unur_invcdf_extremeI() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_extremeI( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.zeta;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_extremeI() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_extremeI( UNUR_DISTR *distr )
{
  /* log of normalization constant */
  LOGNORMCONSTANT = log(DISTR.theta);


  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }

  /* else */
  DISTR.area = ( _unur_cdf_extremeI( DISTR.domain[1],distr) 
		 - _unur_cdf_extremeI( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
  
} /* end of _unur_upd_area_extremeI() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_extremeI( UNUR_DISTR *distr, const double *params, int n_params )
{

  /* check number of parameters for distribution */
  if (n_params < 0) n_params = 0;
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  if (n_params > 0)
    CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter theta */
  if (n_params == 2 && theta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"theta <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form: none */

  /* default parameters */
  DISTR.zeta  = 0.;
  DISTR.theta = 1.;

  /* copy optional parameters */
  switch (n_params) {
  case 2:
    DISTR.theta = theta;
  case 1:
    DISTR.zeta = zeta;
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
} /* end of _unur_set_params_extremeI() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_extremeI( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_EXTREME_I;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  /* DISTR.init = _unur_stdgen_extremeI_init; */

  /* functions */
  DISTR.pdf    = _unur_pdf_extremeI;    /* pointer to PDF                */
  DISTR.dpdf   = _unur_dpdf_extremeI;   /* pointer to derivative of PDF  */
  DISTR.cdf    = _unur_cdf_extremeI;    /* pointer to CDF                */
  DISTR.invcdf = _unur_invcdf_extremeI; /* pointer to inverse CDF        */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );
                
  /* set parameters for distribution */
  if (_unur_set_params_extremeI(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant */
  LOGNORMCONSTANT = log(DISTR.theta);

  /* domain */
  DISTR.domain[0] = -INFINITY;       /* left boundary  */
  DISTR.domain[1] = INFINITY;        /* right boundary */

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.zeta;
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_extremeI;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_extremeI; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_extremeI; /* funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_extremeI() */

/*---------------------------------------------------------------------------*/
#undef c    
#undef alpha
#undef zeta 
#undef DISTR
/*---------------------------------------------------------------------------*/
