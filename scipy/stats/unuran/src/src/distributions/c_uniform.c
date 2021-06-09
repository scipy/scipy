/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_uniform.c                                                  *
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
 *  distr: Uniform distribution [3; ch.26, p.276]                            *
 *                                                                           *
 *  pdf:       f(x) = 1 / (b-a)                                              *
 *  domain:    a <= x <= b                                                   *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  a      ... location                                               *
 *     1:  b (>a) ... location                                               *
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
static const char distr_name[] = "uniform";

/* parameters */
#define a  params[0]
#define b  params[1]

#define DISTR distr->data.cont
/* #define NORMCONSTANT (distr->data.cont.norm_constant) */

/* function prototypes                                                       */
static double _unur_pdf_uniform( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_uniform( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_uniform( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_uniform( double x, const UNUR_DISTR *distr );
static double _unur_cdf_uniform( double x, const UNUR_DISTR *distr );
static double _unur_invcdf_uniform( double u, const UNUR_DISTR *distr );

static int _unur_upd_mode_uniform( UNUR_DISTR *distr );
static int _unur_upd_area_uniform( UNUR_DISTR *distr );
static int _unur_set_params_uniform( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_uniform( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;

  if (x < a || x > b)
    return 0.;
  /* else */
  return 1./(b-a);

  /*    return ((x < a || x > b) ? 0. : 1./(b-a)); */
} /* end of _unur_pdf_uniform() */

/*---------------------------------------------------------------------------*/

double
_unur_logpdf_uniform( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;

  if (x < a || x > b)
    return -INFINITY;
  /* else */
  return -log(b-a);

} /* end of _unur_logpdf_uniform() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_uniform( double x ATTRIBUTE__UNUSED, const UNUR_DISTR *distr ATTRIBUTE__UNUSED )
{ 
  return 0.;
} /* end of _unur_dpdf_uniform() */

/*---------------------------------------------------------------------------*/

double
_unur_dlogpdf_uniform( double x ATTRIBUTE__UNUSED, const UNUR_DISTR *distr ATTRIBUTE__UNUSED )
{ 
  return 0.;
} /* end of _unur_dlogpdf_uniform() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_uniform( double x, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;

  /* standardize */
  x = (x-a) / (b-a);

  if (x<=0.)
    return 0.;
  if (x>=1.)
    return 1.;
  return x;

} /* end of _unur_cdf_uniform() */

/*---------------------------------------------------------------------------*/

double
_unur_invcdf_uniform( double U, const UNUR_DISTR *distr )
{ 
  register const double *params = DISTR.params;
  double X;

  X = (DISTR.n_params==0) ? U : a + U * (b - a);
  return X;
} /* end of _unur_invcdf_uniform() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_uniform( UNUR_DISTR *distr )
{
  DISTR.mode = (DISTR.a + DISTR.b) / 2.;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_uniform() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_uniform( UNUR_DISTR *distr )
{
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }

  /* else */
  DISTR.area = ( _unur_cdf_uniform( DISTR.domain[1],distr) 
		 - _unur_cdf_uniform( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
  
} /* end of _unur_upd_area_uniform() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_uniform( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 0) n_params = 0;
  if (n_params == 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  if (n_params > 0)
    CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameters a and b */
  if (n_params == 2 && (a >= b)) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"a >= b");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form: none */

  /* default parameters */
  DISTR.a = 0.;
  DISTR.b = 1.;

  /* copy optional parameters */
  if (n_params == 2) {
    DISTR.a = a;
    DISTR.b = b;
  }

  /* store total number of parameters */
  DISTR.n_params = 2;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = DISTR.a;      /* left boundary  */
    DISTR.domain[1] = DISTR.b;      /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_uniform() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_uniform( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_UNIFORM;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  /* DISTR.init = _unur_stdgen_uniform_init; */

  /* functions */
  DISTR.pdf     = _unur_pdf_uniform;     /* pointer to PDF               */
  DISTR.logpdf  = _unur_logpdf_uniform;  /* pointer to logPDF            */
  DISTR.dpdf    = _unur_dpdf_uniform;    /* pointer to derivative of PDF */
  DISTR.dlogpdf = _unur_dlogpdf_uniform; /* pointer to deriv. of logPDF  */
  DISTR.cdf     = _unur_cdf_uniform;     /* pointer to CDF               */
  DISTR.invcdf  = _unur_invcdf_uniform;  /* pointer to inverse CDF       */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA );
                
  /* set parameters for distribution */
  if (_unur_set_params_uniform(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* normalization constant: none */

  /* mode and area below p.d.f. */
  DISTR.mode = (DISTR.a + DISTR.b) / 2.;
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_uniform;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_uniform; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_uniform; /* funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_uniform() */

/*---------------------------------------------------------------------------*/
#undef a
#undef b
#undef DISTR
/*---------------------------------------------------------------------------*/
