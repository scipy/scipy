/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_slash.c                                                    *
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
 *  distr: Slash distribution [2; ch.12, p.63]                               *
 *                                                                           *
 *  pdf:       f(x) =  (1 - exp(-x^2/2)) / x^2    for x != 0                 *
 *  pdf:       f(x) =  1 / 2                      for x == 0                 *
 *  domain:    -infinity < x < infinity                                      *
 *  constant:  1 / sqrt(2 * pi)                                              *
 *                                                                           *
 *  parameters: none                                                         *
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

static const char distr_name[] = "slash";

/* parameters */
/* none */

#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)

/* function prototypes                                                       */
static double _unur_pdf_slash( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_slash( double x, const UNUR_DISTR *distr );
/*  static double _unur_cdf_slash( double x, const UNUR_DISTR *distr ); */

static int _unur_upd_mode_slash( UNUR_DISTR *distr );
static int _unur_set_params_slash( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_slash(double x, const UNUR_DISTR *distr)
{
  if (_unur_iszero(x))
    return (0.5 * NORMCONSTANT);
  else
    return ((1. - exp(-x*x/2.)) / (x*x) * NORMCONSTANT);
} /* end of _unur_pdf_slash() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_slash(double x, const UNUR_DISTR *distr ATTRIBUTE__UNUSED)
{ 
  register double xsq = x * x;

  if (_unur_iszero(x))
    return 0.;
  else
    return (NORMCONSTANT * ((-2. + exp(-xsq/2.) * (2. + xsq)) / (xsq * x)));

} /* end of _unur_dpdf_slash() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_slash( UNUR_DISTR *distr )
{
  DISTR.mode = 0.;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_slash() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_slash( UNUR_DISTR *distr, const double *params ATTRIBUTE__UNUSED, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params > 0)
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");

  /* copy parameters for standard form: none */
  /* default parameters: none */
  /* copy optional parameters: none */

  /* store number of parameters */
  DISTR.n_params = 0;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -INFINITY;       /* left boundary  */
    DISTR.domain[1] = INFINITY;        /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_slash() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_slash( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_SLASH;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = _unur_stdgen_slash_init;
   
  /* functions */
  DISTR.pdf  = _unur_pdf_slash;   /* pointer to PDF               */
  DISTR.dpdf = _unur_dpdf_slash;  /* pointer to derivative of PDF */
  /* DISTR.cdf  = _unur_cdf_slash;   pointer to CDF               */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
  		 UNUR_DISTR_SET_MODE   | 
  		 UNUR_DISTR_SET_PDFAREA );
                
  /* set parameters for distribution */
  if (_unur_set_params_slash(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* normalization constant */
  NORMCONSTANT = 1. / (M_SQRT2 * M_SQRTPI);

  /* mode and area below p.d.f. */
  DISTR.mode = 0.;
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_slash;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_slash;   /* funct for computing mode */
  /* DISTR.upd_area  = _unur_upd_area_slash;   funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_slash() */

/*---------------------------------------------------------------------------*/
#undef nu
#undef DISTR
/*---------------------------------------------------------------------------*/
