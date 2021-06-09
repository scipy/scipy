/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_negativebinomial.c                                         *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *   [1] N.L. Johnson, S. Kotz and A.W. Kemp                                 *
 *       Univariate Discrete Distributions,                                  *
 *       2nd edition                                                         *
 *       John Wiley & Sons, Inc., New York, 1992                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  distr: Negative Binomial distribution  [1; ch.5.1, p.200]                *
 *                                                                           *
 *  pmf:       p(k) = (k+r-1 \choose r-1) * p^r * (1-p)^k                    *
 *  domain:    0 <= k < infinity                                             *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  0 < p < 1                                                         *
 *     1:      r > 0                                                         *
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
#include <distr/discr.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "negativebinomial";

/*---------------------------------------------------------------------------*/
/* parameters */
#define p  params[0]
#define r  params[1]

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.discr
#define LOGNORMCONSTANT (distr->data.discr.norm_constant)

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */

static double _unur_pmf_negativebinomial( int k, const UNUR_DISTR *distr );
#ifdef _unur_SF_cdf_negativebinomial
static double _unur_cdf_negativebinomial( int k, const UNUR_DISTR *distr ); 
#endif
#ifdef _unur_SF_invcdf_negativebinomial
static int    _unur_invcdf_negativebinomial( double u, const UNUR_DISTR *distr ); 
#endif

static int _unur_upd_mode_negativebinomial( UNUR_DISTR *distr );
static int _unur_upd_sum_negativebinomial( UNUR_DISTR *distr );
static int _unur_set_params_negativebinomial( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pmf_negativebinomial(int k, const UNUR_DISTR *distr)
{ 
  const double *params = DISTR.params;

  if (k<0)
    return 0.;

  else
    return exp( k*log(1-p) 
		+ _unur_SF_ln_gamma(k+r) - _unur_SF_ln_gamma(k+1.) - LOGNORMCONSTANT ) ;

} /* end of _unur_pmf_negativebinomial() */

/*---------------------------------------------------------------------------*/
#ifdef _unur_SF_cdf_negativebinomial

double
_unur_cdf_negativebinomial(int k, const UNUR_DISTR *distr)
{ 
  const double *params = DISTR.params;

  if (k<0)
    return 0.;

  else
    return _unur_SF_cdf_negativebinomial(k,r,p);

} /* end of _unur_cdf_negativebinomial() */

#endif
/*---------------------------------------------------------------------------*/
#ifdef _unur_SF_invcdf_negativebinomial

int
_unur_invcdf_negativebinomial(double u, const UNUR_DISTR *distr)
{ 
  const double *params = DISTR.params;
  double x;

  x = _unur_SF_invcdf_negativebinomial(u,r,p);
  return ((x>=INT_MAX) ? INT_MAX : ((int) x));
} /* end of _unur_invcdf_negativebinomial() */

#endif
/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_negativebinomial( UNUR_DISTR *distr )
{
  if (DISTR.r > 1.) {
    /* mode = floor( (r-1) * (1-p) / p ) */
    /* (we add a guard against round-off errors */
    DISTR.mode = (int) ((1.+UNUR_EPSILON) * (DISTR.r - 1.) * (1. - DISTR.p) / DISTR.p);
  }
  else { /* r <= 1. */
    DISTR.mode = 0;
  }

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_negativebinomial() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_sum_negativebinomial( UNUR_DISTR *distr )
{
  /* log of normalization constant */
  LOGNORMCONSTANT = - DISTR.r * log(DISTR.p) + _unur_SF_ln_gamma(DISTR.r);

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.sum = 1.;
    return UNUR_SUCCESS;
  }
  
#ifdef _unur_SF_cdf_negativebinomial
  /* else */
  DISTR.sum = ( _unur_cdf_negativebinomial( DISTR.domain[1],distr) 
		 - _unur_cdf_negativebinomial( DISTR.domain[0]-1,distr) );
  return UNUR_SUCCESS;
#else
  return UNUR_ERR_DISTR_REQUIRED;
#endif

} /* end of _unur_upd_sum_negativebinomial() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_negativebinomial( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameters */
  if (p <= 0. || p >= 1. || r <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"p <= 0 || p >= 1 || r <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form */
  DISTR.p = p;
  DISTR.r = r;

  /* default parameters: none */
  /* copy optional parameters: none */

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain: [0, infinity] */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0;           /* left boundary  */
    DISTR.domain[1] = INT_MAX;     /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_negativebinomial() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_negativebinomial( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_discr_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_NEGATIVEBINOMIAL;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  /*    DISTR.init = _unur_stdgen_negativebinomial_init; */
   
  /* functions */
  DISTR.pmf  = _unur_pmf_negativebinomial;   /* pointer to PMF */
#ifdef _unur_SF_cdf_negativebinomial
  DISTR.cdf  = _unur_cdf_negativebinomial;   /* pointer to CDF */
#endif
#ifdef _unur_SF_invcdf_negativebinomial
  DISTR.invcdf = _unur_invcdf_negativebinomial;  /* pointer to inverse CDF */
#endif

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PMFSUM |
		 UNUR_DISTR_SET_MODE );
                
  /* set parameters for distribution */
  if (_unur_set_params_negativebinomial(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant */
  _unur_upd_sum_negativebinomial( distr );

  /* mode and sum over PMF */
  _unur_upd_mode_negativebinomial(distr);
  DISTR.sum = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_negativebinomial;

  /* function for updating derived parameters */
  DISTR.upd_mode = _unur_upd_mode_negativebinomial; /* funct for computing mode */
  DISTR.upd_sum  = _unur_upd_sum_negativebinomial;  /* funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_negativebinomial() */

/*---------------------------------------------------------------------------*/
#undef p
#undef r
#undef DISTR
/*---------------------------------------------------------------------------*/
