/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_binomial.c                                                 *
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
 *  distr: Binomial distribution  [1; ch.3, p.105]                           *
 *                                                                           *
 *  pmf:       p(k) = (n \choose k) * p^k * (1-p)^(n-k)                      *
 *  domain:    0 <= k <= n                                                   *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  n >= 1                                                            *
 *     1:  0 < p < 1                                                         *
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

static const char distr_name[] = "binomial";

/*---------------------------------------------------------------------------*/
/* parameters */
#define n  params[0]
#define p  params[1]

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.discr
/* #define LOGNORMCONSTANT (distr->data.discr.norm_constant) */

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */
static double _unur_pmf_binomial( int k, const UNUR_DISTR *distr );
static double _unur_cdf_binomial( int k, const UNUR_DISTR *distr ); 
#ifdef _unur_SF_invcdf_binomial
static int    _unur_invcdf_binomial( double u, const UNUR_DISTR *distr ); 
#endif

static int _unur_upd_mode_binomial( UNUR_DISTR *distr );
static int _unur_upd_sum_binomial( UNUR_DISTR *distr );
static int _unur_set_params_binomial( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pmf_binomial(int k, const UNUR_DISTR *distr)
{ 
  const double *params = DISTR.params;

  if ( k<0 || k>(n+0.5) )
    return 0.;

  else
    return exp( k * log(p) + (n-k) * log(1.-p) +
		_unur_SF_ln_factorial(n) - _unur_SF_ln_factorial(k) - _unur_SF_ln_factorial(n-k) ) ;

} /* end of _unur_pmf_binomial() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_binomial(int k, const UNUR_DISTR *distr)
{ 
  const double *params = DISTR.params;

  if (k<0)
    return 0.;

  if (k==0)
    return exp(n*(log(1.-p)));

  if(k>(n-0.5))
    return 1.;

  /* else */
  return(_unur_SF_incomplete_beta(1.-p, n-k, k+1.));

} /* end of _unur_cdf_binomial() */

/*---------------------------------------------------------------------------*/
#ifdef _unur_SF_invcdf_binomial

int
_unur_invcdf_binomial(double u, const UNUR_DISTR *distr)
{ 
  const double *params = DISTR.params;
  double x;

  x = _unur_SF_invcdf_binomial(u,n,p);
  return ((x>=INT_MAX) ? INT_MAX : ((int) x));
} /* end of _unur_invcdf_binomial() */

#endif
/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_binomial( UNUR_DISTR *distr )
{
  DISTR.mode = (int) ((DISTR.n + 1) * DISTR.p);

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_binomial() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_sum_binomial( UNUR_DISTR *distr )
{
  /* log of normalization constant: none */

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.sum = 1.;
    return UNUR_SUCCESS;
  }

  /* else */
  DISTR.sum = ( _unur_cdf_binomial( DISTR.domain[1],distr) 
		 - _unur_cdf_binomial( DISTR.domain[0]-1,distr) );
  return UNUR_SUCCESS;

} /* end of _unur_upd_sum_binomial() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_binomial( UNUR_DISTR *distr, const double *params, int n_params )
{
  int nh;

  /* check number of parameters for distribution */
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameters */
  if (p <= 0. || p >= 1. || n <= 0.) { 
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"p <= 0 || p >= 1 || n <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form */

  /* round parameter n */
  nh = (int)(n+0.5);
  if(fabs(nh-n)>0.001)
    _unur_warning(distr_name,UNUR_ERR_DISTR_DOMAIN,"n was rounded to the closest integer value");
  DISTR.n = nh;
  DISTR.p = p;

  /* default parameters: none */
  /* copy optional parameters: none */

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain: [0, n] */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0;        /* left boundary  */
    DISTR.domain[1] = nh;       /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_binomial() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_binomial( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_discr_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_BINOMIAL;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = _unur_stdgen_binomial_init;
   
  /* functions */
  DISTR.pmf  = _unur_pmf_binomial;   /* pointer to PMF */
  DISTR.cdf  = _unur_cdf_binomial;   /* pointer to CDF */
#ifdef _unur_SF_invcdf_binomial
  DISTR.invcdf = _unur_invcdf_binomial;  /* pointer to inverse CDF */
#endif

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PMFSUM |
		 UNUR_DISTR_SET_MODE );

  /* set parameters for distribution */
  if (_unur_set_params_binomial(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant: none */

  /* mode and sum over PMF */
  DISTR.mode = (int) ((DISTR.n + 1) * DISTR.p);
  DISTR.sum = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_binomial;

  /* function for updating derived parameters */
  DISTR.upd_mode = _unur_upd_mode_binomial; /* funct for computing mode */
  DISTR.upd_sum  = _unur_upd_sum_binomial;  /* funct for computing area */
                
  /* return pointer to object */
  return distr;

} /* end of unur_distr_binomial() */

/*---------------------------------------------------------------------------*/
#undef p
#undef r
#undef DISTR
/*---------------------------------------------------------------------------*/
