/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_poisson.c                                                  *
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
 *  distr: Poisson distribution  [1; ch.4, p.151]                            *
 *                                                                           *
 *  pmf:       p(k) = theta^k / k!                                           *
 *  domain:    0 <= k < infinity                                             *
 *  constant:  exp( -theta )                                                 *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  theta > 0  ... shape                                              *
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

static const char distr_name[] = "poisson";

/* parameters */
#define theta  params[0]

#define DISTR distr->data.discr
/*  #define LOGNORMCONSTANT (distr->data.discr.norm_constant) */

/*---------------------------------------------------------------------------*/

/* function prototypes                                                       */
static double _unur_pmf_poisson( int k, const UNUR_DISTR *distr );
static double _unur_cdf_poisson( int k, const UNUR_DISTR *distr );      
#ifdef _unur_SF_invcdf_binomial
static int    _unur_invcdf_poisson( double u, const UNUR_DISTR *distr ); 
#endif

static int _unur_upd_mode_poisson( UNUR_DISTR *distr );
static int _unur_upd_sum_poisson( UNUR_DISTR *distr );
static int _unur_set_params_poisson( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pmf_poisson(int k, const UNUR_DISTR *distr)
{ 
  if (k>=0)
    return exp( -DISTR.theta + k * log(DISTR.theta) - _unur_SF_ln_factorial(k) );
  else
    return 0.;
} /* end of _unur_pmf_poisson() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_poisson(int k, const UNUR_DISTR *distr)
{ 
  if (k>=0)
    return (1.-_unur_SF_incomplete_gamma(DISTR.theta,k+1.));
  else
    return 0.;
} /* end of _unur_cdf_poisson() */

/*---------------------------------------------------------------------------*/
#ifdef _unur_SF_invcdf_poisson

int
_unur_invcdf_poisson(double u, const UNUR_DISTR *distr)
{ 
  const double *params = DISTR.params;
  double x;

  x = _unur_SF_invcdf_poisson(u,theta);
  return ((x>=INT_MAX) ? INT_MAX : ((int) x));
} /* end of _unur_invcdf_poisson() */

#endif
/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_poisson( UNUR_DISTR *distr )
{
  DISTR.mode = (int) DISTR.theta;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_poisson() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_sum_poisson( UNUR_DISTR *distr )
{
  /* log normalization constant: none */

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.sum = 1.;
    return UNUR_SUCCESS;
  }
  
  /* else */
  DISTR.sum = ( _unur_cdf_poisson( DISTR.domain[1],distr) 
		 - _unur_cdf_poisson( DISTR.domain[0]-1,distr) );
  return UNUR_SUCCESS;
} /* end of _unur_upd_sum_poisson() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_poisson( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 1) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 1; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter theta */
  if (theta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"theta <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form */
  DISTR.theta = theta;

  /* default parameters: none */
  /* copy optional parameters: none */

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain: [0, inifinty] */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0;           /* left boundary  */
    DISTR.domain[1] = INT_MAX;     /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_poisson() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_poisson( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_discr_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_POISSON;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = _unur_stdgen_poisson_init;
   
  /* functions */
  DISTR.pmf  = _unur_pmf_poisson;   /* pointer to PMF */
  DISTR.cdf  = _unur_cdf_poisson;   /* pointer to CDF */
#ifdef _unur_SF_invcdf_poisson
  DISTR.invcdf = _unur_invcdf_poisson;  /* pointer to inverse CDF */
#endif

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PMFSUM |
		 UNUR_DISTR_SET_MODE );
                
  /* set parameters for distribution */
  if (_unur_set_params_poisson(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* log normalization constant: none */

  /* mode and sum over PMF */
  DISTR.mode = (int) DISTR.theta;
  DISTR.sum = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_poisson;

  /* function for updating derived parameters */
  DISTR.upd_mode = _unur_upd_mode_poisson; /* funct for computing mode */
  DISTR.upd_sum  = _unur_upd_sum_poisson;  /* funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_poisson() */

/*---------------------------------------------------------------------------*/
#undef theta
#undef DISTR
/*---------------------------------------------------------------------------*/
