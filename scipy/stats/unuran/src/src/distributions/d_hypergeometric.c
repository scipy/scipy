/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_hypergeometric.c                                           *
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
 *  distr: Hypergeometric distribution  [1; ch.6, p.237]                     *
 *                                                                           *
 *  pmf:       p(k) = (M \choose k) * (N-M \choose n-k) / (N \choose n)      *
 *  domain:    max(0,n-N+M) <= k <= min(n,M)                                 *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  N   >= M                                                          *
 *     1:  M   >= 1                                                          *
 *     2:  n   n >= 1                                                        *
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

static const char distr_name[] = "hypergeometric";

/*---------------------------------------------------------------------------*/
/* parameters */
#define N  params[0]
#define M  params[1]
#define n  params[2]

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.discr
#define LOGNORMCONSTANT (distr->data.discr.norm_constant)

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */
static double _unur_pmf_hypergeometric( int k, const UNUR_DISTR *distr );
#ifdef _unur_SF_cdf_hypergeometric
static double _unur_cdf_hypergeometric( int k, const UNUR_DISTR *distr ); 
#endif
#ifdef _unur_SF_invcdf_hypergeometric
static int    _unur_invcdf_hypergeometric( double u, const UNUR_DISTR *distr ); 
#endif

static int _unur_upd_mode_hypergeometric( UNUR_DISTR *distr );
static int _unur_upd_sum_hypergeometric( UNUR_DISTR *distr );
static int _unur_set_params_hypergeometric( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pmf_hypergeometric(int k, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if ( k < _unur_max(0,(n-N+M-0.5)) || k > _unur_min(n,M)+0.5 ) 
    return 0.;

  else
    return exp( LOGNORMCONSTANT - _unur_SF_ln_factorial(k) - _unur_SF_ln_factorial(M-k) -
                _unur_SF_ln_factorial(n-k) - _unur_SF_ln_factorial(N-M-n+k) );

} /* end of _unur_pmf_hypergeometric() */

/*---------------------------------------------------------------------------*/
#ifdef _unur_SF_cdf_hypergeometric

double
_unur_cdf_hypergeometric(int k, const UNUR_DISTR *distr)
{ 
  const double *params = DISTR.params;

  return _unur_SF_cdf_hypergeometric(k,N,M,n);
} /* end of _unur_cdf_hypergeometric() */

#endif
/*---------------------------------------------------------------------------*/
#ifdef _unur_SF_invcdf_hypergeometric

int
_unur_invcdf_hypergeometric(double u, const UNUR_DISTR *distr)
{ 
  const double *params = DISTR.params;
  double x;

  x = _unur_SF_invcdf_hypergeometric(u,N,M,n);
  return ((x>=INT_MAX) ? INT_MAX : ((int) x));
} /* end of _unur_invcdf_hypergeometric() */

#endif
/*---------------------------------------------------------------------------*/
int
_unur_upd_mode_hypergeometric( UNUR_DISTR *distr )
{
  DISTR.mode = (int) ( (DISTR.n + 1) * (DISTR.M + 1.) / (DISTR.N + 2.) );

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_hypergeometric() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_sum_hypergeometric( UNUR_DISTR *distr )
{
  register double *params = DISTR.params;

  /* log of normalization constant: none */
  LOGNORMCONSTANT = _unur_SF_ln_factorial(M) + _unur_SF_ln_factorial(N-M) + _unur_SF_ln_factorial(n) +
    _unur_SF_ln_factorial(N-n) - _unur_SF_ln_factorial(N);

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.sum = 1.;
    return UNUR_SUCCESS;
  }
  
#ifdef _unur_SF_cdf_hypergeometric
  /* else */
  DISTR.sum = ( _unur_cdf_hypergeometric( DISTR.domain[1],distr) 
		 - _unur_cdf_hypergeometric( DISTR.domain[0]-1,distr) );
  return UNUR_SUCCESS;
#else
  return UNUR_ERR_DISTR_REQUIRED;
#endif

} /* end of _unur_upd_sum_hypergeometric() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_hypergeometric( UNUR_DISTR *distr, const double *params, int n_params )
{
  int nh;

  /* check number of parameters for distribution */
  if (n_params < 3) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameters */
  if (M <= 0. || N <=0. || n <= 0. || n >= N || M >= N ) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"M, N, n must be > 0 and n<N M<N");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form */

  nh = (int)(N+0.5);
  if(fabs(nh-N)>0.001)
    _unur_warning(distr_name,UNUR_ERR_DISTR_DOMAIN,"n was rounded to the closest integer value");
  DISTR.N = nh; 

  nh = (int)(M+0.5);
  if(fabs(nh-M)>0.001)
    _unur_warning(distr_name,UNUR_ERR_DISTR_DOMAIN,"n was rounded to the closest integer value");
  DISTR.M = nh; 

  nh = (int)(n+0.5);
  if(fabs(nh-n)>0.001)
    _unur_warning(distr_name,UNUR_ERR_DISTR_DOMAIN,"n was rounded to the closest integer value");
  DISTR.n = nh; 

  /* default parameters: none */
  /* copy optional parameters: none */

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = (int) (_unur_max(0,(DISTR.n - DISTR.N + DISTR.M + 0.5)));  /* left boundary  */
    DISTR.domain[1] = (int) (_unur_min(DISTR.n, DISTR.M) + 0.5);                 /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_hypergeometric() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_hypergeometric( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_discr_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_HYPERGEOMETRIC;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = _unur_stdgen_hypergeometric_init;
   
  /* functions */
  DISTR.pmf  = _unur_pmf_hypergeometric;   /* pointer to PMF */
#ifdef _unur_SF_cdf_hypergeometric
  DISTR.cdf  = _unur_cdf_hypergeometric;   /* pointer to CDF */
#endif
#ifdef _unur_SF_invcdf_hypergeometric
  DISTR.invcdf = _unur_invcdf_hypergeometric;  /* pointer to inverse CDF */
#endif

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PMFSUM |
		 UNUR_DISTR_SET_MODE );
                
  /* set parameters for distribution */
  if (_unur_set_params_hypergeometric(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant */
  _unur_upd_sum_hypergeometric( distr );

  /* mode and sum over PMF */
  _unur_upd_mode_hypergeometric(distr);
  DISTR.sum = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_hypergeometric;

  /* function for updating derived parameters */
  DISTR.upd_mode = _unur_upd_mode_hypergeometric; /* funct for computing mode */
  DISTR.upd_sum  = _unur_upd_sum_hypergeometric;  /* funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_hypergeometric() */

/*---------------------------------------------------------------------------*/
#undef N
#undef M
#undef n
#undef DISTR
/*---------------------------------------------------------------------------*/
