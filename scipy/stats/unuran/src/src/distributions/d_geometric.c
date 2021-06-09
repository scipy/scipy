/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_geometric.c                                                *
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
 *  distr: Geometric distribution  [1; ch.5.2, p.201]                        *
 *                                                                           *
 *  pmf:       p(k) = p * (1-p)^k                                            *
 *  domain:    0 <= k < infinity                                             *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  0 < p < 1   ... shape                                             *
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
#include <distr/discr.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "geometric";

/* parameters */
#define p  params[0]

#define DISTR distr->data.discr
/* #define NORMCONSTANT (distr->data.discr.norm_constant) */

/*---------------------------------------------------------------------------*/

/* function prototypes                                                       */
static double _unur_pmf_geometric( int k, const UNUR_DISTR *distr );
static double _unur_cdf_geometric( int k, const UNUR_DISTR *distr ); 
static int    _unur_invcdf_geometric( double u, const UNUR_DISTR *distr ); 

static int _unur_upd_mode_geometric( UNUR_DISTR *distr );
static int _unur_upd_sum_geometric( UNUR_DISTR *distr );
static int _unur_set_params_geometric( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pmf_geometric(int k, const UNUR_DISTR *distr)
{ 
  return ((k<0) ? 0. : DISTR.p * pow( 1. - DISTR.p, (double)k ));
} /* end of _unur_pmf_geometric() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_geometric(int k, const UNUR_DISTR *distr)
{ 
  return ((k<0) ? 0. : (1. - pow(1. - DISTR.p, k+1.)) );
} /* end of _unur_cdf_geometric() */

/*---------------------------------------------------------------------------*/

int
_unur_invcdf_geometric(double u, const UNUR_DISTR *distr)
{ 
  double x;

  if (_unur_isone(DISTR.p))
    return 0;
  
  /* else */
  x = ceil(log1p(-u) / log1p(-DISTR.p) - 1.);
  
  return ((x>=INT_MAX) ? INT_MAX : ((int) x));
} /* end of _unur_invcdf_geometric() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_geometric( UNUR_DISTR *distr )
{
  DISTR.mode = 0;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0] || DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = (DISTR.domain[0]<0) ? 0 : DISTR.domain[0];

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_geometric() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_sum_geometric( UNUR_DISTR *distr )
{
  /* normalization constant: none */
  
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.sum = 1.;
    return UNUR_SUCCESS;
  }
  
  /* else */
  DISTR.sum = ( _unur_cdf_geometric( DISTR.domain[1],distr) 
		 - _unur_cdf_geometric( DISTR.domain[0]-1,distr) );
  return UNUR_SUCCESS;

} /* end of _unur_upd_sum_geometric() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_geometric( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 1) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 1; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter p */
  if (p <= 0. || p >= 1.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"p <= 0 || p >= 1");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form */
  DISTR.p = p;

  /* default parameters: none */
  /* copy optional parameters: none */

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain: [0,inifinity] */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0;           /* left boundary  */
    DISTR.domain[1] = INT_MAX;     /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_geometric() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_geometric( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_discr_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_GEOMETRIC;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  /* DISTR.init = _unur_stdgen_geometric_init; */
   
  /* functions */
  DISTR.pmf     = _unur_pmf_geometric;    /* pointer to PMF */
  DISTR.cdf     = _unur_cdf_geometric;    /* pointer to CDF */
  DISTR.invcdf  = _unur_invcdf_geometric; /* pointer to invsere of CDF */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE | 
		 UNUR_DISTR_SET_PMFSUM );
                
  /* set parameters for distribution */
  if (_unur_set_params_geometric(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant: none */

  /* mode and sum over PMF */
  DISTR.mode = 0;
  DISTR.sum = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_geometric;

  /* function for updating derived parameters */
  DISTR.upd_mode = _unur_upd_mode_geometric; /* funct for computing mode */
  DISTR.upd_sum  = _unur_upd_sum_geometric; /* funct for computing sum */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_geometric() */

/*---------------------------------------------------------------------------*/
#undef p
#undef DISTR
/*---------------------------------------------------------------------------*/
