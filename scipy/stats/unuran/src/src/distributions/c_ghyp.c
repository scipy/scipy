/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_ghyp.c                                                     *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  distr: Generalized hyperbolic distribution                               *
 *                                                                           *
 *  pdf:   f(x) = (delta^2+(x-mu)^2)^(1/2*(lambda-1/2)) * exp(beta*(x-mu))   *
 *                  * K_{lambda-1/2}(alpha*sqrt(delta^2+(x-mu)^2))           * 
 *                                                                           *
 *  domain:   infinity < x < infinity                                        *
 *                                                                           *
 *  constant: ( (gamma/delta)^lambda )                                       *
 *            / ( sqrt(2*pi) * alpha^(lambda-1/2) * K_{lambda}(delta*gamma) )
 *                                                                           *
 *             [gamma = sqrt(alpha^2 - beta^2)                        ]      *
 *             [K_theta(.) ... modified Bessel function of second kind]      *
 *                                                                           *
 *  parameters: 4                                                            *
 *     0 : lambda        ... shape                                           *
 *     1 : alpha >|beta| ... shape                                           *
 *     2 : beta          ... shape (asymmetry)                               *
 *     3 : delta > 0     ... scale                                           *
 *     4 : mu            ... location                                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   We use the Rmath library for computing the Bessel function K_n.         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2009-2010 Wolfgang Hoermann and Josef Leydold             *
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

static const char distr_name[] = "ghyp";

/* parameters */
#define lambda  params[0]    /* shape */
#define alpha   params[1]    /* shape */
#define beta    params[2]    /* shape (asymmetry) */
#define delta   params[3]    /* scale */
#define mu      params[4]    /* location */

#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)

/* function prototypes                                                       */
#ifdef _unur_SF_bessel_k
static double _unur_pdf_ghyp( double x, const UNUR_DISTR *distr );
/* static double _unur_dpdf_ghyp( double x, const UNUR_DISTR *distr ); */
/* static double _unur_cdf_ghyp( double x, const UNUR_DISTR *distr ); */
#endif

static int _unur_upd_center_ghyp( UNUR_DISTR *distr );
static double _unur_normconstant_ghyp( const double *params, int n_params );
static int _unur_set_params_ghyp( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

#ifdef _unur_SF_bessel_k
double
_unur_pdf_ghyp(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;
  double tmp = delta*delta + (x-mu)*(x-mu);

  /* f(x) = (delta^2+(x-mu)^2)^(1/2*(lambda-1/2)) * exp(beta*(x-mu))    */
  /*          * K_{lambda-1/2}(alpha*sqrt(delta^2+(x-mu)^2))            */

  return ( NORMCONSTANT 
	   * pow( tmp, 0.5*lambda-0.25 ) 
	   * exp(beta*(x-mu))
	   * _unur_SF_bessel_k( alpha * sqrt(tmp), lambda-0.5 ) );

} /* end of _unur_pdf_ghyp() */
#endif

/*---------------------------------------------------------------------------*/

int
_unur_upd_center_ghyp( UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;

  /* we simply use parameter 'mu' */
  DISTR.center = mu;

  /* an alternative approach would be the mean of the distribution:          */
  /* double gamma = sqrt(alpha*alpha-beta*beta);                             */
  /* DISTR.center = ( mu                                                     */
  /* 		   + ( delta*beta * bessel_k( delta*gamma, lambda+1, 1) )    */
  /* 		   / ( gamma * bessel_k( delta*gamma, lambda, 1) ) );        */

  /* center must be in domain */
  if (DISTR.center < DISTR.domain[0]) 
    DISTR.center = DISTR.domain[0];
  else if (DISTR.center > DISTR.domain[1]) 
    DISTR.center = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_center_ghyp() */

/*---------------------------------------------------------------------------*/

double
_unur_normconstant_ghyp(const double *params ATTRIBUTE__UNUSED, int n_params ATTRIBUTE__UNUSED)
{ 
#ifdef _unur_SF_bessel_k
  double gamm = sqrt(alpha*alpha-beta*beta);

  /* ( (gamma/delta)^lambda ) / ( sqrt(2*pi) * alpha^(lambda-1/2) * K_{lambda}(delta*gamma) ) */

  return ( pow(gamm/delta, lambda ) 
	   / ( (M_SQRTPI*M_SQRT2) * pow(alpha, lambda-0.5)
	       * _unur_SF_bessel_k( delta*gamm, lambda ) ) );
#else
  return 1.;
#endif
} /* end of _unur_normconstant_ghyp() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_ghyp( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 5) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 5) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 5; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter omega */
  if (delta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"delta <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  if (alpha <= fabs(beta)) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"alpha <= |beta|");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form */
  DISTR.mu = mu;
  DISTR.alpha = alpha;
  DISTR.beta = beta;
  DISTR.delta = delta;
  DISTR.lambda = lambda;

  /* default parameters: none */

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -INFINITY;   /* left boundary  */
    DISTR.domain[1] = INFINITY;    /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_ghyp() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_ghyp( const double *params, int n_params)
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_GHYP;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  /* DISTR.init = _unur_stdgen_ghyp_init; */
   
  /* functions */
#ifdef _unur_SF_bessel_k
  DISTR.pdf     = _unur_pdf_ghyp;     /* pointer to PDF                  */
#endif

  /* indicate which parameters are set */
#ifdef _unur_SF_bessel_k
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_CENTER );
		 /* UNUR_DISTR_SET_PDFAREA ); */
#else
  distr->set = ( UNUR_DISTR_SET_DOMAIN | UNUR_DISTR_SET_STDDOMAIN );
#endif
                
  /* set parameters for distribution */
  if (_unur_set_params_ghyp(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* normalization constant */
  NORMCONSTANT = _unur_normconstant_ghyp(params,n_params);

  /* we need the center of the distribution */
  if (_unur_upd_center_ghyp(distr)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* mode and area below p.d.f. */
  /* DISTR.mode = ? */
  /* DISTR.area = ? */

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_ghyp;

  /* function for updating derived parameters */
  /* DISTR.upd_mode  = _unur_upd_mode_ghyp; /\* funct for computing mode *\/ */
  /* DISTR.upd_area  = _unur_upd_area_ghyp; funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_ghyp() */

/*---------------------------------------------------------------------------*/
#undef mu
#undef alpha
#undef beta
#undef delta
#undef lambda
#undef DISTR
/*---------------------------------------------------------------------------*/
