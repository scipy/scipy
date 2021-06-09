/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_hyperbolic.c                                               *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  distr: Hyperbolic distribution                                           *
 *                                                                           *
 *  pdf:   f(x) = exp( -alpha * sqrt(delta^2 + (x - mu)^2) + beta*(x-mu) )   *
 *  domain:    infinity < x < infinity                                       *
 *  constant:  gamma / (2 * alpha * delta * K_1(delta * gamma)               *
 *          [gamma = sqrt(alpha^2 - beta^2)                        ]         *
 *          [K_theta(.) ... modified Bessel function of second (third) kind] *
 *                                                                           *
 *  parameters: 4                                                            *
 *     0 : alpha >|beta|     ... shape (tail)                                *
 *     1 : beta              ... shape (asymmetry)                           *
 *     2 : delta > 0         ... scale                                       *
 *     3 : mu                ... location                                    *
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
#include <distr/cont.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "hyperbolic";

/* parameters */
#define alpha  params[0]    /* shape (tail) */
#define beta   params[1]    /* shape (asymmetry) */
#define delta  params[2]    /* scale */
#define mu     params[3]    /* location */

#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)

/* function prototypes                                                       */
static double _unur_pdf_hyperbolic( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_hyperbolic( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_hyperbolic( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_hyperbolic( double x, const UNUR_DISTR *distr );
/* static double _unur_cdf_hyperbolic( double x, const UNUR_DISTR *distr ); */

static int _unur_upd_mode_hyperbolic( UNUR_DISTR *distr );
static double _unur_normconstant_hyperbolic( const double *params, int n_params );
static int _unur_set_params_hyperbolic( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_hyperbolic(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  return (NORMCONSTANT * exp(-alpha * sqrt(delta*delta + (x-mu)*(x-mu)) + beta*(x-mu) ) );
} /* end of _unur_pdf_hyperbolic() */

/*---------------------------------------------------------------------------*/

double
_unur_logpdf_hyperbolic(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  return (-alpha * sqrt(delta*delta + (x-mu)*(x-mu)) + beta*(x-mu) + log(NORMCONSTANT) );
} /* end of _unur_logpdf_hyperbolic() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_hyperbolic(double x, const UNUR_DISTR *distr)
{ 
  return (NORMCONSTANT * _unur_pdf_hyperbolic(x,distr) * _unur_dlogpdf_hyperbolic(x,distr));
} /* end of _unur_dpdf_hyperbolic() */

/*---------------------------------------------------------------------------*/

double
_unur_dlogpdf_hyperbolic(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  return (beta - (alpha*(x-mu))/sqrt(delta*delta + (x-mu)*(x-mu)) + log(NORMCONSTANT));
} /* end of _unur_dlogpdf_hyperbolic() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_hyperbolic( UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;

  DISTR.mode =
    mu + delta*beta / sqrt(alpha*alpha - beta*beta);

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_hyperbolic() */

/*---------------------------------------------------------------------------*/

double
_unur_normconstant_hyperbolic(const double *params ATTRIBUTE__UNUSED, int n_params ATTRIBUTE__UNUSED)
{ 
#ifdef _unur_SF_bessel_k
  double gamm = sqrt(alpha*alpha-beta*beta);

  /*
   *  constant:  gamma / (2 * alpha * delta * K_1(delta * gamma)               *
   *             [gamma = sqrt(alpha^2 - beta^2)                        ]      *
   *             [K_theta(.) ... modified Bessel function of second kind]      *
   */

  return ( gamm / ( 2 * alpha * delta * _unur_SF_bessel_k(delta*gamm, 1) ) );
#else
  return 1.;
#endif
} /* end of _unur_normconstant_hyperbolic() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_hyperbolic( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 4) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 4) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 4; }
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

  /* default parameters: none */

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -INFINITY;   /* left boundary  */
    DISTR.domain[1] = INFINITY;    /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_hyperbolic() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_hyperbolic( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_HYPERBOLIC;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  /* DISTR.init = _unur_stdgen_hyperbolic_init; */
   
  /* functions */
  DISTR.pdf     = _unur_pdf_hyperbolic;     /* pointer to PDF                  */
  DISTR.logpdf  = _unur_logpdf_hyperbolic;  /* pointer to logPDF               */
  DISTR.dpdf    = _unur_dpdf_hyperbolic;    /* pointer to derivative of PDF    */
  DISTR.dlogpdf = _unur_dlogpdf_hyperbolic; /* pointer to derivative of logPDF */
  DISTR.cdf  = NULL;                 /* _unur_cdf_hyperbolic; pointer to CDF   */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   );
		 /* UNUR_DISTR_SET_PDFAREA ); */
                
  /* set parameters for distribution */
  if (_unur_set_params_hyperbolic(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* normalization constant */
  NORMCONSTANT = _unur_normconstant_hyperbolic(params,n_params);

  /* mode and area below p.d.f. */
  _unur_upd_mode_hyperbolic(distr);
  /*    DISTR.area = ? */

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_hyperbolic;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_hyperbolic; /* funct for computing mode */
  /* DISTR.upd_area  = _unur_upd_area_hyperbolic; funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_hyperbolic() */

/*---------------------------------------------------------------------------*/
#undef mu
#undef alpha
#undef beta
#undef delta
#undef DISTR
/*---------------------------------------------------------------------------*/
