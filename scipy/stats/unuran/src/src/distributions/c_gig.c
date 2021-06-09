/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_gig.c                                                      *
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
 *  distr: Generalized Inverse Gaussian (GIG) distribution [2; ch.15, p.284] *
 *                                                                           *
 *  pdf:       f(x) = x^(theta-1) * exp( -1/2 * omega * (x/eta + eta/x))     *
 *  domain:    0 < x < infinity                                              *
 *  constant:  1 / (2 * eta^theta K_theta(omega))                            *
 *             [K_theta(.) ... modified Bessel function of third kind]       *
 * K_theta(x) = 1/2 * int_-inf^inf  cosh(theta*u) * exp(-x*cosh(u)) du       *
 *                                                                           *
 *                              inf                                          *
 *                               -                                           *
 *                          1   |                                            *
 *         K_theta(omega) = -   | cosh(theta*u) * exp(-omega*cosh(u)) du     *
 *                          2   |                                            *
 *                             -                                             *
 *                            -inf                                           *
 *                                                                           *
 *                                                                           *
 *  parameters: 3                                                            *
 *     0:  theta            ... shape                                        *
 *     1:  omega > 0        ... scale                                        *
 *     2:  eta   > 0   (1)  ... shape                                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Reparametrization from "gig2":                                           *
 *     omega = sqrt(psi*chi)                                                 *
 *     eta   = sqrt(chi/psi)                                                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Standard form:                                                           *
 *                                                                           *
 *  pdf:       f(x) = x^(theta-1) * exp( -1/2 * omega * (x + 1/x))           *
 *  domain:    0 < x < infinity                                              *
 *  constant:  2 * K_theta(omega)                                            *
 *                                                                           *
 *  parameters: 2                                                            *
 *     0:  theta       ... shape                                             *
 *     1:  omega > 0   ... scale                                             *
 *                                                                           *
 *     2:  eta   = 1                                                         *
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

static const char distr_name[] = "gig";

/* parameters */
#define theta  params[0]    /* shape */
#define omega  params[1]    /* scale */
#define eta    params[2]    /* shape */

#define DISTR distr->data.cont
/* #define NORMCONSTANT (distr->data.cont.norm_constant) */

/* function prototypes                                                       */
static double _unur_pdf_gig( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_gig( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_gig( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_gig( double x, const UNUR_DISTR *distr );
/* static double _unur_cdf_gig( double x, const UNUR_DISTR *distr ); */

static int _unur_upd_mode_gig( UNUR_DISTR *distr );
static int _unur_set_params_gig( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_gig(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (x <= 0.)
    /* out of support */
    return 0.;

  return (exp( (theta-1.) * log(x) - 0.5 * omega * (x/eta + eta/x) ));

} /* end of _unur_pdf_gig() */

/*---------------------------------------------------------------------------*/

double
_unur_logpdf_gig(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (x <= 0.)
    /* out of support */
    return -INFINITY;

  return ( (theta-1.) * log(x) - 0.5 * omega * (x/eta + eta/x) );

} /* end of _unur_logpdf_gig() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_gig(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (x <= 0.)
    /* out of support */
    return 0.;

  return ( exp( (theta-3.) * log(x) - 0.5 * omega * (x/eta + eta/x) )
	   * (eta*eta*omega + 2.*eta*(theta-1.)*x - omega*x*x) / (2*eta) );

} /* end of _unur_dpdf_gig() */


/*---------------------------------------------------------------------------*/

double
_unur_dlogpdf_gig(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (x <= 0.)
    /* out of support */
    return 0.;

  return ( -0.5*(omega*(1/eta - eta/(x*x))) + (theta-1.)/x );

} /* end of _unur_dlogpdf_gig() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_gig( UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;

  DISTR.mode =
    (eta*(-1. + sqrt(omega*omega + (theta-1.)*(theta-1.)) + theta))/omega;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_gig() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_gig( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter omega */
  if (omega <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"omega <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* check parameter eta */
  if (n_params == 3 && eta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"eta <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }


  /* copy parameters for standard form */
  DISTR.theta = theta;
  DISTR.omega = omega;

  /* default parameters */
  DISTR.eta  = 1.;

  /* copy optional parameters */
  if (n_params > 2)
    DISTR.eta = eta;
  n_params = 3;

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0.;          /* left boundary  */
    DISTR.domain[1] = INFINITY;    /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_gig() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_gig( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_GIG;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = _unur_stdgen_gig_init;
   
  /* functions */
  DISTR.pdf     = _unur_pdf_gig;     /* pointer to PDF                  */
  DISTR.logpdf  = _unur_logpdf_gig;  /* pointer to logPDF               */
  DISTR.dpdf    = _unur_dpdf_gig;    /* pointer to derivative of PDF    */
  DISTR.dlogpdf = _unur_dlogpdf_gig; /* pointer to derivative of logPDF */
  DISTR.cdf  = NULL;                 /* _unur_cdf_gig; pointer to CDF   */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   );
		 /* UNUR_DISTR_SET_PDFAREA ); */
                
  /* set parameters for distribution */
  if (_unur_set_params_gig(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant */
  /*    DISTR.LOGNORMCONSTANT = ? */

  /* mode and area below p.d.f. */
  _unur_upd_mode_gig(distr);
  /*    DISTR.area = ? */

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_gig;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_gig; /* funct for computing mode */
  /* DISTR.upd_area  = _unur_upd_area_gig; funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_gig() */

/*---------------------------------------------------------------------------*/
#undef theta
#undef omega
#undef eta
#undef DISTR
/*---------------------------------------------------------------------------*/
