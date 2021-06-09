/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_gig2.c                                                     *
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
 *  Alternative parametrization.                                             *
 *                                                                           *
 *  pdf:       f(x) = x^(theta-1) * exp( -1/2 * (chi/x + psi*x))             *
 *  domain:    0 < x < infinity                                              *
 *  constant:  (psi/chi)^(theta/2) / (2*K_theta(sqrt(psi*chi)))              *
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
 *     0:  theta         ... shape                                           *
 *     1:  psi   > 0     ... shape                                           *
 *     2:  chi   > 0     ... shape                                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Reparametrization towards "gig":                                         *
 *     omega = sqrt(psi*chi)                                                 *
 *     eta   = sqrt(chi/psi)                                                 *
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

#ifdef USE_EXPERIMENTAL_CODE
#include <gsl/gsl_integration.h>
#endif

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "gig2";

/* parameters */
#define theta  params[0]    /* shape */
#define psi    params[1]    /* shape */
#define chi    params[2]    /* shape */

#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)

/* function prototypes                                                       */
static double _unur_pdf_gig2( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_gig2( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_gig2( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_gig2( double x, const UNUR_DISTR *distr );
#ifdef USE_EXPERIMENTAL_CODE
static double _unur_cdf_gig2( double x, const UNUR_DISTR *distr );
#endif

static int _unur_upd_mode_gig2( UNUR_DISTR *distr );
static double _unur_normconstant_gig2( const double *params, int n_params );
static int _unur_set_params_gig2( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_gig2(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (x <= 0.)
    /* out of support */
    return 0.;

  return (NORMCONSTANT * exp( (theta-1.) * log(x) - 0.5 * (chi/x + psi*x) ));

} /* end of _unur_pdf_gig2() */

/*---------------------------------------------------------------------------*/

double
_unur_logpdf_gig2(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (x <= 0.)
    /* out of support */
    return -INFINITY;

  return ( (theta-1.) * log(x) - 0.5 * (chi/x + psi*x) + log(NORMCONSTANT) );

} /* end of _unur_logpdf_gig2() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_gig2(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (x <= 0.)
    /* out of support */
    return 0.;

  return ( NORMCONSTANT * 0.5 * exp( (theta-3.) * log(x) - (chi + psi*x*x)/(2*x) )
	   * (chi - x*(2 - 2*theta + psi*x)) );

} /* end of _unur_dpdf_gig2() */


/*---------------------------------------------------------------------------*/

double
_unur_dlogpdf_gig2(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (x <= 0.)
    /* out of support */
    return 0.;

  return ( -0.5*(psi - chi/(x*x)) + (theta-1.)/x  + log(NORMCONSTANT) ) ;

} /* end of _unur_dlogpdf_gig2() */

/*---------------------------------------------------------------------------*/

#ifdef USE_EXPERIMENTAL_CODE

#ifndef _unur_SF_bessel_k
#error run ./configure with flag --with-Rmath
#endif 

double
_unur_cdf_gig2(double x, const UNUR_DISTR *distr)
{
  double epsabs = 1.e-13;
  double epsrel = 1.e-13;
  double result, abserr;
  size_t limit = 1000;
  gsl_integration_workspace *work;

  if (x <= 0.)
    /* out of support */
    return 0.;

  /* use GSL function */
  work = gsl_integration_workspace_alloc (limit);
  gsl_function F;
  F.function = _unur_pdf_gig2;
  F.params = distr;

  gsl_integration_qag (&F, 0, x, epsabs, epsrel, limit, /* key = */ GSL_INTEG_GAUSS61,
		       work, &result, &abserr);

  /* The integration rule is determined by the value of 'key', 
   * which should be chosen from the following symbolic names,
   * 
   *   GSL_INTEG_GAUSS15  (key = 1)
   *   GSL_INTEG_GAUSS21  (key = 2)
   *   GSL_INTEG_GAUSS31  (key = 3)
   *   GSL_INTEG_GAUSS41  (key = 4)
   *   GSL_INTEG_GAUSS51  (key = 5)
   *   GSL_INTEG_GAUSS61  (key = 6)
   */

  gsl_integration_workspace_free (work);

  return result;

} /* end of _unur_cdf_gig2() */

#endif

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_gig2( UNUR_DISTR *distr )
{
  register const double *params = DISTR.params;

  DISTR.mode =
    ((theta-1.)+sqrt((theta-1.)*(theta-1.) + psi*chi)) / psi;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_gig2() */

/*---------------------------------------------------------------------------*/

double
_unur_normconstant_gig2(const double *params ATTRIBUTE__UNUSED, int n_params ATTRIBUTE__UNUSED)
{ 
#ifdef _unur_SF_bessel_k

  /*
   *  constant:  (psi/chi)^(theta/2) / (2*K_theta(sqrt(psi*chi)))
   */

  return ( pow(psi/chi, theta/2.) / (2. * _unur_SF_bessel_k(sqrt(psi*chi),theta)) );
#else
  return 1.;
#endif
} /* end of _unur_normconstant_gig2() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_gig2( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 3) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); 
    return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter omega */
  if (psi <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"psi <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* check parameter eta */
  if (chi <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"chi <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }


  /* copy parameters for standard form */
  DISTR.theta = theta;
  DISTR.psi = psi;
  DISTR.chi = chi;

  /* default parameters: none */

  /* copy optional parameters: none */

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0.;          /* left boundary  */
    DISTR.domain[1] = INFINITY;    /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_gig2() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_gig2( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_GIG2;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  /* DISTR.init = _unur_stdgen_gig2_init; */
   
  /* functions */
  DISTR.pdf     = _unur_pdf_gig2;     /* pointer to PDF                  */
  DISTR.logpdf  = _unur_logpdf_gig2;  /* pointer to logPDF               */
  DISTR.dpdf    = _unur_dpdf_gig2;    /* pointer to derivative of PDF    */
  DISTR.dlogpdf = _unur_dlogpdf_gig2; /* pointer to derivative of logPDF */
#ifdef USE_EXPERIMENTAL_CODE
  DISTR.cdf     = _unur_cdf_gig2;     /* pointer to CDF                  */
#endif

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   );
		 /* UNUR_DISTR_SET_PDFAREA ); */
                
  /* set parameters for distribution */
  if (_unur_set_params_gig2(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* normalization constant */
  NORMCONSTANT = _unur_normconstant_gig2(params,n_params);

  /* mode and area below p.d.f. */
  _unur_upd_mode_gig2(distr);
  /*    DISTR.area = ? */

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_gig2;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_gig2; /* funct for computing mode */
  /* DISTR.upd_area  = _unur_upd_area_gig2; funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_gig2() */

/*---------------------------------------------------------------------------*/
#undef theta
#undef psi
#undef chi
#undef DISTR
/*---------------------------------------------------------------------------*/
