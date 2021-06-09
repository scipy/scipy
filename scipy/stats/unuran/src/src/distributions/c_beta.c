/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_beta.c                                                     *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *   [3] N.L. Johnson, S. Kotz and N. Balakrishnan                           *
 *       Continuous Univariate Distributions,                                *
 *       Volume 2, 2nd edition                                               *
 *       John Wiley & Sons, Inc., New York, 1995                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  distr: Beta distribution [3; ch.25, p.210]                               *
 *                                                                           *
 *  pdf:       f(x) = (x-a)^(p-1) * (b-x)^(q-1)                              *
 *  domain:    a < x < b                                                     *
 *  constant:  1 / ( Beta(p,q) * (b-a)^(p+q-1) )                             *
 *                                                                           *
 *  parameters: 4                                                            *
 *     0:  p > 0        ... shape                                            *
 *     1:  q > 0        ... shape                                            *
 *     2:  a      (0)   ... location                                         *
 *     3:  b >a   (1)   ... location                                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:       f(x) = x^(p-1) * (1-x)^(q-1)                                  *
 *  domain:    0 < x < 1                                                     *
 *  constant:  1 / Beta(p,q)                                                 *
 *                                                                           *
 *  parameters: 2                                                            *
 *     0:  p > 0    ... shape                                                *
 *     1:  q > 0    ... shape                                                *
 *                                                                           *
 *     2:  a = 0                                                             *
 *     3:  b = 1                                                             *
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

static const char distr_name[] = "beta";

/*---------------------------------------------------------------------------*/
/* parameters */
#define p  params[0]
#define q  params[1]
#define a  params[2]
#define b  params[3]

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)

/*---------------------------------------------------------------------------*/

/* function prototypes                                                       */
static double _unur_pdf_beta( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_beta( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_beta( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_beta( double x, const UNUR_DISTR *distr );
static double _unur_cdf_beta( double x, const UNUR_DISTR *distr );
#ifdef _unur_SF_invcdf_beta
static double _unur_invcdf_beta( double x, const UNUR_DISTR *distr );
#endif

static int _unur_upd_mode_beta( UNUR_DISTR *distr );
static int _unur_upd_area_beta( UNUR_DISTR *distr );
inline static double _unur_lognormconstant_beta( const double *params, int n_params );
static int _unur_set_params_beta( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_beta(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (DISTR.n_params > 2)
    /* standardize */
    x = (x-a) / (b-a);

  /* standard form */

  if (x > 0. && x < 1.)
    return exp((p-1.)*log(x) + (q-1.)*log(1.-x) - LOGNORMCONSTANT);

  if ((_unur_iszero(x) && _unur_isone(p)) || (_unur_isone(x) && _unur_isone(q)))
    return exp(-LOGNORMCONSTANT);

  if ((_unur_iszero(x) && p<1.) || (_unur_isone(x) && q<1.))
    return INFINITY;

  /* out of support */
  return 0.;

} /* end of _unur_pdf_beta() */

/*---------------------------------------------------------------------------*/

double
_unur_logpdf_beta(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (DISTR.n_params > 2)
    /* standardize */
    x = (x-a) / (b-a);

  /* standard form */

  if (x > 0. && x < 1.)
    return ((p-1.)*log(x) + (q-1.)*log(1.-x) - LOGNORMCONSTANT);

  if ((_unur_iszero(x) && _unur_isone(p)) 
      || (_unur_isone(x) && _unur_isone(q)))
    return (-LOGNORMCONSTANT);

  if ((_unur_iszero(x) && p<1.) || (_unur_isone(x) && q<1.))
    return INFINITY;

  /* out of support */
  return -INFINITY;

} /* end of _unur_logpdf_beta() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_beta(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (DISTR.n_params > 2) {
    /* standardize */
    x = (x-a) / (b-a);
  }

  /* standard form */

  if (x > 0. && x < 1.)
    return (exp((p-2.)*log(x) + (q-2.)*log(1.-x) - LOGNORMCONSTANT) * ( (p-1.)*(1.-x) - (q-1.)*x ) / (b-a) );

  if (_unur_iszero(x) && _unur_isone(p))
    return (1.-q)*exp(-LOGNORMCONSTANT)/(b-a);

  if (_unur_iszero(x) && _unur_isfsame(p,2.))
    return exp(-LOGNORMCONSTANT)/(b-a);

  if (_unur_iszero(x) && p<2.)
    return (p>1. ? INFINITY : -INFINITY);

  /*   if (_unur_iszero(x) && p>2.) */
  /*     return 0.; */

  if (_unur_isone(x) && _unur_isone(q))
    return (p-1.)*exp(-LOGNORMCONSTANT)/(b-a);

  if (_unur_isone(x) && _unur_isfsame(q,2.))
    return -exp(-LOGNORMCONSTANT)/(b-a);

  if (_unur_isone(x) && q<2.)
    return (q>1. ? -INFINITY : INFINITY);

  /*   if (_unur_isone(x) && q>2.) */
  /*     return 0.; */

  /* out of support */
  return 0.;

} /* end of _unur_dpdf_beta() */

/*---------------------------------------------------------------------------*/

double
_unur_dlogpdf_beta(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (DISTR.n_params > 2) {
    /* standardize */
    x = (x-a) / (b-a);
  }

  /* standard form */

  if (x > 0. && x < 1.)
    return (((p-1.)/x - (q-1.)/(1.-x)) / (b-a));

  if (_unur_iszero(x) && p<1.)
    return -INFINITY;

  if (_unur_iszero(x) && _unur_isone(p))
    return (-(q-1.)/((1.-x)*(b-a)));

  if (_unur_iszero(x) && p>1.)
    return INFINITY;

  if (_unur_isone(x) && q<1.)
    return INFINITY;

  if (_unur_isone(x) && _unur_isone(q))
    return ((p-1.)/(b-a));

  if (_unur_isone(x) && q>1.)
    return -INFINITY;

  /* out of support */
  return 0.;

} /* end of _unur_dlogpdf_beta() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_beta(double x, const UNUR_DISTR *distr)
{
  register const double *params = DISTR.params;

  if (DISTR.n_params > 2)
    /* standardize */
    x = (x-a) / (b-a);

  /* standard form */

  /* out of support of p.d.f.? */
  if (x <= 0.) return 0.;
  if (x >= 1.) return 1.;

  return _unur_SF_incomplete_beta(x,p,q);

} /* end of _unur_cdf_beta() */

/*---------------------------------------------------------------------------*/

#ifdef _unur_SF_invcdf_beta

double
_unur_invcdf_beta(double x, const UNUR_DISTR *distr)
{
  const double *params = DISTR.params;
  
  if (DISTR.n_params == 2)
    return _unur_SF_invcdf_beta(x,p,q);
  else
    return (a + _unur_SF_invcdf_beta(x,p,q))*(b-a);

} /* end of _unur_invcdf_beta() */
#endif

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_beta( UNUR_DISTR *distr )
{
  register double *params = DISTR.params;

  if (p <= 1. && q > 1.)
    DISTR.mode = 0.;              /* left limit of domain */

  else if (p > 1. && q <= 1.)
    DISTR.mode = 1.;              /* right limit of domain */

  else if (p > 1. && q > 1.)
    DISTR.mode = (p - 1.) / (p + q - 2.);
  
  else {
    /* p.d.f. is not unimodal */
    DISTR.mode = INFINITY;
    return UNUR_ERR_DISTR_PROP;
  }

  if (DISTR.n_params > 2)
    DISTR.mode = DISTR.mode * (b - a) + a;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_beta() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_beta( UNUR_DISTR *distr )
{
  /* log of normalization constant */
  LOGNORMCONSTANT = _unur_lognormconstant_beta(DISTR.params,DISTR.n_params);
  
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  
  /* else */
  DISTR.area = ( _unur_cdf_beta( DISTR.domain[1],distr) 
		 - _unur_cdf_beta( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;

} /* end of _unur_upd_area_beta() */

/*---------------------------------------------------------------------------*/

double
_unur_lognormconstant_beta(const double *params, int n_params)
{ 
  if (n_params > 2)
    /* non-standard form */
    /* log( Beta(p,q) * (b-a) ) */
    return (_unur_SF_ln_gamma(p) + _unur_SF_ln_gamma(q) - _unur_SF_ln_gamma(p+q) + log(b-a) );

  else
    /* standard form */
    /* log( Beta(p,q) ) */
    return (_unur_SF_ln_gamma(p) + _unur_SF_ln_gamma(q) - _unur_SF_ln_gamma(p+q));

} /* end of _unur_lognormconstant_beta() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_beta( UNUR_DISTR *distr, const double *params, int n_params )
{

  /* check number of parameters for distribution */
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params == 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"");
    n_params = 2; }
  if (n_params > 4) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 4; }
  CHECK_NULL(params,UNUR_ERR_NULL);


  /* check parameters p and q */
  if (p <= 0. || q <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"p <= 0 or q <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* check parameters a and b */
  if (n_params > 2 && a >= b) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"a >= b");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form */
  DISTR.p = p;
  DISTR.q = q;

  /* copy optional parameters */
  if (n_params > 2) {
    DISTR.a = a;
    DISTR.b = b;
  }
  else { /* or use defaults */
    DISTR.a = 0.;      /* default for a */
    DISTR.b = 1.;      /* default for b */
  }

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = DISTR.a; /* left boundary  */
    DISTR.domain[1] = DISTR.b; /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_beta() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_beta( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_BETA;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = _unur_stdgen_beta_init;

  /* functions */
  DISTR.pdf     = _unur_pdf_beta;     /* pointer to PDF                  */
  DISTR.logpdf  = _unur_logpdf_beta;  /* pointer to logPDF               */
  DISTR.dpdf    = _unur_dpdf_beta;    /* pointer to derivative of PDF    */
  DISTR.dlogpdf = _unur_dlogpdf_beta; /* pointer to derivative of logPDF */
  DISTR.cdf     = _unur_cdf_beta;     /* pointer to CDF                  */
#ifdef _unur_SF_invcdf_beta
  DISTR.invcdf  = _unur_invcdf_beta;  /* pointer to inverse CDF          */
#endif

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA |
		 UNUR_DISTR_SET_MODE );

  /* set parameters for distribution */
  if (_unur_set_params_beta(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant */
  LOGNORMCONSTANT = _unur_lognormconstant_beta(DISTR.params,DISTR.n_params);

  /* mode and area below p.d.f. */
  _unur_upd_mode_beta( distr );
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_beta;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_beta; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_beta; /* funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_beta() */

/*---------------------------------------------------------------------------*/
#undef p
#undef q
#undef a
#undef b
#undef DISTR
/*---------------------------------------------------------------------------*/
