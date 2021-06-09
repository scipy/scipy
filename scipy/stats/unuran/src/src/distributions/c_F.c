/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_F.c                                                        *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *   [3] N.L. Johnson, S. Kotz, and N. Balakrishnan                          *
 *       Continuous Univariate Distributions,                                *
 *       Volume 2, 2nd edition                                               *
 *       John Wiley & Sons, Inc., New York, 1995                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  distr: F distribution   [3; ch.27, p.332]                                *
 *                                                                           *
 *  pdf:       f(x) = (x^{{nu_1}/2-1}) / ((1+nu_1/nu_2 x)^{(nu_1+nu_2)/2})   *
 *  domain:    0 < x < infinity                                              *
 *  constant:  \frac{(nu_1/nu_2)^{{nu_1}/2}/Beta({nu_1}/2,{nu_2}/2)          *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  nu_1 > 0   ... shape (degrees of freedom)                         *
 *     1:  nu_2 > 0   ... shape (degrees of freedom)                         *
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

static const char distr_name[] = "F";

/*---------------------------------------------------------------------------*/
/* parameters */
#define nua  params[0]
#define nub  params[1]

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */
static double _unur_pdf_F( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_F( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_F( double x, const UNUR_DISTR *distr );
static double _unur_dlogpdf_F( double x, const UNUR_DISTR *distr );
static double _unur_cdf_F( double x, const UNUR_DISTR *distr );
#ifdef _unur_SF_invcdf_F
static double _unur_invcdf_F( double x, const UNUR_DISTR *distr );
#endif

static int _unur_upd_mode_F( UNUR_DISTR *distr );
static int _unur_upd_area_F( UNUR_DISTR *distr );
inline static double _unur_lognormconstant_F( const double *params, int n_params );
static int _unur_set_params_F( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_F(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (x < 0.)
    /* out of support */
    return 0.;

  else if (_unur_iszero(x)) {
    if (nua < 2.)
      return INFINITY;
    else if (_unur_isfsame(nua,2.))
      return exp(-LOGNORMCONSTANT);
    else
      return 0.;
  }

  else 
    return exp ((nua/2. - 1.)*log(x) - 0.5*(nua + nub)*log(1. + x * nua / nub) - LOGNORMCONSTANT);

} /* end of _unur_pdf_F() */

/*---------------------------------------------------------------------------*/

double
_unur_logpdf_F(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (x < 0.)
    /* out of support */
    return -INFINITY;

  else if (_unur_iszero(x)) {
    if (nua < 2.)
      return INFINITY;
    else if (_unur_isfsame(nub,2.))
      return -LOGNORMCONSTANT;
    else
      return -INFINITY;
  }

  else 
    return ((nua/2. - 1.)*log(x) - 0.5*(nua + nub)*log(1. + x * nua / nub) - LOGNORMCONSTANT);

} /* end of _unur_logpdf_F() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_F(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (x < 0.)
    /* out of support */
    return 0.;

  else if (_unur_iszero(x)) {
    if (nua < 2.)
      return -INFINITY;
    else if (_unur_isfsame(nub,2.))
      return -(2.+nub)/nub * exp(-LOGNORMCONSTANT);
    else
      return 0.;
  }

  else
    return _unur_pdf_F(x,distr) * _unur_dlogpdf_F(x,distr);

} /* end of _unur_dpdf_F() */

/*---------------------------------------------------------------------------*/

double
_unur_dlogpdf_F(double x, const UNUR_DISTR *distr)
{ 
  register const double *params = DISTR.params;

  if (x < 0.)
    /* out of support */
    return 0.;

  else if (_unur_iszero(x)) {
    if (nua < 2.)
      return -INFINITY;
    else if (_unur_isfsame(nub,2.))
      return -(2.+nub)/nub;
    else
      return INFINITY;
  }

  else
    return ((nua/2.-1.)/x - nua*(nua+nub)/(2.*nub)/(1.+x*nua/nub));

} /* end of _unur_dlogpdf_F() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_F(double x, const UNUR_DISTR *distr)
{ 
#ifdef _unur_SF_cdf_F
  return _unur_SF_cdf_F(x,DISTR.nua,DISTR.nub);
#else
  const double *params = DISTR.params;

  if (x <= 0.)
    /* out of support of p.d.f. */
    return 0.;

  if (nua * x > nub)
    return 1. - _unur_SF_incomplete_beta(nub / (nub + nua * x), nub/2., nua/2.);
  else
    return _unur_SF_incomplete_beta(nua * x / (nub + nua * x), nua/2., nub/2.);
#endif
} /* end of _unur_cdf_chisquare() */

/*---------------------------------------------------------------------------*/

#ifdef _unur_SF_invcdf_F
double
_unur_invcdf_F(double x, const UNUR_DISTR *distr)
{
  return _unur_SF_invcdf_F(x,DISTR.nua,DISTR.nub);
} /* end of _unur_invcdf_F() */
#endif

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_F( UNUR_DISTR *distr )
{
  if (DISTR.nua >= 2.)
    DISTR.mode = ((DISTR.nua - 2.) * DISTR.nub) / (DISTR.nua * (DISTR.nub + 2.));
  else
    DISTR.mode = 0.;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_F() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_F( UNUR_DISTR *distr )
{
  /* log of normalization constant */
  LOGNORMCONSTANT = _unur_lognormconstant_F(DISTR.params,DISTR.n_params);
  
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }

  /* else */
  DISTR.area = ( _unur_cdf_F( DISTR.domain[1],distr)
		 - _unur_cdf_F( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
  
} /* end of _unur_upd_area_F() */

/*---------------------------------------------------------------------------*/

double
_unur_lognormconstant_F(const double *params, int n_params ATTRIBUTE__UNUSED)
{ 
  /* log( Beta(nu1/2, nu2/2) ) - (nu1/2) * log(nu1 / nu2) */
  return ((_unur_SF_ln_gamma(nua/2.) + _unur_SF_ln_gamma(nub/2.) - _unur_SF_ln_gamma((nua+nub)/2.))
	  - 0.5 * nua * log(nua/nub));
} /* end of _unur_lognormconstant_F() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_F( UNUR_DISTR *distr, const double *params, int n_params )
{

  /* check number of parameters for distribution */
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameters nu1 and nu2 */
  if (nua <= 0. || nub <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"nu <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form */
  DISTR.nua = nua;
  DISTR.nub = nub;

  /* copy optional parameters: none */

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0.;          /* left boundary  */
    DISTR.domain[1] = INFINITY;    /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_F() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_F( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_F;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = NULL;
   
  /* functions */
  DISTR.pdf     = _unur_pdf_F;           /* pointer to PDF                  */
  DISTR.logpdf  = _unur_logpdf_F;        /* pointer to logPDF               */
  DISTR.dpdf    = _unur_dpdf_F;          /* pointer to derivative of PDF    */
  DISTR.dlogpdf = _unur_dlogpdf_F;       /* pointer to derivative of logPDF */
  DISTR.cdf     = _unur_cdf_F;           /* pointer to CDF                  */
#ifdef _unur_SF_invcdf_student
  DISTR.invcdf  = _unur_invcdf_F;        /* pointer to inverse CDF          */
#endif

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA |
		 UNUR_DISTR_SET_MODE );
                
  /* set parameters for distribution */
  if (_unur_set_params_F(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant */
  LOGNORMCONSTANT = _unur_lognormconstant_F(DISTR.params,DISTR.n_params);

  /* mode and area below p.d.f. */
  _unur_upd_mode_F( distr );
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_F;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_F;   /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_F;   /* funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_F() */

/*---------------------------------------------------------------------------*/
#undef nu1
#undef nu2
#undef DISTR
/*---------------------------------------------------------------------------*/
