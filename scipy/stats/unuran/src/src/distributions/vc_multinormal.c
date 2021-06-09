/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      vc_multinormal.c                                              *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *   [5] S. Kotz, N. Balakrishnan, and N.L. Johnson                          *
 *       Continuous Multivariate Distributions,                              *
 *       Volume 1: Models and Applications                                   *
 *       John Wiley & Sons, Inc., New York, 2000                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  distr: Multinormal distribution [5; ch.45, p.107]                        *
 *                                                                           *
 *  pdf:       f(x) = exp( -1/2 * (x-mu)^t . Sigma^(-1) . (x-mu) )           * 
 *  domain:    Reals^(dim)                                                   *
 *  constant:  1 / ( (2 pi)^(dim/2) * sqrt(det(Sigma)) )                     *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  mean  ... mu      (0-vector)                                      *
 *     1:  covar ... Sigma   (identity matrix)                               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:       f(x) = exp( -1/2 * x^t . x )                                  *
 *  domain:    Reals^(dim)                                                   *
 *  constant:  1/(2 pi)^(dim/2)                                              *
 *                                                                           *
 *  parameters:                                                              *
 *     none                                                                  *
 *                                                                           *
 *     mean  = (0,...,0)  ... 0-vector                                       *
 *     covar = identity matrix                                               *
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
#include <distr/distr.h>
#include <distr/cvec.h>
#include <specfunct/unur_specfunct_source.h>
#include <utils/matrix_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"

#ifdef USE_DEPRECATED_CODE
#  include <distr/deprecated_distr.h>
#endif

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "multinormal";

/*---------------------------------------------------------------------------*/
/* parameters */

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cvec
#define LOGNORMCONSTANT (distr->data.cvec.norm_constant)

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */

static double _unur_pdf_multinormal( const double *x, UNUR_DISTR *distr );
static double _unur_logpdf_multinormal( const double *x, UNUR_DISTR *distr );
static int _unur_dlogpdf_multinormal( double *result, const double *x, UNUR_DISTR *distr );
static double _unur_pdlogpdf_multinormal( const double *x, int coord, UNUR_DISTR *distr );

static int _unur_upd_mode_multinormal( UNUR_DISTR *distr );
static int _unur_upd_volume_multinormal( UNUR_DISTR *distr );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_multinormal( const double *x, UNUR_DISTR *distr )
{ 
  return exp( _unur_logpdf_multinormal( x, distr ) );
} /* end of _unur_pdf_multinormal() */

/*---------------------------------------------------------------------------*/

double
_unur_logpdf_multinormal( const double *x, UNUR_DISTR *distr )
{ 
#define idx(a,b) ((a)*dim+(b))

  int i,j, dim;
  double *mean;
  const double *covar_inv; 
  
  double xx; /* argument used in the evaluation of exp(-xx/2) */
  double cx; /* element of multiplication of covariance matrix and x */
  
  dim = distr->dim;
  
  if (DISTR.mean == NULL) {
    if (DISTR.covar != NULL) {
      _unur_warning(distr->name,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    }
    /* standard form */
    xx=0.;
    for (i=0; i<dim; i++) { xx += x[i]*x[i]; }
    return (-xx/2. + LOGNORMCONSTANT);  
  }

  /* mean vector */
  mean = DISTR.mean;

  /* get inverse of covariance matrix */
  covar_inv = unur_distr_cvec_get_covar_inv(distr);
  if (covar_inv==NULL) 
    /* inverse of covariance matrix not available */
    return INFINITY;

  xx=0.; /* resetting exponential function argument */
  for (i=0; i<dim; i++) {
    cx=0.; 
    /* multiplication of inverse covariance matrix and (x-mean) */
    for (j=0; j<dim; j++) {
      cx += covar_inv[idx(i,j)] * (x[j]-mean[j]);
    }
    xx += (x[i]-mean[i])*cx;
  }
  
  return (-xx/2. + LOGNORMCONSTANT);

#undef idx
} /* end of _unur_logpdf_multinormal() */

/*---------------------------------------------------------------------------*/

int
_unur_dlogpdf_multinormal( double *result, const double *x, UNUR_DISTR *distr )
{
#define idx(a,b) ((a)*dim+(b))

  int i,j, dim;
  double *mean;
  const double *covar_inv;
    
  dim = distr->dim;
  mean = DISTR.mean;

  /* get inverse of covariance matrix */
  covar_inv = unur_distr_cvec_get_covar_inv(distr);
  if (covar_inv==NULL) 
    /* inverse of covariance matrix not available */
    return UNUR_FAILURE;

  for (i=0; i<dim; i++) {
    result[i] = 0.;
    for (j=0; j<dim; j++) 
      result[i] += -0.5 * (x[j]-mean[j]) * (covar_inv[idx(i,j)]+covar_inv[idx(j,i)]);
  }
  
  return UNUR_SUCCESS; 

#undef idx
} /* end of _unur_dlogpdf_multinormal() */

/*---------------------------------------------------------------------------*/

double
_unur_pdlogpdf_multinormal( const double *x, int coord, UNUR_DISTR *distr )
{
#define idx(a,b) ((a)*dim+(b))

  int j, dim;
  double *mean;
  const double *covar_inv;
  double result;
    
  dim = distr->dim;
  mean = DISTR.mean;

  /* check arguments */
  if (coord < 0 || coord >= dim) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DOMAIN,"invalid coordinate");
    return INFINITY;
  }

  /* get inverse of covariance matrix */
  covar_inv = unur_distr_cvec_get_covar_inv(distr);
  if (covar_inv==NULL) 
    /* inverse of covariance matrix not available */
    return INFINITY;

  result = 0.;
  for (j=0; j<dim; j++)
    result += -0.5 * (x[j]-mean[j]) * (covar_inv[idx(coord,j)]+covar_inv[idx(j,coord)]);
  
  return result;

#undef idx
} /* end of _unur_pdlogpdf_multinormal() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_multinormal( UNUR_DISTR *distr )
{
  /* TODO: checking if mode is inside domain */

  if (DISTR.mode == NULL) 
    DISTR.mode = _unur_xmalloc( distr->dim * sizeof(double) );
  memcpy( DISTR.mode, DISTR.mean, distr->dim * sizeof(double) );

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_multinormal() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_volume_multinormal( UNUR_DISTR *distr )
{
  /* TODO: checking for modified domain */

  double det_covar;

  /* log of normalization constant */
  det_covar = (DISTR.covar == NULL) 
    ? 1. : _unur_matrix_determinant(distr->dim, DISTR.covar);
  LOGNORMCONSTANT = - ( distr->dim * log(2 * M_PI) + log(det_covar) ) / 2.;

  return UNUR_SUCCESS;
} /* end of _unur_upd_volume_multinormal() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Make distribution object                                               **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_multinormal( int dim, const double *mean, const double *covar )
   
/*   
   @code{UNUR_DISTR *unur_distr_multinormal(int dim, const double *mean, const double *covar)}
   creates a distribution object for the multinormal distribution with
   @var{dim} components. @var{mean} is an array of size @var{dim}.
   A NULL pointer for @var{mean} is interpreted as the zero
   vector (0,@dots{},0).
   @var{covar} is an array of size @var{dim}x@var{dim} and holds the
   covariance matrix, where the rows of the matrix are stored
   consecutively in this array. The NULL pointer can be used
   instead the identity matrix.
   If @var{covar} is not a valid covariance matrix (i.e., not positive
   definite) then no distribution object is created and NULL is returned.

   For standard form of the distribution use the null vector for @var{mean} and 
   the identity matrix for @var{covar}.
*/

{
  struct unur_distr *distr;
  double det_covar; /* determinant of covariance matrix */

  /* get new (empty) distribution object */
  distr = unur_distr_cvec_new(dim);

  /* check new parameter for generator */
  if (distr == NULL) {
    /* error: dim < 1 */
    return NULL;
  }

  /* set distribution id */
  distr->id = UNUR_DISTR_MNORMAL;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = _unur_stdgen_multinormal_init;

  /* copy (and check) parameters */
  if ((unur_distr_cvec_set_mean(distr,mean)!=UNUR_SUCCESS) ||
      (unur_distr_cvec_set_covar(distr,covar)!=UNUR_SUCCESS) ) {
    unur_distr_free( distr );
    return NULL;
  }

  /* functions */
  DISTR.pdf      = _unur_pdf_multinormal;       /* pointer to PDF */
  DISTR.logpdf   = _unur_logpdf_multinormal;    /* pointer to logPDF */
  DISTR.dpdf     = _unur_distr_cvec_eval_dpdf_from_dlogpdf;  /* pointer to gradient of PDF */
  DISTR.dlogpdf  = _unur_dlogpdf_multinormal;    /* pointer to gradient of logPDF */
  DISTR.pdpdf    = _unur_distr_cvec_eval_pdpdf_from_pdlogpdf;  /* pointer to part. deriv. of PDF */
  DISTR.pdlogpdf = _unur_pdlogpdf_multinormal;  /* pointer to partial derivative of logPDF */

#ifdef USE_DEPRECATED_CODE
  {
    /* set standardized marginal distributions */
    struct unur_distr *stdmarginal = unur_distr_normal(NULL,0);
    unur_distr_cvec_set_stdmarginals(distr,stdmarginal);
    unur_distr_free(stdmarginal);
  }
#endif

  /* copy other parameters of distribution */
  /* none */

  /* number of other parameters */
  /* DISTR.n_params = 0;  ... default */

  /* domain */

  /* log of normalization constant */
  det_covar = (DISTR.covar == NULL) ? 1. : _unur_matrix_determinant(dim, DISTR.covar);
  LOGNORMCONSTANT = - ( distr->dim * log(2 * M_PI) + log(det_covar) ) / 2.;

  /* mode */
  DISTR.mode = _unur_xmalloc( distr->dim * sizeof(double) );
  memcpy( DISTR.mode, DISTR.mean, distr->dim * sizeof(double) );

  /* volume below p.d.f. */
  DISTR.volume = 1.; 

  /* function for updating derived parameters */
  DISTR.upd_mode   = _unur_upd_mode_multinormal;   /* funct for computing mode   */
  DISTR.upd_volume = _unur_upd_volume_multinormal; /* funct for computing volume */
                
  /* indicate which parameters are set (additional to mean and covariance) */
  distr->set |= ( UNUR_DISTR_SET_STDDOMAIN |
		  UNUR_DISTR_SET_PDFVOLUME |
		  UNUR_DISTR_SET_MODE );

  /* return pointer to object */
  return distr;

} /* end of unur_distr_multinormal() */

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/
