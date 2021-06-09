/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      vc_copula.c                                                  *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  distr: Copula                                                            *
 *                                                                           *
 *  domain:    (0,1)^(dim)                                                   *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  rankcorr ... rank correlation matrix  (identity matrix)           *
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
/* #include <specfunct/unur_specfunct_source.h> */
/* #include <utils/matrix_source.h> */
#include "unur_distributions.h"
#include "unur_distributions_source.h"
/* #include "unur_stddistr.h" */

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "copula";

/*---------------------------------------------------------------------------*/
/* parameters */

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cvec
/* #define LOGNORMCONSTANT (distr->data.cvec.norm_constant) */

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */

/*     static double _unur_pdf_copula( const double *x, UNUR_DISTR *distr ); */
/*     static double _unur_logpdf_copula( const double *x, UNUR_DISTR *distr ); */
/*     static int _unur_dlogpdf_copula( double *result, const double *x, UNUR_DISTR *distr ); */
/** TODO:
    static int _unur_upd_mode_copula( UNUR_DISTR *distr );
    static int _unur_upd_area_copula( UNUR_DISTR *distr );
**/

/*---------------------------------------------------------------------------*/

/* double */
/* _unur_pdf_copula( const double *x, UNUR_DISTR *distr ) */
/* {  */
/*   return 1; */
/* } /\* end of _unur_pdf_copula() *\/ */

/*---------------------------------------------------------------------------*/

/* int */
/* _unur_dpdf_copula( double *result, const double *x, UNUR_DISTR *distr ) */
/* { */
/*   return UNUR_SUCCESS;  */
/* } /\* end of _unur_dlogpdf_copula() *\/ */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Make distribution object                                               **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_copula( int dim, const double *rankcorr )
/*
   @code{UNUR_DISTR *unur_distr_copula(int dim, const double *rankcorr)}
   creates a distribution object for a copula with @var{dim} components. 
   @var{rankcorr} is an array of size @var{dim}x@var{dim} and holds the
   rank correlation matrix (Spearman's correlation), where the rows of
   the matrix are stored consecutively in this array. The NULL pointer
   can be used instead the identity matrix.

   If @var{covar} is not a valid rank correlation matrix (i.e., not positive
   definite) then no distribution object is created and NULL is returned.
*/   
{
  struct unur_distr *distr;
  struct unur_distr *marginal;

  /* get new (empty) distribution object */
  distr = unur_distr_cvec_new(dim);

  /* check new parameter for generator */
  if (distr == NULL) {
    /* error: dim < 1 */
    return NULL;
  }

  /* set distribution id */
  distr->id = UNUR_DISTR_COPULA;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = NULL;

  /* copy (and check) parameters */
  if ( unur_distr_cvec_set_rankcorr(distr,rankcorr)!=UNUR_SUCCESS ) {
    unur_distr_free( distr );
    return NULL;
  }

  /* functions: none */

  /* set marginal distributions */
  marginal = unur_distr_uniform(NULL,0);
  unur_distr_cvec_set_marginals(distr,marginal);
  unur_distr_free(marginal);

  /* copy other parameters of distribution */
  /* none */

  /* number of other parameters */
  /* DISTR.n_params = 0;  ... default */

  /* domain */

  /* normalization constant */

  /* mode */

  /* volume below p.d.f. */

  /* indicate which parameters are set (additional to mean and covariance) */
  distr->set |= ( 0 );

  /* return pointer to object */
  return distr;

} /* end of unur_distr_copula() */

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/
