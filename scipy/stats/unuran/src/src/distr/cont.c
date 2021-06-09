/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      cont.c                                                       *
 *                                                                           *
 *   manipulate univariate continuous distribution objects                   *
 *                                                                           *
 *   return:                                                                 *
 *     UNUR_SUCCESS ... on success                                           *
 *     error code   ... on error                                             *
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
#include <distributions/unur_stddistr.h>
#include <parser/functparser_source.h>
#include <utils/fmax_source.h>
#include "distr_source.h"
#include "distr.h"
#include "cont.h"

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cont

/* for derived distributions (e.g. order statistics):
   data of underlying distributions */
#define BASE  distr->base->data.cont

/*---------------------------------------------------------------------------*/

static double _unur_distr_cont_eval_pdf_tree( double x, const struct unur_distr *distr );
static double _unur_distr_cont_eval_logpdf_tree( double x, const struct unur_distr *distr );
static double _unur_distr_cont_eval_dpdf_tree( double x, const struct unur_distr *distr );
static double _unur_distr_cont_eval_dlogpdf_tree( double x, const struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* evaluate function tree for (derivative of) (log) PDF.                     */
/*---------------------------------------------------------------------------*/

static double _unur_distr_cont_eval_cdf_tree( double x, const struct unur_distr *distr );
static double _unur_distr_cont_eval_logcdf_tree( double x, const struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* evaluate function tree for (log) CDF.                                     */
/*---------------------------------------------------------------------------*/

static double _unur_distr_cont_eval_hr_tree( double x, const struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* evaluate function tree for HR.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_distr_cont_free( struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* destroy distribution object.                                              */
/*---------------------------------------------------------------------------*/

static int _unur_distr_cont_find_mode( struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* find mode of unimodal univariate PDF numerically                          */
/*---------------------------------------------------------------------------*/

static double _unur_aux_pdf(double x, void *p);
/*---------------------------------------------------------------------------*/
/* Auxiliary function used in the computation of the mode                    */
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** univariate continuous distributions                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_cont_new( void )
     /*----------------------------------------------------------------------*/
     /* create a new (empty) distribution object                             */
     /* type: univariate continuous with given p.d.f.                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   none                                                               */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to distribution object                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  register struct unur_distr *distr;
  int i;

  /* get empty distribution object */
  distr = _unur_distr_generic_new();
  if (!distr) return NULL;

  /* set magic cookie */
  COOKIE_SET(distr,CK_DISTR_CONT);

  /* set type of distribution */
  distr->type = UNUR_DISTR_CONT;

  /* set id to generic distribution */
  distr->id = UNUR_DISTR_GENERIC;

  /* dimension of random vector */
  distr->dim = 1;   /* univariant */

  /* destructor */
  distr->destroy = _unur_distr_cont_free;

  /* clone */
  distr->clone = _unur_distr_cont_clone;

  /* set defaults                                                            */
  DISTR.pdf       = NULL;          /* pointer to PDF                         */
  DISTR.dpdf      = NULL;          /* pointer to derivative of PDF           */
  DISTR.logpdf    = NULL;          /* pointer to logPDF                      */
  DISTR.dlogpdf   = NULL;          /* pointer to derivative of logPDF        */
  DISTR.cdf       = NULL;          /* pointer to CDF                         */
  DISTR.logcdf    = NULL;          /* pointer to logCDF                      */
  DISTR.invcdf    = NULL;          /* pointer to inverse CDF                 */
  DISTR.hr        = NULL;          /* pointer to HR                          */

  DISTR.init      = NULL;          /* pointer to special init routine        */

  /* initialize parameters of the p.d.f.                                     */
  DISTR.n_params  = 0;             /* number of parameters of the pdf        */  
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++)
    DISTR.params[i] = 0.;

  /* initialize parameter vectors of the PDF                                 */
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++) {
    DISTR.n_param_vec[i] = 0;
    DISTR.param_vecs[i] = NULL;
  }  

  DISTR.norm_constant = 1.;        /* (log of) normalization constant for p.d.f.
				      (initialized to avoid accidently floating
				      point exception                        */

  DISTR.mode       = INFINITY;     /* location of mode (default: not known)  */
  DISTR.center     = 0.;           /* location of center (default: not given)*/
  DISTR.area       = 1.;           /* area below PDF (default: not known)    */

  DISTR.trunc[0] = DISTR.domain[0] = -INFINITY; /* left boundary of domain   */
  DISTR.trunc[1] = DISTR.domain[1] = INFINITY;  /* right boundary of domain  */

  DISTR.set_params = NULL;         /* funct for setting parameters and domain*/
  DISTR.upd_mode   = _unur_distr_cont_find_mode; /* funct for computing mode */
  DISTR.upd_area   = NULL;         /* funct for computing area               */

  DISTR.pdftree    = NULL;         /* pointer to function tree for PDF       */
  DISTR.dpdftree   = NULL;         /* pointer to function tree for dPDF      */
  DISTR.logpdftree = NULL;         /* pointer to function tree for logPDF    */
  DISTR.dlogpdftree = NULL;        /* pointer to function tree for dlogPDF   */
  DISTR.cdftree    = NULL;         /* pointer to function tree for CDF       */
  DISTR.logcdftree = NULL;         /* pointer to function tree for logCDF    */
  DISTR.hrtree     = NULL;         /* pointer to function tree for HR        */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_cont_new() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
_unur_distr_cont_clone( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* copy (clone) distribution object                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to source distribution object                    */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of distribution object                            */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
#define CLONE clone->data.cont

  struct unur_distr *clone;
  int i;

  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );

  /* allocate memory */
  clone = _unur_xmalloc( sizeof(struct unur_distr) );
  
  /* copy distribution object into clone */
  memcpy( clone, distr, sizeof( struct unur_distr ) );

  /* copy function trees into generator object (when there is one) */
  CLONE.pdftree  = (DISTR.pdftree)  ? _unur_fstr_dup_tree(DISTR.pdftree)  : NULL;
  CLONE.dpdftree = (DISTR.dpdftree) ? _unur_fstr_dup_tree(DISTR.dpdftree) : NULL;
  CLONE.logpdftree  = (DISTR.logpdftree)  ? _unur_fstr_dup_tree(DISTR.logpdftree)  : NULL;
  CLONE.dlogpdftree = (DISTR.dlogpdftree) ? _unur_fstr_dup_tree(DISTR.dlogpdftree) : NULL;
  CLONE.cdftree  = (DISTR.cdftree)  ? _unur_fstr_dup_tree(DISTR.cdftree)  : NULL;
  CLONE.logcdftree  = (DISTR.logcdftree)  ? _unur_fstr_dup_tree(DISTR.logcdftree)  : NULL;
  CLONE.hrtree   = (DISTR.hrtree)   ? _unur_fstr_dup_tree(DISTR.hrtree)   : NULL;
 
  /* clone of parameter arrays */  
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++) {
    CLONE.n_param_vec[i] = DISTR.n_param_vec[i];
    if (DISTR.param_vecs[i]) {
      CLONE.param_vecs[i] = _unur_xmalloc( DISTR.n_param_vec[i] * sizeof(double) );
      memcpy( CLONE.param_vecs[i], DISTR.param_vecs[i], DISTR.n_param_vec[i] * sizeof(double) );
    }
  }  

  /* Remark:
     The external object to which DISTR.extobj is pointing to is not
     copied. Only the pointer is copied.
     Thus is is in the reponsibility of the user of the 
     unur_distr_cont_set_extobj() call to handle this situation
     correctly!
  */
  
  /* copy user name for distribution */
  if (distr->name_str) {
    size_t len = strlen(distr->name_str) + 1;
    clone->name_str = _unur_xmalloc(len);
    memcpy( clone->name_str, distr->name_str, len );
    clone->name = clone->name_str;
  }
  
  /* for a derived distribution we also have to copy the underlying */
  /* distribution object                                            */
  if (distr->base != NULL) {
    clone->base = _unur_distr_clone(distr->base);
  }

  return clone;

#undef CLONE
} /* end of _unur_distr_cont_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_distr_cont_free( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* free distribution object                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* check arguments */
  if( distr == NULL ) /* nothing to do */
    return;
  _unur_check_distr_object( distr, CONT, RETURN_VOID );

  /* parameter arrays */
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++)
    if (DISTR.param_vecs[i]) free( DISTR.param_vecs[i] );
  
  /* function trees */
  if (DISTR.pdftree)  _unur_fstr_free(DISTR.pdftree);
  if (DISTR.dpdftree) _unur_fstr_free(DISTR.dpdftree);
  if (DISTR.logpdftree)  _unur_fstr_free(DISTR.logpdftree);
  if (DISTR.dlogpdftree) _unur_fstr_free(DISTR.dlogpdftree);
  if (DISTR.cdftree)  _unur_fstr_free(DISTR.cdftree);
  if (DISTR.logcdftree)  _unur_fstr_free(DISTR.logcdftree);
  if (DISTR.hrtree)   _unur_fstr_free(DISTR.hrtree);

  /* derived distribution */
  if (distr->base) _unur_distr_free(distr->base);

  /* user name for distribution */
  if (distr->name_str) free(distr->name_str);
  
  COOKIE_CLEAR(distr);
  free( distr );
} /* end of _unur_distr_cont_free() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_pdf( struct unur_distr *distr, UNUR_FUNCT_CONT *pdf )
     /*----------------------------------------------------------------------*/
     /* set PDF of distribution                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   pdf   ... pointer to PDF                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, pdf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  /* we do not allow overwriting a pdf */
  if (DISTR.pdf != NULL || DISTR.logpdf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of PDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.pdf = pdf;
  return UNUR_SUCCESS;

} /* end of unur_distr_cont_set_pdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_dpdf( struct unur_distr *distr, UNUR_FUNCT_CONT *dpdf )
     /*----------------------------------------------------------------------*/
     /* set derivative of PDF of distribution                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   dpdf   ... pointer to derivative of PDF                            */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, dpdf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  
  /* we do not allow overwriting a dpdf */
  if (DISTR.dpdf != NULL || DISTR.dlogpdf != NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of dPDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.dpdf = dpdf;
  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_dpdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_logpdf( struct unur_distr *distr, UNUR_FUNCT_CONT *logpdf )
     /*----------------------------------------------------------------------*/
     /* set logPDF of distribution                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   logpdf ... pointer to logPDF                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, logpdf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  /* we do not allow overwriting a pdf */
  if (DISTR.pdf != NULL || DISTR.logpdf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of logPDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.logpdf = logpdf;
  DISTR.pdf = _unur_distr_cont_eval_pdf_from_logpdf;

  return UNUR_SUCCESS;

} /* end of unur_distr_cont_set_logpdf() */

/*---------------------------------------------------------------------------*/

double
_unur_distr_cont_eval_pdf_from_logpdf( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate PDF of distribution at x                                    */
     /* wrapper when only logPDF is given                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for pdf                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   PDF(x)                                                             */
     /*----------------------------------------------------------------------*/
{
  if (DISTR.logpdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return exp(_unur_cont_logPDF(x,distr));
} /* end of _unur_distr_cont_eval_pdf_from_logpdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_dlogpdf( struct unur_distr *distr, UNUR_FUNCT_CONT *dlogpdf )
     /*----------------------------------------------------------------------*/
     /* set derivative of logPDF of distribution                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr   ... pointer to distribution object                         */
     /*   dlogpdf ... pointer to derivative of logPDF                        */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, dlogpdf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  
  /* we do not allow overwriting a dpdf */
  if (DISTR.dpdf != NULL || DISTR.dlogpdf != NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of dlogPDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.dlogpdf = dlogpdf;
  DISTR.dpdf = _unur_distr_cont_eval_dpdf_from_dlogpdf;

  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_dlogpdf() */

/*---------------------------------------------------------------------------*/

double
_unur_distr_cont_eval_dpdf_from_dlogpdf( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate derivative of PDF of distribution at x                      */
     /* wrapper when only derivative of logPDF is given                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x      ... argument for dPDF                                       */
     /*   distr  ... pointer to distribution object                          */
     /*                                                                      */
     /* return:                                                              */
     /*   dPDF(x)   ... on success                                           */
     /*   INFINITY  ... on error                                             */
     /*----------------------------------------------------------------------*/
{
  if (DISTR.logpdf == NULL || DISTR.dlogpdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return exp(_unur_cont_logPDF(x,distr)) * _unur_cont_dlogPDF(x,distr);
} /* end of _unur_distr_cont_eval_dpdf_from_dlogpdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_cdf( struct unur_distr *distr, UNUR_FUNCT_CONT *cdf )
     /*----------------------------------------------------------------------*/
     /* set CDF of distribution                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   cdf   ... pointer to CDF                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, cdf,UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  
  /* we do not allow overwriting a cdf */
  if (DISTR.cdf != NULL || DISTR.logcdf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of CDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.cdf = cdf;
  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_cdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_invcdf( struct unur_distr *distr, UNUR_FUNCT_CONT *invcdf )
     /*----------------------------------------------------------------------*/
     /* set inverse CDF of distribution                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   invcdf ... pointer to inverse CDF                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, invcdf,UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  
  /* we do not allow overwriting an inverse cdf */
  if (DISTR.invcdf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of inverse CDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.invcdf = invcdf;
  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_invcdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_logcdf( struct unur_distr *distr, UNUR_FUNCT_CONT *logcdf )
     /*----------------------------------------------------------------------*/
     /* set logCDF of distribution                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   logcdf ... pointer to logCDF                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, logcdf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  /* we do not allow overwriting a cdf */
  if (DISTR.cdf != NULL || DISTR.logcdf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of logCDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.logcdf = logcdf;
  DISTR.cdf = _unur_distr_cont_eval_cdf_from_logcdf;

  return UNUR_SUCCESS;

} /* end of unur_distr_cont_set_logcdf() */

/*---------------------------------------------------------------------------*/

double
_unur_distr_cont_eval_cdf_from_logcdf( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate CDF of distribution at x                                    */
     /* wrapper when only logCDF is given                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for cdf                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   CDF(x)                                                             */
     /*----------------------------------------------------------------------*/
{
  if (DISTR.logcdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return exp(_unur_cont_logCDF(x,distr));
} /* end of _unur_distr_cont_eval_cdf_from_logcdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_hr( struct unur_distr *distr, UNUR_FUNCT_CONT *hr )
     /*----------------------------------------------------------------------*/
     /* set HR of distribution                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   hr    ... pointer to HR                                            */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, hr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  
  /* we do not allow overwriting a cdf */
  if (DISTR.hr != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of HR not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.hr = hr;
  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_hr() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_pdfstr( struct unur_distr *distr, const char *pdfstr )
     /*----------------------------------------------------------------------*/
     /* set PDF and its derivative of distribution via a string interface    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   pdfstr ... string that describes function term of PDF              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, pdfstr, UNUR_ERR_NULL );

  /* When the user calls unur_distr_cont_set_cdfstr() before this function   */
  /* then the PDF and dPDF probably are already set and we get an error.     */
  /* This might happen in particular for the String API.                     */
  /* To avoid confusion we remove these settings.                            */
  if ( DISTR.pdftree || DISTR.logpdftree ) {
    if (DISTR.pdftree)  _unur_fstr_free(DISTR.pdftree);
    if (DISTR.dpdftree) _unur_fstr_free(DISTR.dpdftree);
    if (DISTR.logpdftree)  _unur_fstr_free(DISTR.logpdftree);
    if (DISTR.dlogpdftree) _unur_fstr_free(DISTR.dlogpdftree);
    DISTR.pdf = NULL;
    DISTR.dpdf = NULL;
    DISTR.logpdf = NULL;
    DISTR.dlogpdf = NULL;
  }

  /* we do not allow overwriting a PDF */
  if (DISTR.pdf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of PDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* parse PDF string */
  if ( (DISTR.pdftree = _unur_fstr2tree(pdfstr)) == NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.pdf  = _unur_distr_cont_eval_pdf_tree;

  /* make derivative */
  if ( (DISTR.dpdftree = _unur_fstr_make_derivative(DISTR.pdftree)) == NULL )
    return UNUR_ERR_DISTR_DATA;
  DISTR.dpdf = _unur_distr_cont_eval_dpdf_tree;

  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_pdfstr() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_logpdfstr( struct unur_distr *distr, const char *logpdfstr )
     /*----------------------------------------------------------------------*/
     /* set logPDF and its derivative of distribution via a string interface */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr     ... pointer to distribution object                       */
     /*   logpdfstr ... string that describes function term of logPDF        */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, logpdfstr, UNUR_ERR_NULL );

  /* When the user calls unur_distr_cont_set_cdfstr() before this function   */
  /* then the PDF and dPDF probably are already set and we get an error.     */
  /* This might happen in particular for the String API.                     */
  /* To avoid confusion we remove these settings.                            */
  if ( DISTR.pdftree || DISTR.logpdftree ) {
    if (DISTR.pdftree)  _unur_fstr_free(DISTR.pdftree);
    if (DISTR.dpdftree) _unur_fstr_free(DISTR.dpdftree);
    if (DISTR.logpdftree)  _unur_fstr_free(DISTR.logpdftree);
    if (DISTR.dlogpdftree) _unur_fstr_free(DISTR.dlogpdftree);
    DISTR.pdf = NULL;
    DISTR.dpdf = NULL;
    DISTR.logpdf = NULL;
    DISTR.dlogpdf = NULL;
  }

  /* we do not allow overwriting a PDF */
  if (DISTR.pdf != NULL || DISTR.logpdf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of logPDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* parse logPDF string */
  if ( (DISTR.logpdftree = _unur_fstr2tree(logpdfstr)) == NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.logpdf  = _unur_distr_cont_eval_logpdf_tree;
  DISTR.pdf = _unur_distr_cont_eval_pdf_from_logpdf;

  /* make derivative */
  if ( (DISTR.dlogpdftree = _unur_fstr_make_derivative(DISTR.logpdftree)) == NULL )
    return UNUR_ERR_DISTR_DATA;
  DISTR.dlogpdf = _unur_distr_cont_eval_dlogpdf_tree;
  DISTR.dpdf = _unur_distr_cont_eval_dpdf_from_dlogpdf;

  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_logpdfstr() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_cdfstr( struct unur_distr *distr, const char *cdfstr )
     /*----------------------------------------------------------------------*/
     /* set CDF of distribution via a string interface                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   cdfstr ... string that describes function term of CDF              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, cdfstr, UNUR_ERR_NULL );

  /* we do not allow overwriting a CDF */
  if (DISTR.cdf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of CDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* parse string */
  if ( (DISTR.cdftree = _unur_fstr2tree(cdfstr)) == NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.cdf  = _unur_distr_cont_eval_cdf_tree;
  
  /* make derivatives (if necessary) */
  /* We do not compute derivatives, if these strings are alreasy set by an   */
  /* unur_distr_cont_set_pdfstr() or unur_distr_cont_set_logpdfstr() call.   */
  /* We also do not return an error code if computation of derivatives       */
  /* did not work.                                                           */
  if (DISTR.pdftree == NULL)
    if ( (DISTR.pdftree = _unur_fstr_make_derivative(DISTR.cdftree)) != NULL )
      DISTR.pdf = _unur_distr_cont_eval_pdf_tree;
  if (DISTR.dpdftree == NULL)
    if ( (DISTR.dpdftree = _unur_fstr_make_derivative(DISTR.pdftree)) != NULL )
      DISTR.dpdf = _unur_distr_cont_eval_dpdf_tree;

  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_cdfstr() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_logcdfstr( struct unur_distr *distr, const char *logcdfstr )
     /*----------------------------------------------------------------------*/
     /* set logCDF and its derivative of distribution via a string interface */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr     ... pointer to distribution object                       */
     /*   logcdfstr ... string that describes function term of logCDF        */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, logcdfstr, UNUR_ERR_NULL );

  /* we do not allow overwriting a CDF */
  if (DISTR.cdf != NULL || DISTR.logcdf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of logCDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* parse logCDF string */
  if ( (DISTR.logcdftree = _unur_fstr2tree(logcdfstr)) == NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.logcdf  = _unur_distr_cont_eval_logcdf_tree;
  DISTR.cdf = _unur_distr_cont_eval_cdf_from_logcdf;

  /* [ we do not make derivatives ] */

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_logcdfstr() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_hrstr( struct unur_distr *distr, const char *hrstr )
     /*----------------------------------------------------------------------*/
     /* set HR of distribution via a string interface                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   hrstr  ... string that describes function term of CDF              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, hrstr, UNUR_ERR_NULL );

  /* we do not allow overwriting a CDF */
  if (DISTR.hr != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of CDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* parse string */
  if ( (DISTR.hrtree = _unur_fstr2tree(hrstr)) == NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;
  }

  /* set evaluation function */
  DISTR.hr  = _unur_distr_cont_eval_hr_tree;

  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_hrstr() */

/*---------------------------------------------------------------------------*/

double
_unur_distr_cont_eval_pdf_tree( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate function tree for PDF.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for PDF                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   PDF at x                                                           */
     /*----------------------------------------------------------------------*/
{
  return ((DISTR.pdftree) ? _unur_fstr_eval_tree(DISTR.pdftree,x) : INFINITY);
} /* end of _unur_distr_cont_eval_pdf_tree() */

/*---------------------------------------------------------------------------*/

double
_unur_distr_cont_eval_logpdf_tree( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate function tree for logPDF.                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for logPDF                                      */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   logPDF at x                                                        */
     /*----------------------------------------------------------------------*/
{
  return ((DISTR.logpdftree) ? _unur_fstr_eval_tree(DISTR.logpdftree,x) : INFINITY);
} /* end of _unur_distr_cont_eval_logpdf_tree() */

/*---------------------------------------------------------------------------*/

double
_unur_distr_cont_eval_dpdf_tree( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate function tree for derivative of PDF.                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for derivative of PDF                           */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   derivative of PDF at x                                             */
     /*----------------------------------------------------------------------*/
{
  return ((DISTR.dpdftree) ? _unur_fstr_eval_tree(DISTR.dpdftree,x) : INFINITY);
} /* end of _unur_distr_cont_eval_dpdf_tree() */

/*---------------------------------------------------------------------------*/


double
_unur_distr_cont_eval_dlogpdf_tree( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate function tree for derivative of logPDF.                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for derivative of logPDF                        */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   derivative of logPDF at x                                          */
     /*----------------------------------------------------------------------*/
{
  return ((DISTR.dlogpdftree) ? _unur_fstr_eval_tree(DISTR.dlogpdftree,x) : INFINITY);
} /* end of _unur_distr_cont_eval_dpdf_tree() */

/*---------------------------------------------------------------------------*/

double
_unur_distr_cont_eval_cdf_tree( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate function tree for CDF.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for CDF                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   CDF at x                                                           */
     /*----------------------------------------------------------------------*/
{
  return ((DISTR.cdftree) ? _unur_fstr_eval_tree(DISTR.cdftree,x) : INFINITY);
} /* end of _unur_distr_cont_eval_cdf_tree() */

/*---------------------------------------------------------------------------*/

double
_unur_distr_cont_eval_logcdf_tree( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate function tree for logCDF.                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for logCDF                                      */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   logCDF at x                                                        */
     /*----------------------------------------------------------------------*/
{
  return ((DISTR.logcdftree) ? _unur_fstr_eval_tree(DISTR.logcdftree,x) : INFINITY);
} /* end of _unur_distr_cont_eval_logcdf_tree() */

/*---------------------------------------------------------------------------*/

double
_unur_distr_cont_eval_hr_tree( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate function tree for HR.                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for CDF                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   CDF at x                                                           */
     /*----------------------------------------------------------------------*/
{
  return ((DISTR.hrtree) ? _unur_fstr_eval_tree(DISTR.hrtree,x) : INFINITY);
} /* end of _unur_distr_cont_eval_hr_tree() */

/*---------------------------------------------------------------------------*/

char *
unur_distr_cont_get_pdfstr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get PDF string that is given via the string interface                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to resulting string.                                       */
     /*                                                                      */
     /* comment:                                                             */
     /*   This string should be freed when it is not used any more.          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  _unur_check_NULL( NULL, DISTR.pdftree, NULL );

  /* make and return string */
  return _unur_fstr_tree2string(DISTR.pdftree,"x","PDF",TRUE);
} /* end of unur_distr_cont_get_pdfstr() */

/*---------------------------------------------------------------------------*/

char *
unur_distr_cont_get_dpdfstr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get string for derivative of PDF that is given via string interface  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to resulting string.                                       */
     /*                                                                      */
     /* comment:                                                             */
     /*   This string should be freed when it is not used any more.          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  _unur_check_NULL( NULL, DISTR.dpdftree, NULL );

  /* make and return string */
  return _unur_fstr_tree2string(DISTR.dpdftree,"x","dPDF",TRUE);
} /* end of unur_distr_cont_get_dpdfstr() */

/*---------------------------------------------------------------------------*/

char *
unur_distr_cont_get_logpdfstr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get logPDF string that is given via the string interface             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to resulting string.                                       */
     /*                                                                      */
     /* comment:                                                             */
     /*   This string should be freed when it is not used any more.          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  _unur_check_NULL( NULL, DISTR.logpdftree, NULL );

  /* make and return string */
  return _unur_fstr_tree2string(DISTR.logpdftree,"x","logPDF",TRUE);
} /* end of unur_distr_cont_get_logpdfstr() */

/*---------------------------------------------------------------------------*/

char *
unur_distr_cont_get_dlogpdfstr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get string for derivative of logPDF that is given via string API     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to resulting string.                                       */
     /*                                                                      */
     /* comment:                                                             */
     /*   This string should be freed when it is not used any more.          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  _unur_check_NULL( NULL, DISTR.dlogpdftree, NULL );

  /* make and return string */
  return _unur_fstr_tree2string(DISTR.dlogpdftree,"x","dlogPDF",TRUE);
} /* end of unur_distr_cont_get_dlogpdfstr() */

/*---------------------------------------------------------------------------*/

char *
unur_distr_cont_get_cdfstr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get CDF string that is given via the string interface                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to resulting string.                                       */
     /*                                                                      */
     /* comment:                                                             */
     /*   This string should be freed when it is not used any more.          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  _unur_check_NULL( NULL, DISTR.cdftree, NULL );

  /* make and return string */
  return _unur_fstr_tree2string(DISTR.cdftree,"x","CDF",TRUE);
} /* end of unur_distr_cont_get_cdfstr() */

/*---------------------------------------------------------------------------*/

char *
unur_distr_cont_get_logcdfstr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get logCDF string that is given via the string interface             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to resulting string.                                       */
     /*                                                                      */
     /* comment:                                                             */
     /*   This string should be freed when it is not used any more.          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  _unur_check_NULL( NULL, DISTR.logcdftree, NULL );

  /* make and return string */
  return _unur_fstr_tree2string(DISTR.logcdftree,"x","logCDF",TRUE);
} /* end of unur_distr_cont_get_logcdfstr() */

/*---------------------------------------------------------------------------*/

char *
unur_distr_cont_get_hrstr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get HR string that is given via the string interface                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to resulting string.                                       */
     /*                                                                      */
     /* comment:                                                             */
     /*   This string should be freed when it is not used any more.          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );
  _unur_check_NULL( NULL, DISTR.hrtree, NULL );

  /* make and return string */
  return _unur_fstr_tree2string(DISTR.hrtree,"x","HR",TRUE);
} /* end of unur_distr_cont_get_hrstr() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_CONT *
unur_distr_cont_get_pdf( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to PDF of distribution                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to PDF                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );

  return DISTR.pdf;
} /* end of unur_distr_cont_get_pdf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_CONT *
unur_distr_cont_get_dpdf( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to derivative of PDF of distribution                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to derivative of PDF                                       */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );

  return DISTR.dpdf;
} /* end of unur_distr_cont_get_dpdf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_CONT *
unur_distr_cont_get_logpdf( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to logPDF of distribution                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to logPDF                                                  */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );

  return DISTR.logpdf;
} /* end of unur_distr_cont_get_logpdf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_CONT *
unur_distr_cont_get_dlogpdf( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to derivative of logPDF of distribution                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to derivative of logPDF                                       */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );

  return DISTR.dlogpdf;
} /* end of unur_distr_cont_get_dlogpdf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_CONT *
unur_distr_cont_get_cdf( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to CDF of distribution                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to CDF                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );

  return DISTR.cdf;
} /* end of unur_distr_cont_get_cdf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_CONT *
unur_distr_cont_get_invcdf( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to inverse CDF of distribution                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to inverse CDF                                             */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );

  return DISTR.invcdf;
} /* end of unur_distr_cont_get_invcdf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_CONT *
unur_distr_cont_get_logcdf( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to logCDF of distribution                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to logCDF                                                  */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );

  return DISTR.logcdf;
} /* end of unur_distr_cont_get_logcdf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_CONT *
unur_distr_cont_get_hr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to HR of distribution                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to HR                                                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );

  return DISTR.hr;
} /* end of unur_distr_cont_get_hr() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cont_eval_pdf( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate PDF of distribution at x                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for pdf                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pdf(x)                                                             */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  if (DISTR.pdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_cont_PDF(x,distr);
} /* end of unur_distr_cont_eval_pdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cont_eval_dpdf( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate derivative of PDF of distribution at x                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for dpdf                                        */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   (pdf(x))'                                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  if (DISTR.dpdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_cont_dPDF(x,distr);
} /* end of unur_distr_cont_eval_dpdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cont_eval_logpdf( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate logPDF of distribution at x                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for logPDF                                      */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   logPDF(x)                                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  if (DISTR.logpdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_cont_logPDF(x,distr);
} /* end of unur_distr_cont_eval_logpdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cont_eval_dlogpdf( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate derivative of logPDF of distribution at x                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for dlogPDF                                     */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   (pdf(x))'                                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  if (DISTR.dlogpdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_cont_dlogPDF(x,distr);
} /* end of unur_distr_cont_eval_dlogpdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cont_eval_cdf( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate CDF of distribution at x                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for CDF                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   cdf(x)                                                             */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  if (DISTR.cdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_cont_CDF(x,distr);
} /* end of unur_distr_cont_eval_cdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cont_eval_invcdf( double u, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate inverse CDF of distribution at u                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   u     ... argument for inverse CDF                                 */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   invcdf(u)                                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  if (DISTR.invcdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  if (u<=0.)
    return DISTR.domain[0];
  if (u>=1.)
    return DISTR.domain[1];
  else
    return _unur_cont_invCDF(u,distr);

} /* end of unur_distr_cont_eval_invcdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cont_eval_logcdf( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate logCDF of distribution at x                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for logCDF                                      */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   logCDF(x)                                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  if (DISTR.logcdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_cont_logCDF(x,distr);
} /* end of unur_distr_cont_eval_logcdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cont_eval_hr( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate HR of distribution at x                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for cdf                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   hr(x)                                                              */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  if (DISTR.hr == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_cont_HR(x,distr);
} /* end of unur_distr_cont_eval_hr() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_pdfparams( struct unur_distr *distr, const double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* set array of parameters for distribution                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   params   ... list of arguments                                     */
     /*   n_params ... number of arguments                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  if (n_params>0) _unur_check_NULL(distr->name,params,UNUR_ERR_NULL);

  /* first check number of new parameter for the distribution */
  if (n_params < 0 || n_params > UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"");
    return UNUR_ERR_DISTR_NPARAMS;
  }

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* even if the set routine fails, the derived parameters are
     marked as unknown. but this is o.k. since in this case something
     has been wrong. */

  /* use special routine for setting parameters
     (if there is one) */

  if (distr->base && BASE.set_params) 
    return (BASE.set_params(distr->base,params,n_params));

  if (DISTR.set_params)
    return (DISTR.set_params(distr,params,n_params));

  /* otherwise simply copy parameters */

  if (distr->base) {
    BASE.n_params = n_params;
    if (n_params) memcpy( BASE.params, params, n_params*sizeof(double) );
  }

  else {
    DISTR.n_params = n_params;
    if (n_params) memcpy( DISTR.params, params, n_params*sizeof(double) );
  }

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_pdfparams() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_get_pdfparams( const struct unur_distr *distr, const double **params )
     /*----------------------------------------------------------------------*/
     /* get number of pdf parameters and sets pointer to array params[] of   */
     /* parameters                                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   params   ... pointer to list of arguments                          */
     /*                                                                      */
     /* return:                                                              */
     /*   number of pdf parameters                                           */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CONT, 0 );

  if (distr->base) {
    /* for derived distributions (e.g. order statistics)
       the parameters for the underlying distributions are returned */
    *params = (BASE.n_params) ? BASE.params : NULL;
    return BASE.n_params;
  }
  else {
    *params = (DISTR.n_params) ? DISTR.params : NULL;
    return DISTR.n_params;
  }

} /* end of unur_distr_cont_get_pdfparams() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_pdfparams_vec( struct unur_distr *distr, int par, const double *param_vec, int n_param_vec )
     /*----------------------------------------------------------------------*/
     /* set vector array parameters for distribution                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   par      ... which parameter is set                                */
     /*   param_vec   ... parameter array with number `par'                  */
     /*   n_param_vec ... length of parameter array                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  /* check new parameter for distribution */
  if (par < 0 || par >= UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"invalid parameter position");
    return UNUR_ERR_DISTR_NPARAMS;
  }

  if (param_vec != NULL) {
    /* allocate memory */
    DISTR.param_vecs[par] = _unur_xrealloc( DISTR.param_vecs[par], n_param_vec * sizeof(double) );
    /* copy parameters */
    memcpy( DISTR.param_vecs[par], param_vec, n_param_vec*sizeof(double) );
    /* set length of array */
    DISTR.n_param_vec[par] = n_param_vec;
  }
  else {
    if (DISTR.param_vecs[par]) free(DISTR.param_vecs[par]);
    DISTR.param_vecs[par] = NULL;
    DISTR.n_param_vec[par] = 0;
  }

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_pdfparams_vec() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_get_pdfparams_vec( const struct unur_distr *distr, int par, const double **param_vecs )
     /*----------------------------------------------------------------------*/
     /* get number of PDF parameters and sets pointer to array params[] of   */
     /* parameters                                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   par      ... which parameter is read                               */
     /*   params   ... pointer to parameter array with number `par'          */
     /*                                                                      */
     /* return:                                                              */
     /*   length of parameter array with number `par'                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CONT, 0 );

  /* check new parameter for distribution */
  if (par < 0 || par >= UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"invalid parameter position");
    *param_vecs = NULL;
    return 0;
  }
  
  *param_vecs = DISTR.param_vecs[par];

  return (*param_vecs) ? DISTR.n_param_vec[par] : 0;
} /* end of unur_distr_cont_get_pdfparams_vec() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_domain( struct unur_distr *distr, double left, double right )
     /*----------------------------------------------------------------------*/
     /* set the left and right borders of the domain of the distribution     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   the new boundary points may be +/- INFINITY                        */
     /*----------------------------------------------------------------------*/
{
  int unsigned is_set = 0u;

  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  /* check new parameter for distribution */
  if (left >= right) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }

  /* we have to deal with the mode */
  if ( distr->set & UNUR_DISTR_SET_MODE ) {
    is_set |= UNUR_DISTR_SET_MODE;
    /* mode has been set and new domain might be subset of old domain */
    /* we assume the density is unimodal and thus monotone on
       either side of the mode!! */
    if ( DISTR.mode < left)       DISTR.mode = left;
    else if ( DISTR.mode > right) DISTR.mode = right;
  }

  /* we have to deal with the mode */
  if ( distr->set & UNUR_DISTR_SET_CENTER ) {
    is_set |= UNUR_DISTR_SET_CENTER;
    if ( DISTR.center < left)       DISTR.center = left;
    else if ( DISTR.center > right) DISTR.center = right;
  }

  /* store data */
  DISTR.trunc[0] = DISTR.domain[0] = left;
  DISTR.trunc[1] = DISTR.domain[1] = right;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_DOMAIN;

  /* if distr is an object for a standard distribution, this     */
  /* is not the original domain of it. (not a "standard domain") */
  /* However, since we have changed the domain, we assume        */
  /* that this is not a truncated distribution.                  */
  /* At last we have to mark all derived parameters as unknown.  */
  distr->set &= ~(UNUR_DISTR_SET_STDDOMAIN |
		  UNUR_DISTR_SET_TRUNCATED | 
		  UNUR_DISTR_SET_MASK_DERIVED );
  distr->set |= is_set;

  if (distr->base) {
    /* for derived distributions (e.g. order statistics)
       we also set the domain for the underlying distribution */
    BASE.trunc[0] = BASE.domain[0] = left;
    BASE.trunc[1] = BASE.domain[1] = right;
    distr->base->set &= ~(UNUR_DISTR_SET_STDDOMAIN |
			  UNUR_DISTR_SET_TRUNCATED | 
			  UNUR_DISTR_SET_MASK_DERIVED );
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_distr_cont_set_domain() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_get_domain( const struct unur_distr *distr, double *left, double *right )
     /*----------------------------------------------------------------------*/
     /* set the left and right borders of the domain of the distribution     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   if no boundaries have been set +/- INFINITY is returned.           */
     /*----------------------------------------------------------------------*/
{
  /* in case of error the boundaries are set to +/- INFINITY */
  *left = -INFINITY;
  *right = INFINITY;

  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  /* o.k. */
  *left  = DISTR.domain[0];
  *right = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of unur_distr_cont_get_domain() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_get_truncated( const struct unur_distr *distr, double *left, double *right )
     /*----------------------------------------------------------------------*/
     /* set the left and right borders of the truncated domain of the        */
     /* distribution                                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   if no boundaries have been set +/- INFINITY is returned.           */
     /*----------------------------------------------------------------------*/
{
  /* in case of error the boundaries are set to +/- INFINITY */
  *left = -INFINITY;
  *right = INFINITY;

  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  /* o.k. */
  *left  = (distr->set & UNUR_DISTR_SET_TRUNCATED) ? DISTR.trunc[0] : DISTR.domain[0];
  *right = (distr->set & UNUR_DISTR_SET_TRUNCATED) ? DISTR.trunc[1] : DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of unur_distr_cont_get_truncated() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_mode( struct unur_distr *distr, double mode )
     /*----------------------------------------------------------------------*/
     /* set mode of distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   mode  ... mode of p.d.f.                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  /* check whether mode is inside the domain */
  if (mode < DISTR.domain[0] || mode > DISTR.domain[1]) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"mode not in domain");
    return UNUR_ERR_DISTR_SET;
  }

  /* store data */
  DISTR.mode = mode;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_MODE;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_mode() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_cont_upd_mode( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* (re-) compute mode of distribution (if possible)                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  if (DISTR.upd_mode == NULL) {
    /* no function to compute mode available */
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }

  /* compute mode */
  if ((DISTR.upd_mode)(distr)==UNUR_SUCCESS) {
    /* changelog */
    distr->set |= UNUR_DISTR_SET_MODE;
    return UNUR_SUCCESS;
  }
  else {
    /* computing of mode failed */
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }

} /* end of unur_distr_cont_upd_mode() */
  
/*---------------------------------------------------------------------------*/

double
unur_distr_cont_get_mode( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get mode of distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   mode of distribution                                               */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  /* mode known ? */
  if ( !(distr->set & UNUR_DISTR_SET_MODE) ) {
    /* try to compute mode */
    if (DISTR.upd_mode == NULL) {
      /* no function to compute mode available */
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"mode");
      return INFINITY;
    }
    else {
      /* compute mode */
      if (unur_distr_cont_upd_mode(distr)!=UNUR_SUCCESS) {
	/* finding mode not successfully */
	_unur_error(distr->name,UNUR_ERR_DISTR_GET,"mode");
	return INFINITY;
      }
    }
  }

  return DISTR.mode;

} /* end of unur_distr_cont_get_mode() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_center( struct unur_distr *distr, double center )
     /*----------------------------------------------------------------------*/
     /* set center of distribution                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   center ... center of PDF                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  DISTR.center = center;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_CENTER;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_center() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cont_get_center( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get center of distribution                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   center of distribution                                             */
     /*----------------------------------------------------------------------*/
{
  double center;

  /* check arguments */
  _unur_check_NULL( NULL, distr, 0. );
  _unur_check_distr_object( distr, CONT, 0. );

  /* center given */
  if ( distr->set & UNUR_DISTR_SET_CENTER )
    center = DISTR.center;

  /* else try mode */
  else if ( distr->set & UNUR_DISTR_SET_MODE ) 
    center = DISTR.mode;

  /* else unknown: use default */
  else
    center = 0.;

  /* return result 0 */
  return center;

} /* end of unur_distr_cont_get_center() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_pdfarea( struct unur_distr *distr, double area )
     /*----------------------------------------------------------------------*/
     /* set area below p.d.f.                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   area  ... area below p.d.f.                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  /* check new parameter for distribution */
  if (area <= 0.) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"pdf area <= 0");
    return UNUR_ERR_DISTR_SET;
  }

  DISTR.area = area;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_PDFAREA;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_distr_cont_set_pdfarea() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_cont_upd_pdfarea( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* (re-) compute area below p.d.f. of distribution (if possible)        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  if (DISTR.upd_area == NULL) {
    /* no function to compute mode available */
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }

  /* compute area */
  if (((DISTR.upd_area)(distr)!=UNUR_SUCCESS) || DISTR.area <= 0.) {
    /* computing of area failed */
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"upd area <= 0");
    DISTR.area = 1.;   /* avoid possible floating point exceptions */
    distr->set &= ~UNUR_DISTR_SET_PDFAREA;
    return UNUR_ERR_DISTR_SET;
  }

  /* changelog */
  distr->set |= UNUR_DISTR_SET_PDFAREA;

  return UNUR_SUCCESS;
} /* end of unur_distr_cont_upd_pdfarea() */
  
/*---------------------------------------------------------------------------*/

double
unur_distr_cont_get_pdfarea( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get area below p.d.f. of distribution                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   area below p.d.f. of distribution                                  */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  /* area known ? */
  if ( !(distr->set & UNUR_DISTR_SET_PDFAREA) ) {
    /* try to compute area */
    if ( unur_distr_cont_upd_pdfarea(distr) != UNUR_SUCCESS ) {
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"area");
      return INFINITY;
    }
  }

  return DISTR.area;

} /* end of unur_distr_cont_get_pdfarea() */

/*****************************************************************************/

double
_unur_aux_pdf(double x, void *p) 
     /*----------------------------------------------------------------------*/
     /* Auxiliary function used in the computation of the mode               */
     /*----------------------------------------------------------------------*/	
{
  struct unur_distr *distr = p;
  return (DISTR.pdf(x, distr));
} 

/*---------------------------------------------------------------------------*/

int
_unur_distr_cont_find_mode( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
{
  struct unur_funct_generic pdf;  /* density function to be maximized */
  double mode;

  /* check arguments */
  CHECK_NULL( distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  if (DISTR.pdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_GET,"PDF required for finding mode numerically"); 
    return UNUR_ERR_DISTR_DATA;
  }

  /* set parameters of pdf */
  pdf.f = _unur_aux_pdf;
  pdf.params = distr;

  /* compute mode; use DISTR.center as first guess */
  mode = _unur_util_find_max( pdf, DISTR.domain[0], DISTR.domain[1], DISTR.center );

  /* check result */
  if (_unur_isfinite(mode)){
    /* mode successfully computed */
    DISTR.mode = mode;
    /* changelog */
    distr->set |= UNUR_DISTR_SET_MODE | UNUR_DISTR_SET_MODE_APPROX ; 
    /* o.k. */
    return UNUR_SUCCESS;
  }

  else {
    /* computing mode did not work */
    /* (we do not change mode entry in distribution object) */
    return UNUR_ERR_DISTR_DATA;
  }

} /* end of _unur_distr_cont_find_mode() */

/*---------------------------------------------------------------------------*/

int
_unur_distr_cont_find_center( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* search for an appropriate point for center.                          */
     /* if such a point is found, then it is stored in 'distr'.              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
#define PDF(x)     _unur_cont_PDF((x),(distr))     /* call to PDF */
#define logPDF(x)  _unur_cont_logPDF((x),(distr))  /* call to logPDF */

  double center, fc;   /* given center and PDF at center */
  double x, fx;        /* auxiliary point and PDF at point */
  int i,d;

  /* given center of distribution */
  center = DISTR.center;
  fc = (DISTR.logpdf!=NULL) ? exp(logPDF(center)) : PDF(center);

  /* check given center */
  if (fc>0. && _unur_isfinite(fc)) return UNUR_SUCCESS;

  /* search */
  for (d=0; d<2; d++) {
    /* left and right boundary of truncated domain, resp. */
    x = DISTR.trunc[d];

    /* Do x and the given center differ? */
    if (_unur_FP_equal(center,x))
      continue;

    /* search from boundary towards given center */
    for (i=0; i<50; i++) {
      x = _unur_arcmean(x,center);
      fx = (DISTR.logpdf!=NULL) ? exp(logPDF(x)) : PDF(x);

      if (fx>0. && _unur_isfinite(fx)) {
	/* successfull */
  	DISTR.center = x;
	/* changelog */
	distr->set |= UNUR_DISTR_SET_CENTER | UNUR_DISTR_SET_CENTER_APPROX ;
	/* o.k. */
  	return UNUR_SUCCESS;
      }
    }
  }

  /* could not find appropriate point */
  return UNUR_FAILURE;

#undef PDF
#undef logPDF
} /* end of _unur_distr_cont_find_center() */


/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_distr_cont_debug( const struct unur_distr *distr, const char *genid )
     /*----------------------------------------------------------------------*/
     /* write info about distribution into LOG file                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   genid ... pointer to generator id                                  */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  int i;

  /* check arguments */
  CHECK_NULL(distr,RETURN_VOID);
  COOKIE_CHECK(distr,CK_DISTR_CONT,RETURN_VOID);

  LOG = unur_get_stream();

  /* is this a derived distribution */
  if (distr->base) {
    switch (distr->id) {
    case UNUR_DISTR_CORDER:
      _unur_distr_corder_debug(distr,genid);
      return;
    case UNUR_DISTR_CXTRANS:
      _unur_distr_cxtrans_debug(distr,genid);
      return;
    case UNUR_DISTR_CONDI:
      _unur_distr_condi_debug(distr,genid);
      return;
    default:
      /* nothing to do */
      fprintf(LOG,"%s: derived distribution.\n",genid);
      fprintf(LOG,"%s:\n",genid);
    }
  }

  fprintf(LOG,"%s: distribution:\n",genid);
  fprintf(LOG,"%s:\ttype = continuous univariate distribution\n",genid);
  fprintf(LOG,"%s:\tname = %s\n",genid,distr->name);

  fprintf(LOG,"%s:\tPDF with %d argument(s)\n",genid,DISTR.n_params);
  for( i=0; i<DISTR.n_params; i++ )
      fprintf(LOG,"%s:\t\tparam[%d] = %g\n",genid,i,DISTR.params[i]);

  fprintf(LOG,"%s:\tfunctions: ",genid);
  if (DISTR.cdf) fprintf(LOG,"CDF ");
  if (DISTR.logcdf) fprintf(LOG,"logCDF ");
  if (DISTR.pdf) fprintf(LOG,"PDF ");
  if (DISTR.logpdf) fprintf(LOG,"logPDF ");
  if (DISTR.dpdf) fprintf(LOG,"dPDF ");
  if (DISTR.dlogpdf) fprintf(LOG,"dlogPDF ");
  if (DISTR.hr) fprintf(LOG,"HR ");
  fprintf(LOG,"\n");

  if (distr->set & UNUR_DISTR_SET_MODE)
    fprintf(LOG,"%s:\tmode = %g\n",genid,DISTR.mode);
  else
    fprintf(LOG,"%s:\tmode unknown\n",genid);
  
  fprintf(LOG,"%s:\tcenter = ",genid);
  if (distr->set & UNUR_DISTR_SET_CENTER)
    fprintf(LOG,"%g%s\n", DISTR.center,
	    (distr->set & UNUR_DISTR_SET_CENTER_APPROX) ? " [guess]" : "");
  else
    fprintf(LOG,"%g [default]\n", unur_distr_cont_get_center(distr));

  fprintf(LOG,"%s:\tdomain = (%g, %g)",genid,DISTR.domain[0],DISTR.domain[1]);
  _unur_print_if_default(distr,UNUR_DISTR_SET_DOMAIN);

  fprintf(LOG,"\n%s:\tarea below p.d.f. = %g",genid,DISTR.area);
  _unur_print_if_default(distr,UNUR_DISTR_SET_PDFAREA);
  fprintf(LOG,"\n%s:\n",genid);

} /* end of _unur_distr_cont_debug() */

/*---------------------------------------------------------------------------*/
#endif    /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/
