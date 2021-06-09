/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      condi.c                                                      *
 *                                                                           *
 *   manipulate univariate continuous full conditional distributions         *
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
#include <utils/matrix_source.h>
#include "distr.h"
#include "cont.h"
#include "condi.h"
#include "distr_source.h"

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "conditional";

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cvec    /* underlying (base) distribution          */
#define CONDI condi->data.cont    /* conditional distribution                */

#define iK           0            /* index of variable to be variated        */
#define K            (condi->data.cont.params[iK])

#define iPOSITION    0            /* condition (position of point)           */
#define iDIRECTION   1            /* direction of variation                  */
#define iXARG        2            /* working array for storing point x       */
#define iGRADF       3            /* working array for storing grad of PDF   */

#define POSITION     (condi->data.cont.param_vecs[iPOSITION])
#define DIRECTION    (condi->data.cont.param_vecs[iDIRECTION])
#define XARG         (condi->data.cont.param_vecs[iXARG])
#define GRADF        (condi->data.cont.param_vecs[iGRADF])

/*---------------------------------------------------------------------------*/

/* function prototypes                                                       */
static double _unur_pdf_condi( double x, const struct unur_distr *condi );
static double _unur_dpdf_condi( double x, const struct unur_distr *condi );

static double _unur_logpdf_condi( double x, const struct unur_distr *condi );
static double _unur_dlogpdf_condi( double x, const struct unur_distr *condi );

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** univariate continuous conditional distribution of multivaraiate         **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_condi_new( const struct unur_distr *distr, const double *pos, const double *dir, int k )
     /*----------------------------------------------------------------------*/
     /* Create an object for full conditional distribution .                 */
     /*   either for k-th variable and the other variables fixed to pos;     */
     /*   or     for pos + t*dir.                                            */
     /* The first one is chosen when dir == NULL.                            */
     /* `distr' must be a pointer to a multivariate continuous distribution. */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to multivariate continuous distribution.         */
     /*   pos   ... values for constant variables                            */
     /*   dir   ... direction (or NULL)                                      */
     /*   k     ... variable (must be in range 0 ... dim-1)                  */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to conditional distribution object                         */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_distr *condi;
  double *ar;

  /* check arguments */
  _unur_check_NULL( distr_name, distr, NULL );

  /* check distribution */
  if (distr->type != UNUR_DISTR_CVEC) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CVEC,NULL);

  /* check parameters pos and k */
  _unur_check_NULL( distr_name, pos, NULL );
  if ( dir==NULL && (k<0 || k>=distr->dim)) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,"k < 0 or k >= dim");
    return NULL;
  }

  /* get distribution object for generic continuous univariate distribution */
  condi = unur_distr_cont_new();
  if (!condi) return NULL;

  /* set id to distribution of conditional distribution */
  condi->id = UNUR_DISTR_CONDI;

  /* name of distribution */
  condi->name = distr_name;

  /* this is a derived distribution */
  /* clone base distribution ... */
  condi->base = _unur_distr_cvec_clone( distr );
  if (!condi->base) { _unur_distr_free(condi); return NULL; }

  /* set parameters for conditional distribution */
  CONDI.n_params = 1;                 /* one parameters: k */
  if ( unur_distr_condi_set_condition( condi, pos, dir, k ) != UNUR_SUCCESS ) {
    _unur_distr_free(condi); return NULL; 
  }

  /* we need two working arrays for computing PDF and dPDF */
  /* we abuse parameter vectors for this purpose */
  ar = _unur_xmalloc( distr->dim * sizeof(double) );
  memset( ar, 0, distr->dim * sizeof(double) );
  if ( (unur_distr_cont_set_pdfparams_vec( condi, iXARG, ar, distr->dim ) != UNUR_SUCCESS) ||
       (unur_distr_cont_set_pdfparams_vec( condi, iGRADF, ar, distr->dim ) != UNUR_SUCCESS) ) {
    _unur_distr_free(condi); free(ar); return NULL; 
  }
  free(ar);

  /* pointer to PDF, its derivative, and CDF */
  if (DISTR.pdf) {
    CONDI.pdf = _unur_pdf_condi;      /* pointer to PDF    */
    if (DISTR.dpdf)
      CONDI.dpdf = _unur_dpdf_condi;  /* derivative of PDF */
  }

  /* pointer to logPDF, its derivative */
  if (DISTR.logpdf) {
    CONDI.logpdf = _unur_logpdf_condi;      /* pointer to logPDF    */
    if (DISTR.dlogpdf)
      CONDI.dlogpdf = _unur_dlogpdf_condi;  /* derivative of logPDF */
  }

  /* return pointer to object */
  return condi;

} /* end of unur_distr_condi_new() */

/*---------------------------------------------------------------------------*/

int
unur_distr_condi_set_condition( struct unur_distr *condi, const double *pos, const double *dir, int k )
     /*----------------------------------------------------------------------*/
     /* Change values of fixed variables to pos and use k-th variable /      */
     /* direction dir.                                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   condi ... pointer to conditional distribution                      */
     /*   pos   ... values for constant variables                            */
     /*   dir   ... direction (or NULL)                                      */
     /*   k     ... variable (must be in range 0 ... dim-1)                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int dim;
  double *domain;

  /* check arguments */
  _unur_check_NULL( distr_name, condi, UNUR_ERR_NULL );
  _unur_check_distr_object( condi, CONT, UNUR_ERR_DISTR_INVALID );

  /* check distribution */
  if (condi->id != UNUR_DISTR_CONDI) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); 
    return UNUR_ERR_DISTR_INVALID; }
  COOKIE_CHECK(condi,CK_DISTR_CONT,UNUR_ERR_COOKIE);

  /* dimension of underlying distribution */
  dim = condi->base->dim;

  /* check parameters pos and k */
  _unur_check_NULL( condi->name, pos, UNUR_ERR_NULL );
  if ( dir== NULL && (k<0 || k>=dim)) {
    _unur_error(condi->name,UNUR_ERR_DISTR_INVALID,"k < 0 or k >= dim");
    return UNUR_ERR_DISTR_INVALID;
  }

  /* set parameters for conditional distribution */
  K = (double) k;
  if ( unur_distr_cont_set_pdfparams_vec( condi, iPOSITION, pos, dim ) != UNUR_SUCCESS ) {
    return UNUR_ERR_DISTR_INVALID;
  }

  if ( unur_distr_cont_set_pdfparams_vec( condi, iDIRECTION, dir, dim ) != UNUR_SUCCESS ) {
    return UNUR_ERR_DISTR_INVALID;
  }

  /* set domain of conditional distribution */
  if ( (domain = condi->base->data.cvec.domainrect) != NULL ) {
    if (dir == NULL) {
      CONDI.trunc[0] = CONDI.domain[0] = domain[2*k];
      CONDI.trunc[1] = CONDI.domain[1] = domain[2*k+1];
    }
    else {
      /* arbitrary direction: we ignore given domain */
      CONDI.trunc[0] = CONDI.domain[0] = -INFINITY;
      CONDI.trunc[1] = CONDI.domain[1] = INFINITY;
    }
  }

  /* changelog */
  condi->set &= ~UNUR_DISTR_SET_MODE; /* mode unknown */

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_condi_set_condition() */

/*---------------------------------------------------------------------------*/

int
unur_distr_condi_get_condition( struct unur_distr *condi, const double **pos, const double **dir, int *k )
     /*----------------------------------------------------------------------*/
     /* Read values of fixed variables, variable k / direction dir.          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   condi ... pointer to conditional distribution                      */
     /*   pos   ... pointer to array to store constant variables             */
     /*   dir   ... pointer to array to store direction                      */
     /*   k     ... pointer to store variable                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( distr_name, condi, UNUR_ERR_NULL );
  _unur_check_distr_object( condi, CONT, UNUR_ERR_DISTR_INVALID );

  /* check distribution */
  if (condi->id != UNUR_DISTR_CONDI) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); 
    return UNUR_ERR_DISTR_INVALID; }
  COOKIE_CHECK(condi,CK_DISTR_CONT,UNUR_ERR_COOKIE);

  /* position in vector */
  *k = (int) K;

  /* store pointer to values */
  *pos = POSITION;

  /* store pointer to direction */
  *dir = DIRECTION;

  /* return position */
  return UNUR_SUCCESS;
} /* end of unur_distr_condi_get_condition() */

/*---------------------------------------------------------------------------*/

const struct unur_distr *
unur_distr_condi_get_distribution( const struct unur_distr *condi )
     /*----------------------------------------------------------------------*/
     /* get pointer to distribution object for underlying distribution       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   condi ... pointer to conditional distribution object               */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to underlying distribution                                 */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( distr_name, condi, NULL );
  _unur_check_distr_object( condi, CONT, NULL );

  /* check distribution */
  if (condi->id != UNUR_DISTR_CONDI) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_INVALID,"");
    return NULL; 
  }

  return condi->base;
} /* end of unur_distr_condi_get_distribution() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** PDF, its derivative and CDF of conditional distribution                 **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

double
_unur_pdf_condi( double x, const struct unur_distr *condi )
     /* 
	PDF(x) = mvpdf(p(x))

	mvpdf(.) ... PDF of underlying multivariate distribution
	p(x)     ... vector with x in k-th position and pos in all others
     */
{ 
  int dim = condi->base->dim;   /* dimension of underlying distribution */
  int k = (int) K;              /* position in vector */
  int i;

  /* set point for multivariate PDF */
  if (DIRECTION==NULL) {  /* use k-th variable */
    memcpy(XARG, POSITION, dim * sizeof(double) );
    XARG[k] = x;  
  }
  else {   /* use direction vector */
    memcpy(XARG, POSITION, dim * sizeof(double) );
    for (i=0; i<dim; i++)
      XARG[i] += x*DIRECTION[i];
  }

  return _unur_cvec_PDF(XARG,condi->base);

} /* end of _unur_pdf_condi() */

/*---------------------------------------------------------------------------*/

double
_unur_logpdf_condi( double x, const struct unur_distr *condi )
     /* 
	logPDF(x) = logmvpdf(p(x))

	mvpdf(.) ... logPDF of underlying multivariate distribution
	p(x)     ... vector with x in k-th position and pos in all others
     */
{ 
  int dim = condi->base->dim;   /* dimension of underlying distribution */
  int k = (int) K;              /* position in vector */
  int i;

  /* set point for multivariate logPDF */
  if (DIRECTION==NULL) {  /* use k-th variable */
    memcpy(XARG, POSITION, dim * sizeof(double) );
    XARG[k] = x;  
  }
  else {   /* use direction vector */
    memcpy(XARG, POSITION, dim * sizeof(double) );
    for (i=0; i<dim; i++)
      XARG[i] += x*DIRECTION[i];
  }

  return _unur_cvec_logPDF(XARG,condi->base);

} /* end of _unur_logpdf_condi() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_condi( double x, const struct unur_distr *condi )
     /* 
	(dPDF/dxk)(x) = (grad mvpdf(p(x)))[k]

	mvpdf(.) ... PDF of underlying multivariate distribution
	p(x)     ... vector with x in k-th position and pos in all others
     */
{
  int dim = condi->base->dim;    /* dimension of underlying distribution */
  int k = (int) K;               /* position in vector */
  int i;
  double df;

  if (DIRECTION==NULL) {  /* use k-th variable */
    /* set point for multivariate PDF */
    memcpy(XARG, POSITION, dim * sizeof(double) );
    XARG[k] = x;  
    if (condi->base->data.cvec.pdpdf) {
      /* we have a pointer to the partial derivative */
      df = _unur_cvec_pdPDF(XARG, k, condi->base);
    }
    else {
      /* we do not have partial derivatives --> have to take coordinate from gradient */
      /* compute gradient */
      _unur_cvec_dPDF(GRADF, XARG, condi->base);
      /* return k-th component */
      df = GRADF[k];
    }
  }

  else {   /* use direction vector */
    memcpy(XARG, POSITION, dim * sizeof(double) );
    for (i=0; i<dim; i++)
      XARG[i] += x*DIRECTION[i];
    _unur_cvec_dPDF(GRADF, XARG, condi->base);
    for (df=0.,i=0; i<dim; i++)
      df += GRADF[i]*DIRECTION[i];
  }  

  return df;

} /* end of _unur_dpdf_condi() */

/*---------------------------------------------------------------------------*/

double
_unur_dlogpdf_condi( double x, const struct unur_distr *condi )
     /* 
	(dlogPDF/dxk)(x) = (grad logmvpdf(p(x)))[k]

	logmvpdf(.) ... PDF of underlying multivariate distribution
	p(x)        ... vector with x in k-th position and pos in all others
     */
{
  int dim = condi->base->dim;    /* dimension of underlying distribution */
  int k = (int) K;               /* position in vector */
  int i;
  double df;

  if (DIRECTION==NULL) {  /* use k-th variable */
    /* set point for multivariate logPDF */
    memcpy(XARG, POSITION, dim * sizeof(double) );
    XARG[k] = x;  
    if (condi->base->data.cvec.pdlogpdf) {
      /* we have a pointer to the partial derivative */
      df = _unur_cvec_pdlogPDF(XARG, k, condi->base);
    }
    else {
      /* we do not have partial derivatives --> have to take coordinate from gradient */
      /* compute gradient */
      _unur_cvec_dlogPDF(GRADF, XARG, condi->base);
      /* return k-th component */
      df = GRADF[k];
    }
  }

  else {   /* use direction vector */
    memcpy(XARG, POSITION, dim * sizeof(double) );
    for (i=0; i<dim; i++)
      XARG[i] += x*DIRECTION[i];
    _unur_cvec_dlogPDF(GRADF, XARG, condi->base);
    for (df=0.,i=0; i<dim; i++)
      df += GRADF[i]*DIRECTION[i];
  }  

  return df;

} /* end of _unur_dlogpdf_condi() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** debug                                                                   **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_distr_condi_debug( const struct unur_distr *condi, const char *genid )
     /*----------------------------------------------------------------------*/
     /* write info about distribution into LOG file                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   condi ... pointer to conditional distribution                      */
     /*   genid ... pointer to generator id                                  */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(condi,RETURN_VOID);
  COOKIE_CHECK(condi,CK_DISTR_CONT,RETURN_VOID);
  CHECK_NULL(condi->base,RETURN_VOID);

  LOG = unur_get_stream();

  /* print data about conditional distribution */
  fprintf(LOG,"%s: distribution:\n",genid);
  fprintf(LOG,"%s:\ttype = full conditional distribution of continuous multivariate distribution\n",genid);
  fprintf(LOG,"%s:\tname = %s\n",genid,condi->name);
  fprintf(LOG,"%s:\n",genid);

  /* number of the variable for which the conditional distribution 
     is to be shown */
  fprintf(LOG,"%s:\tcondition (at time of creation):\n",genid);
  if (DIRECTION==NULL) {
    fprintf(LOG,"%s:\tvariable = %d\n",genid,(int)(K));
    _unur_matrix_print_vector( condi->base->dim, POSITION, "\tpoint =", LOG, genid, "\t   ");
  }
  else {
    _unur_matrix_print_vector( condi->base->dim, DIRECTION, "\tdirection =", LOG, genid, "\t   ");
    _unur_matrix_print_vector( condi->base->dim, POSITION, "\tpoint =", LOG, genid, "\t   ");
  }
  fprintf(LOG,"%s:\n",genid);

  /* domain */
  fprintf(LOG,"%s:\tdomain = (%g, %g)\n",genid,CONDI.domain[0],CONDI.domain[1]);
  fprintf(LOG,"%s:\n",genid);

  /* print data about underlying distribution */
  fprintf(LOG,"%s: Underlying distribution:\n",genid);
  _unur_distr_cvec_debug(condi->base, genid);

} /* end of _unur_distr_condi_debug() */
/*---------------------------------------------------------------------------*/
#endif    /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
