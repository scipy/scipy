/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:  cvemp.c                                                          *
 *                                                                           *
 *   manipulate empirical multivariate continuous distribution objects       *
 *   (i.e. samples)                                                          *
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
#include "distr_source.h"
#include "distr.h"
#include "cvemp.h"

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cvemp

/*---------------------------------------------------------------------------*/

static void _unur_distr_cvemp_free( struct unur_distr *distr );

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** empirical univariate continuous distributions                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_cvemp_new( int dim )
     /*----------------------------------------------------------------------*/
     /* create a new (empty) distribution object                             */
     /* type: empirical multivariate continuous                              */
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

  /* check dimension for new parameter for distribution */
  if (dim < 2) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"dimension < 2");
    return NULL;
  }

  /* get empty distribution object */
  distr = _unur_distr_generic_new();
  if (!distr) return NULL;

  /* set magic cookie */
  COOKIE_SET(distr,CK_DISTR_CVEMP);

  /* set type of distribution */
  distr->type = UNUR_DISTR_CVEMP;

  /* set id to generic distribution */
  distr->id = UNUR_DISTR_GENERIC;

  /* dimension of random vector */
  distr->dim = dim;   /* multivariant */

  /* name for distribution */
  distr->name = "(empirical)";
  distr->name_str = NULL;

  /* destructor */
  distr->destroy = _unur_distr_cvemp_free;

  /* clone */
  distr->clone = _unur_distr_cvemp_clone;

  /* set defaults */

  /* observed sample */
  DISTR.sample    = NULL;    /* sample      */
  DISTR.n_sample  = 0;       /* sample size */
  
  /* return pointer to object */
  return distr;

} /* end of unur_distr_cvemp_new() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
_unur_distr_cvemp_clone( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* copy (clone) distribution object                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to source distribution object                    */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of distribution object                            */
     /*----------------------------------------------------------------------*/
{
#define CLONE clone->data.cvemp

  struct unur_distr *clone;

  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEMP, NULL );

  /* allocate memory */
  clone = _unur_xmalloc( sizeof(struct unur_distr) );
  
  /* copy distribution object into clone */
  memcpy( clone, distr, sizeof( struct unur_distr ) );

  /* copy data about sample into generator object (when there is one) */
  if (DISTR.sample) {
    CLONE.sample = _unur_xmalloc( DISTR.n_sample * distr->dim * sizeof(double) );
    memcpy( CLONE.sample, DISTR.sample, DISTR.n_sample * distr->dim * sizeof(double) );
  }

  /* copy user name for distribution */
  if (distr->name_str) {
    size_t len = strlen(distr->name_str) + 1;
    clone->name_str = _unur_xmalloc(len);
    memcpy( clone->name_str, distr->name_str, len );
    clone->name = clone->name_str;
  }

  return clone;

#undef CLONE
} /* end of _unur_distr_cvemp_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_distr_cvemp_free( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* free distribution object                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  if( distr == NULL ) /* nothing to do */
    return;

  COOKIE_CHECK(distr,CK_DISTR_CVEMP,RETURN_VOID);

  if (DISTR.sample) free( DISTR.sample );

  /* user name for distribution */
  if (distr->name_str) free(distr->name_str);

  COOKIE_CLEAR(distr);
  free( distr );

} /* end of unur_distr_cvemp_free() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvemp_set_data( struct unur_distr *distr, const double *sample, int n_sample )
     /*----------------------------------------------------------------------*/
     /* set observed sample for distribution                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   sample   ... pointer to array of observations                      */
     /*   n_sample ... number of observations (sample size)                  */
     /*                                                                      */
     /* sample must be an array of size dim x n_sample                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEMP, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( distr->name, sample, UNUR_ERR_NULL );

  /* check new parameter for generator */
  if (n_sample <= 0) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"sample size");
    return UNUR_ERR_DISTR_SET;
  }

  /* allocate memory for sample */
  DISTR.sample = _unur_xmalloc( n_sample * distr->dim * sizeof(double) );
  if (!DISTR.sample) return UNUR_ERR_MALLOC;

  /* copy observed sample */
  memcpy( DISTR.sample, sample, n_sample * distr->dim * sizeof(double) );
  DISTR.n_sample = n_sample;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_cvemp_set_sample() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvemp_read_data( struct unur_distr *distr, const char *filename )
     /*----------------------------------------------------------------------*/
     /* Read data from file `filename'.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   filename ... name of data file                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEMP, UNUR_ERR_DISTR_INVALID );

  /* read data from file */
  DISTR.n_sample = _unur_read_data( filename, distr->dim, &(DISTR.sample) );

  /* o.k. ? */
  return (DISTR.n_sample > 0) ? UNUR_SUCCESS : UNUR_ERR_DISTR_DATA;

} /* end of unur_distr_cvemp_read_data() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_cvemp_get_data( const struct unur_distr *distr, const double **sample )
     /*----------------------------------------------------------------------*/
     /* get number of observations and set pointer to array observations     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   sample   ... pointer to observed sampletor                         */
     /*                                                                      */
     /* return:                                                              */
     /*   sample size                                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CVEMP, 0 );

  *sample = (DISTR.sample) ? DISTR.sample : NULL;
  return DISTR.n_sample;

} /* end of unur_distr_cvemp_get_sample() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_distr_cvemp_debug( const struct unur_distr *distr, const char *genid, unsigned printvector )
     /*----------------------------------------------------------------------*/
     /* write info about distribution into LOG file                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   genid ... pointer to generator id                                  */
     /*   printvector ... print observed sample if not 0                     */
     /*----------------------------------------------------------------------*/
{
#define idx(k,l)  (k * distr->dim + l)

  FILE *LOG;
  int i,j;

  /* check arguments */
  CHECK_NULL(distr,RETURN_VOID);
  COOKIE_CHECK(distr,CK_DISTR_CVEMP,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: distribution:\n",genid);
  fprintf(LOG,"%s:\ttype = continuous multivariate distribution (ie. a sample)\n",genid);
  fprintf(LOG,"%s:\tname = %s\n",genid,distr->name);

  fprintf(LOG,"%s:\tdimension = %d\n",genid,distr->dim);

  if (DISTR.n_sample>0) {
    /* observed samples */
    fprintf(LOG,"%s:\tsample size = %d\n",genid,DISTR.n_sample);
    if (printvector) {
      for (i=0; i<DISTR.n_sample; i++) {
	fprintf(LOG,"%s:\t( %.5f",genid,DISTR.sample[idx(i,0)]);
	for (j=1; j<distr->dim; j++) 
	  fprintf(LOG,", %.5f",DISTR.sample[idx(i,j)]);
	fprintf(LOG,")\n");
      }
    }
  }
  fprintf(LOG,"%s:\n",genid);

#undef idx
} /* end of _unur_distr_cvemp_debug() */

/*---------------------------------------------------------------------------*/
#endif    /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
