/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: matr.c                                                            *
 *                                                                           *
 *   manipulate matrix distribution objects                                  *
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
#include <utils/matrix_source.h>
#include <distributions/unur_stddistr.h>
#include "distr_source.h"
#include "distr.h"
#include "matr.h"

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.matr

/*---------------------------------------------------------------------------*/

static void _unur_distr_matr_free( struct unur_distr *distr );

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** matrix distributions                                                    **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_matr_new( int n_rows, int n_cols )
     /*----------------------------------------------------------------------*/
     /* create a new (empty) distribution object                             */
     /* type: matrix distribution                                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   n_rows ... number of rows of matrix                                */
     /*   n_cols ... number of columns of matrix                             */
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
  if (n_rows < 1 || n_cols < 1) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"n_rows or n_cols < 1");
    return NULL;
  }

  /* get empty distribution object */
  distr = _unur_distr_generic_new();
  if (!distr) return NULL;

  /* set magic cookie */
  COOKIE_SET(distr,CK_DISTR_MATR);

  /* set type of distribution */
  distr->type = UNUR_DISTR_MATR;

  /* set id to generic distribution */
  distr->id = UNUR_DISTR_GENERIC;

  /* dimension of random vector */
  DISTR.n_rows = n_rows;
  DISTR.n_cols = n_cols;
  distr->dim = n_rows * n_cols;

  /* destructor */
  distr->destroy = _unur_distr_matr_free;

  /* clone */
  distr->clone = _unur_distr_matr_clone;

  /* set defaults */
  DISTR.init      = NULL;   /* pointer to special init routine (default: none) */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_matr_new() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
_unur_distr_matr_clone( const struct unur_distr *distr )
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
#define CLONE clone->data.matr

  struct unur_distr *clone;

  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, MATR, NULL );

  /* allocate memory */
  clone = _unur_xmalloc( sizeof(struct unur_distr) );
  
  /* copy distribution object into clone */
  memcpy( clone, distr, sizeof( struct unur_distr ) );

  /* copy data about distribution */

  /* copy user name for distribution */
  if (distr->name_str) {
    size_t len = strlen(distr->name_str) + 1;
    clone->name_str = _unur_xmalloc(len);
    memcpy( clone->name_str, distr->name_str, len );
    clone->name = clone->name_str;
  }

  return clone;

#undef CLONE
} /* end of _unur_distr_matr_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_distr_matr_free( struct unur_distr *distr )
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

  COOKIE_CHECK(distr,CK_DISTR_MATR,RETURN_VOID);

  /* user name for distribution */
  if (distr->name_str) free(distr->name_str);

  COOKIE_CLEAR(distr);
  free( distr );

} /* end of unur_distr_matr_free() */

/*---------------------------------------------------------------------------*/

int
unur_distr_matr_get_dim( const struct unur_distr *distr, int *n_rows, int *n_cols )
     /*----------------------------------------------------------------------*/
     /* get number of rows and columns of matrix                             */
     /* return total number of entries                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   n_rows ... to store number of rows                                 */
     /*   n_cols ... to store number of rows                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   total number of entries                                            */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, MATR, 0 );
  CHECK_NULL( n_rows, 0 );
  CHECK_NULL( n_cols, 0 );

  *n_rows = DISTR.n_rows;
  *n_cols = DISTR.n_cols;

  return distr->dim;

} /* end of unur_distr_matr_get_dim() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_distr_matr_debug( const struct unur_distr *distr, const char *genid )
     /*----------------------------------------------------------------------*/
     /* write info about distribution into LOG file                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   genid ... pointer to generator id                                  */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(distr,RETURN_VOID);
  COOKIE_CHECK(distr,CK_DISTR_MATR,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: distribution:\n",genid);
  fprintf(LOG,"%s:\ttype = matrix distribution\n",genid);
  fprintf(LOG,"%s:\tname = %s\n",genid,distr->name);

  fprintf(LOG,"%s:\tdimension = %d x %d   (= %d)\n",genid,DISTR.n_rows,DISTR.n_cols,distr->dim);

  fprintf(LOG,"%s:\n",genid);

} /* end of _unur_distr_matr_debug() */

/*---------------------------------------------------------------------------*/
#endif    /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/

