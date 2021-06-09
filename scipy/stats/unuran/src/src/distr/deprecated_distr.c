/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: deprecated_distr.c                                                *
 *                                                                           *
 *   Deprecated routines                                                     *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   THESE ROUTINES SHOULD NOT BE USED ANY MORE!                             *
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
#include <stdarg.h>
#include <distr/distr_source.h>
#include <distr/distr.h>
#include "deprecated_distr.h"

/*---------------------------------------------------------------------------*/
#ifdef USE_DEPRECATED_CODE
/*---------------------------------------------------------------------------*/

static void _unur_distr_cvec_marginals_free ( struct unur_distr **marginals, int dim );

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#define DISTR distr->data.cvec
/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_set_stdmarginals ( struct unur_distr *distr, struct unur_distr *stdmarginal)
     /*----------------------------------------------------------------------*/
     /* Copy standardized marginal distribution into distribution object.    */
     /* Only one local copy is made.                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr       ... pointer to distribution object                     */
     /*   stdmarginal ... pointer to standardized marginal distrib. object   */
     /*                                                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_distr *clone;
  int i;

  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( distr->name, stdmarginal, UNUR_ERR_NULL );
  _unur_check_distr_object( stdmarginal, CONT, UNUR_ERR_DISTR_INVALID );

  /* first we have to check whether there is already a list of marginal distributions */
  if (DISTR.stdmarginals)
    _unur_distr_cvec_marginals_free(DISTR.stdmarginals, distr->dim);

  /* make copy of standardized marginal distribution object */
  clone = _unur_distr_clone( stdmarginal );

  /* allocate memory for array */
  DISTR.stdmarginals = _unur_xmalloc (distr->dim * sizeof(struct unur_distr *));

  /* copy pointer */
  for (i=0; i<distr->dim; i++)
    DISTR.stdmarginals[i] = clone;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_STDMARGINAL;

  return UNUR_SUCCESS;
} /* end of unur_distr_cvec_set_stdmarginals() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_set_stdmarginal_array ( struct unur_distr *distr, struct unur_distr **stdmarginals)
     /*----------------------------------------------------------------------*/
     /* Copy standardized marginal distributions into distribution object.   */
     /* For each dimension a new copy is made even if the pointer in         */
     /* the array stdmarginals coincide.                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr        ... pointer to distribution object                    */
     /*   stdmarginals ... pointer to array of std. marginal distr. objects  */
     /*                                                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( distr->name, stdmarginals, UNUR_ERR_NULL );

  for (i=0; i<distr->dim; i++) {
    _unur_check_NULL( distr->name, *(stdmarginals+i), UNUR_ERR_NULL );
    _unur_check_distr_object( *(stdmarginals+i), CONT, UNUR_ERR_DISTR_INVALID );
  }
    
  /* first we have to check whether there is already a list of marginal distributions */
  if (DISTR.stdmarginals)
    _unur_distr_cvec_marginals_free(DISTR.stdmarginals, distr->dim);

  /* allocate memory for array */
  DISTR.stdmarginals = _unur_xmalloc (distr->dim * sizeof(struct unur_distr *));

  /* make copy of standardized marginal distribution objects */
  for (i=0; i<distr->dim; i++) 
    DISTR.stdmarginals[i] = _unur_distr_clone( *(stdmarginals+i) );

  /* changelog */
  distr->set |= UNUR_DISTR_SET_STDMARGINAL;

  return UNUR_SUCCESS;
} /* end of unur_distr_cvec_set_stdmarginal_array() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_set_stdmarginal_list ( struct unur_distr *distr, ... )
     /*----------------------------------------------------------------------*/
     /* Copy standardized marginal distributions into distribution object.   */
     /* For each dimenision there must be a pointer to a discribution object.*/
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   ...    ... pointer to array of std. marginal distibution objects   */
     /*                                                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* IMPORTANT:                                                           */
     /*   All marginal distribution objects are destroyed after they have    */
     /*   been copied into the distribution object.                          */
     /*----------------------------------------------------------------------*/
{
  int i;
  int failed = FALSE;
  struct unur_distr *stdmarginal;
  struct unur_distr **stdmarginal_list;
  va_list vargs;

  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );

  /* allocate memory for array */
  stdmarginal_list = _unur_xmalloc (distr->dim * sizeof(struct unur_distr *));
  for (i=0; i<distr->dim; i++) stdmarginal_list[i] = NULL;

  /* make copy of standardized marginal distribution objects */
  va_start(vargs, distr);
  for (i=0; i<distr->dim; i++) {
    stdmarginal = (struct unur_distr *) va_arg(vargs, struct unur_distr *);
    if (stdmarginal) {
      stdmarginal_list[i] = _unur_distr_clone( stdmarginal );
      _unur_distr_free(stdmarginal);
    }
    else {
      failed = TRUE;
    }
  }
  va_end(vargs);

  if (failed) {
    /* some of the pointers are NULL pointers */
    _unur_distr_cvec_marginals_free(stdmarginal_list, distr->dim);
    _unur_error(distr->name ,UNUR_ERR_DISTR_SET,"stdmarginals == NULL");
    return UNUR_ERR_DISTR_SET;
  }

  /* copy list of marginal distributions. However, first we have to check  */
  /* whether there is already a list of marginal distributions.            */
  if (DISTR.stdmarginals)
    _unur_distr_cvec_marginals_free(DISTR.stdmarginals, distr->dim);

  DISTR.stdmarginals = stdmarginal_list;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_STDMARGINAL;

  return UNUR_SUCCESS;
} /* end of unur_distr_cvec_set_stdmarginal_list() */

/*---------------------------------------------------------------------------*/

const struct unur_distr *
unur_distr_cvec_get_stdmarginal( const struct unur_distr *distr, int n )
     /*----------------------------------------------------------------------*/
     /* Get pointer to standardized marginal distribution object.            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   n     ... position of standardized marginal distribution object    */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );

  if (n<=0 || n > distr->dim) {
    _unur_error(distr->name,UNUR_ERR_DISTR_GET,"n not in 1 .. dim");
    return NULL;
  }

  /* mean vector known ? */
  if ( !(distr->set & UNUR_DISTR_SET_STDMARGINAL) ) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_GET,"std marginals");
    return NULL;
  }

  _unur_check_NULL( distr->name, DISTR.stdmarginals, NULL );

  /* return standarized marginal distribution object */
  return (DISTR.stdmarginals[n-1]);
} /* end of unur_distr_cvec_get_stdmarginal() */

/*---------------------------------------------------------------------------*/

void
_unur_distr_cvec_marginals_free ( struct unur_distr **marginals, int dim )
     /*----------------------------------------------------------------------*/
     /* free list of marginal distribution objects                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   marginals ... pointer to list of marginal distribution objects     */
     /*   dim       ... number of marginal distributions                     */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of list of marginal distribution objects          */
     /*----------------------------------------------------------------------*/
{
  int i;

  if (_unur_distr_cvec_marginals_are_equal(marginals,dim)) {
    _unur_distr_free(marginals[0]);
  }

  else {
    for (i=0; i<dim; i++) 
      if (marginals[i]) _unur_distr_free(marginals[i]);
  }

  free (marginals);
} /* end of _unur_distr_cvec_marginals_free() */

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
#endif   /* USE_DEPRECATED_CODE */
/*---------------------------------------------------------------------------*/
