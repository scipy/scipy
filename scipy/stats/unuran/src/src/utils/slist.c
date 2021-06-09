/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: slist.c                                                           *
 *                                                                           *
 *   Handling simple lists                                                   *
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

/*---------------------------------------------------------------------------*/

struct unur_slist *
_unur_slist_new( void )
     /*----------------------------------------------------------------------*/
     /* Make new simple list.                                                */
     /* Append pointer to element to simple list.                            */
     /*                                                                      */
     /* parameters: none                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to new empty simple list                                   */
     /*----------------------------------------------------------------------*/
{
  struct unur_slist *slist;

  /* allocate structure */
  slist = _unur_xmalloc(sizeof(struct unur_slist));
  COOKIE_SET(slist,CK_SLIST);

  /* initialize */

  slist->ptr   = NULL;
  slist->n_ptr = 0;

  return slist;
} /* end of _unur_slist_new() */

/*---------------------------------------------------------------------------*/

int
_unur_slist_length( const struct unur_slist *slist )
     /*----------------------------------------------------------------------*/
     /* Get length if list (number of list entries).                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   slist   ... pointer to simple list                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   number of elements                                                 */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(slist,0);
  COOKIE_CHECK(slist,CK_SLIST,0);

  if (slist->ptr==NULL)
    return 0;

  return (slist->n_ptr);

} /* end of _unur_slist_length() */

/*---------------------------------------------------------------------------*/

void *
_unur_slist_get( const struct unur_slist *slist, int n )
     /*----------------------------------------------------------------------*/
     /* Get pointer to n-th element.                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   slist   ... pointer to simple list                                 */
     /*   n       ... index element                                          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to element                                                 */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(slist,NULL);
  COOKIE_CHECK(slist,CK_SLIST,NULL);

  if (slist->ptr==NULL || n >= slist->n_ptr || n < 0) {
    _unur_warning("list",UNUR_ERR_GENERIC,"element does not exist");
    return NULL;
  }
  
  return (slist->ptr[n]);

} /* end of _unur_slist_get() */

/*---------------------------------------------------------------------------*/

int
_unur_slist_append( struct unur_slist *slist, void *element )
     /*----------------------------------------------------------------------*/
     /* Append pointer to element to simple list.                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   slist   ... pointer to simple list                                 */
     /*   element ... pointer to element to be appended                      */
     /*                                                                      */
     /* return:                                                              */
     /*   success ... UNUR_SUCCESS                                           */
     /*   error   ... error code                                             */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(slist,UNUR_ERR_NULL);
  COOKIE_CHECK(slist,CK_SLIST,UNUR_ERR_COOKIE);

  /* allocate memory for the list of blocks */
  slist->ptr = _unur_xrealloc(slist->ptr,(slist->n_ptr+1)*sizeof(void *));

  /* store allocated element */
  slist->ptr[slist->n_ptr] = element;

  /* update number of allocated elements */
  ++(slist->n_ptr);

  return UNUR_SUCCESS;

} /* end of _unur_slist_append() */

/*---------------------------------------------------------------------------*/

void *
_unur_slist_replace( struct unur_slist *slist, int n, void *element )
     /*----------------------------------------------------------------------*/
     /* Replace (existing) pointer to n-th element by 'element'.             */
     /*                                                                      */
     /* If n-th element does not exist, no new element is inserted!          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   slist   ... pointer to simple list                                 */
     /*   n       ... index element                                          */
     /*   element ... pointer to new element to be inserted                  */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to old element                                             */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  void *old_element; 

  /* check arguments */
  CHECK_NULL(slist,NULL);
  COOKIE_CHECK(slist,CK_SLIST,NULL);

  if (slist->ptr==NULL || n >= slist->n_ptr || n < 0) {
    _unur_warning("list",UNUR_ERR_GENERIC,"element does not exist");
    return NULL;
  }

  old_element = slist->ptr[n];
  slist->ptr[n] = element;
  return old_element;

} /* end of _unur_slist_replace() */

/*---------------------------------------------------------------------------*/

void 
_unur_slist_free( struct unur_slist *slist )
     /*----------------------------------------------------------------------*/
     /* Free all elements and list in simple list.                           */
     /* compute "arctan mean" of two numbers.                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   slist  ... pointer to simple list                                  */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* check arguments */
  if (slist == NULL) return;  /* nothing to do */
  COOKIE_CHECK(slist,CK_SLIST,RETURN_VOID);

  if ( slist->ptr != NULL ) {
    /* free memory blocks */
    for (i=0; i < slist->n_ptr; i++)
      if (slist->ptr[i]) free(slist->ptr[i]); 
    free(slist->ptr);
    slist->ptr = NULL;
  }

  /* free list */
  free (slist);

} /* end of _unur_slist_free() */

/*---------------------------------------------------------------------------*/

