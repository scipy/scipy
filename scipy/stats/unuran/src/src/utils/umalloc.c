/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: umalloc.c                                                         *
 *                                                                           *
 *   allocate memory                                                         *
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

void*
_unur_xmalloc(size_t size)
     /*----------------------------------------------------------------------*/
     /* allocate memory                                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   size ... size of allocated block                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   abort program                                                      */
     /*----------------------------------------------------------------------*/
{
  register void *ptr;

  /* allocate memory */
  ptr = malloc( size );

  /* successful ? */
  if (ptr == NULL) {
    _unur_error(NULL,UNUR_ERR_MALLOC,"");
    exit (EXIT_FAILURE);
  }

  return ptr;

} /* end of _unur_xmalloc() */

/*---------------------------------------------------------------------------*/

void*
_unur_xrealloc(void *ptr, size_t size)
     /*----------------------------------------------------------------------*/
     /* reallocate memory                                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   ptr  ... address of memory block previously allocated by malloc.   */
     /*   size ... size of reallocated block                                 */
     /*                                                                      */
     /* error:                                                               */
     /*   abort program                                                      */
     /*----------------------------------------------------------------------*/
{
  register void *new_ptr;

  /* reallocate memory */
  
  new_ptr = realloc( ptr, size );

  /* successful ? */
  if (new_ptr == NULL) {
    _unur_error(NULL,UNUR_ERR_MALLOC,"");
    exit (EXIT_FAILURE);
  }

  return new_ptr;

} /* end of _unur_xrealloc() */

/*---------------------------------------------------------------------------*/
