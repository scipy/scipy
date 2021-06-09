/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      parser.c                                                     *
 *                                                                           *
 *   common routines for parsing strings.                                    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   REMARK:                                                                 *
 *   The parser always uses the default debugging flag.                      *
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

#include <ctype.h>
#include <unur_source.h>
#include "parser.h"
#include "parser_source.h"

/*---------------------------------------------------------------------------*/

char *
_unur_parser_prepare_string( const char *str )
     /*----------------------------------------------------------------------*/
     /* Prepare string for processing:                                       */
     /*   Make a working copy of the string.                                 */
     /*   Remove all white spaces and convert to lower case letters.         */
     /*   Single quotes (') are substituted with double quotes (")           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   str      ... pointer to string                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to working string.                                         */
     /*                                                                      */
     /* as a side effect, a new string is allocated.                         */
     /*----------------------------------------------------------------------*/
{
  char *tmp, *ptr;
  char *new;       /* pointer to working copy of string */
  size_t len;      /* length of string */

  /* length of string */
  len = strlen(str)+1;
  /* allocate memory for copy */
  new = _unur_xmalloc( len * sizeof(char) );
  /* copy memory */
  ptr = memcpy(new,str,len);

  /* copy characters but skip all white spaces */
  for (tmp = ptr; *tmp != '\0'; tmp++)
    if ( !isspace(*tmp) ) {
      *ptr = tolower(*tmp);
      /* substitute ' with " */
      if (*ptr == '\'') *ptr = '"';
      ptr++;
    }

  /* terminate string */
  *ptr = '\0';

  /* return pointer to working copy */
  return new;

} /* end of _unur_parser_prepare_string() */

/*---------------------------------------------------------------------------*/

