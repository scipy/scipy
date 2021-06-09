/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: string.c                                                          *
 *                                                                           *
 *   routines for handling strings                                           *
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

#include <time.h>
#include <stdarg.h>
#include <ctype.h>

#define MEMBLOCKSIZE     128    /* block size for allocating memory */
#define MAXSTRINGSIZE   1024    /* maximum size of printed string   */

/*---------------------------------------------------------------------------*/

struct unur_string *
_unur_string_new ( void )
     /*----------------------------------------------------------------------*/
     /* Make new string                                                      */
     /*                                                                      */
     /* parameters: none                                                     */
     /*----------------------------------------------------------------------*/
{
  struct unur_string *string;
  
  string = _unur_xmalloc(sizeof(struct unur_string));
  string->text = NULL;
  string->length = 0;
  string->allocated = 0;

  return string;
} /* end of _unur_string_new() */

/*---------------------------------------------------------------------------*/

int
_unur_string_append ( struct unur_string *string, const char *format, ... )
     /*----------------------------------------------------------------------*/
     /* Append to string                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   string ... structure for string                                    */
     /*   format ... format for sprintf()                                    */
     /*   ...    ... (optional) arguments to be be printed                   */
     /*                                                                      */
     /* Important!                                                           */
     /*   The generated string must not be longer than 1023 characters!      */
     /*----------------------------------------------------------------------*/
{
  size_t len;
  va_list ap;

  /* optional list of arguments */
  va_start(ap, format);

  /* Resize the allocated memory if necessary */
  while (string->length + MAXSTRINGSIZE + 1 > string->allocated) {
    string->allocated += MEMBLOCKSIZE;
    string->text = _unur_xrealloc( string->text, (size_t)string->allocated );
  }

  /* print into string */
#if HAVE_DECL_VSNPRINTF
  /* this is a GNU extension */
  len = vsnprintf (string->text+string->length, (size_t)MAXSTRINGSIZE, format, ap);
#else
  /** TODO: this is dangerous, since we have to take care, that
      the generated string text is not longer than MAXSTRINGSIZE-1.  **/
  len = vsprintf (string->text+string->length, format, ap);
  if (len >= MAXSTRINGSIZE) {
    _unur_error("UTIL",UNUR_ERR_SHOULD_NOT_HAPPEN,"string too long");
    exit (-1);   /* fatal error */
  }
#endif

  /* update length of string */
  string->length += len;

  /* close optional list of arguments */
  va_end(ap);

  return UNUR_SUCCESS;
} /* end of _unur_string_append() */

/*---------------------------------------------------------------------------*/

int 
_unur_string_appendtext ( struct unur_string *string, const char *text )
     /*----------------------------------------------------------------------*/
     /* Append text to string                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   string ... structure for stringing string                          */
     /*   text   ... char array with text to be appended                     */
     /*----------------------------------------------------------------------*/
{
  int len;

  /* length of text string */
  len = strlen(text);
  
  /* Resize the allocated memory if necessary */
  while (string->length + len + 1 > string->allocated) {
    string->allocated += MEMBLOCKSIZE;
    string->text = _unur_xrealloc( string->text, (size_t)string->allocated );
  }

  /* copy text into string */
  strncpy( string->text+string->length, text, len+1 );

  /* update length of string */
  string->length += len;

  return UNUR_SUCCESS;
} /* end of _unur_string_appendtext() */

/*---------------------------------------------------------------------------*/

void
_unur_string_free ( struct unur_string *string )
     /*----------------------------------------------------------------------*/
     /* Destroy string (free all memory)                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   string ... structure for stringing string                          */
     /*----------------------------------------------------------------------*/
{
  if (string) {
    if (string->text)  free (string->text);
    free (string);
    string = NULL;
  }
} /* end of _unur_string_free() */

/*---------------------------------------------------------------------------*/

void
_unur_string_clear ( struct unur_string *string )
     /*----------------------------------------------------------------------*/
     /* Clear string (set length of string to 0)                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   string ... structure for stringing string                          */
     /*----------------------------------------------------------------------*/
{
  if (string) {
    string->length = 0;
    *(string->text) = '\0';
  }
} /* end of _unur_string_clear() */

/*---------------------------------------------------------------------------*/

