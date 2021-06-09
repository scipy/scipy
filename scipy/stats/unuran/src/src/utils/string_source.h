/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: string_source.h                                                   *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for strings                *
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
/* Function prototypes                                                       */

/* Make new string                                                           */
struct unur_string * _unur_string_new ( void );

/* Append to string                                                          */
/* Important: The generated string must not be longer than 1023 characters!  */
int _unur_string_append ( struct unur_string *string, const char *format, ... )
     ATTRIBUTE__FORMAT(2,3);

/* Append text to string                                                     */
int _unur_string_appendtext ( struct unur_string *string, const char *text );

/* Destroy string                                                            */
void _unur_string_free ( struct unur_string *string );

/* Clear string (set length of string to 0)                                  */
void _unur_string_clear ( struct unur_string *string );

/*---------------------------------------------------------------------------*/
