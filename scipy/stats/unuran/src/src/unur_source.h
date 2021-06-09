/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_source.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         To be included as first header file in all sources.               *
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
#ifndef UNUR_SOURCE_H_SEEN
#define UNUR_SOURCE_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* config file generated be autoconf                                         */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#else
#  error "config.h" required
#endif

/*---------------------------------------------------------------------------*/
/* compiler switches and defaults                                            */
#include <unuran_config.h>

/*---------------------------------------------------------------------------*/
/* define macros for GCC attributes                                          */

#ifdef __GNUC__
#  define ATTRIBUTE__FORMAT(a,b)   __attribute__ (( __format__ (printf, (a), (b)) ))
#  define ATTRIBUTE__UNUSED        __attribute__ ((unused))
#  define ATTRIBUTE__MALLOC        __attribute__ ((malloc))
#else
#  define ATTRIBUTE__FORMAT(a,b)
#  define ATTRIBUTE__UNUSED
#  define ATTRIBUTE__MALLOC
#endif

/*---------------------------------------------------------------------------*/
/* include standard header files                                             */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_LIMITS_H
#  include <limits.h>
#endif

/*---------------------------------------------------------------------------*/
/* globally used types                                                       */

#include <unur_typedefs.h>
#include <unur_struct.h>

/*---------------------------------------------------------------------------*/
/* Replacement functions                                                     */

#if !HAVE_DECL_LOG1P
#  include <specfunct/unur_specfunct_source.h> 
#endif

/*---------------------------------------------------------------------------*/
/* Utilities used by almost all sources                                      */

/* uniform random number generators */
#include <urng/urng_source.h>

/* magic cookies */
#include <unur_cookies.h>

/* debugging, warnings and error messages */
#include <utils/debug.h>
#include <utils/debug_source.h>
#include <utils/error.h>
#include <utils/error_source.h>
#include <utils/stream.h>
#include <utils/stream_source.h>
#include <utils/unur_errno.h>

/* floating point arithmetic */
#include <utils/unur_fp_source.h>
#include <utils/unur_fp_const_source.h>

/* mathematics */
#include <utils/umath.h>
#include <utils/umath_source.h>
#include <utils/unur_math_source.h>

/* vectors */
#include <utils/vector_source.h>

/* strings */
#include <utils/string_source.h>

/* allocate memory */
#include <utils/umalloc_source.h>

/* simple lists */
#include <utils/slist.h>

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_SOURCE_H_SEEN */
/*---------------------------------------------------------------------------*/
