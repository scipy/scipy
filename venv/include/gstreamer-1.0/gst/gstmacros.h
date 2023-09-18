/* GStreamer
 * Copyright (C) <1999> Erik Walthinsen <omega@cse.ogi.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */
#ifndef __GST_MACROS_H__
#define __GST_MACROS_H__

#include <glib.h>

G_BEGIN_DECLS

#if __GNUC__ > 2 || (__GNUC__ == 2 && __GNUC_MINOR__ > 4)
# define GST_GNUC_CONSTRUCTOR \
  __attribute__ ((constructor))
#else  /* !__GNUC__ */
# define GST_GNUC_CONSTRUCTOR
#endif /* !__GNUC__ */

#if defined (__GNUC__) && !defined (GST_IMPLEMENT_INLINES)
# define GST_INLINE_FUNC extern __inline__
# define GST_CAN_INLINE 1
#elif defined(_MSC_VER)
# define GST_INLINE_FUNC extern __inline
# define GST_CAN_INLINE 1
#else
# define GST_INLINE_FUNC extern
# undef GST_CAN_INLINE
#endif

/* MSVC defines 'restrict' as a keyword and not a define */
#if (!defined(__STDC_VERSION__) || __STDC_VERSION__ < 199901L) && !defined(restrict)
#  if defined(__GNUC__) && __GNUC__ >= 4
#    define restrict __restrict__
#  elif defined(_MSC_VER)
#    define restrict __restrict
#  else
#    define restrict
#  endif
#endif

G_END_DECLS

#endif /* __GST_MACROS_H__ */


