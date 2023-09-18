/* GStreamer
 * Copyright (C) 2010 Tim-Philipp Müller <tim centricular net>
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

#ifndef __GST_MATH_COMPAT_H__
#define __GST_MATH_COMPAT_H__

#ifndef __GI_SCANNER__

/* This header is not included automatically via gst/gst.h, you need to
 * include it explicitly if you need it. */

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES /* so MSVC defines M_PI etc. */
#endif
#include <math.h>

#include <glib.h>

G_BEGIN_DECLS

/* http://en.wikipedia.org/wiki/Math.h */

#define __GST_MATH_COMPAT_NEED_RINT
#define __GST_MATH_COMPAT_NEED_RINTF
#define __GST_MATH_COMPAT_NEED_ISNAN

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#undef __GST_MATH_COMPAT_NEED_RINT
#undef __GST_MATH_COMPAT_NEED_RINTF
#undef __GST_MATH_COMPAT_NEED_ISNAN
#endif

#if defined(_POSIX_VERSION) && _POSIX_VERSION >= 200112L
#undef __GST_MATH_COMPAT_NEED_RINT
#undef __GST_MATH_COMPAT_NEED_RINTF
#undef __GST_MATH_COMPAT_NEED_ISNAN
#endif

#ifndef M_PI
#define M_PI G_PI
#endif

#ifndef M_PI_2
#define M_PI_2 G_PI_2
#endif

#ifndef M_PI_4
#define M_PI_4 G_PI_4
#endif

static inline double
__gst_math_compat_rint (double x)
{
  return floor (x + 0.5);
}

static inline float
__gst_math_compat_rintf (float x)
{
  return floorf (x + 0.5);
}

static inline gboolean
__gst_math_compat_isnan (double x)
{
  return x != x;
}

#if defined (__GST_MATH_COMPAT_NEED_RINT) && !defined (rint)
#define rint(x) __gst_math_compat_rint(x)
#endif

#if defined (__GST_MATH_COMPAT_NEED_RINTF) && !defined (rintf)
#define rintf(x) __gst_math_compat_rintf(x)
#endif

#if defined (__GST_MATH_COMPAT_NEED_ISNAN) && !defined (isnan)
#define isnan(x) __gst_math_compat_isnan (x)
#endif

#ifndef NAN
#if G_BYTE_ORDER == G_BIG_ENDIAN
#define __GST_NAN_BYTES        { 0x7f, 0xc0, 0, 0 }
#elif G_BYTE_ORDER == G_LITTLE_ENDIAN
#define __GST_NAN_BYTES        { 0, 0, 0xc0, 0x7f }
#endif
static union {
  unsigned char __c[4];
  float __d;
} __gst_nan_union G_GNUC_UNUSED = {
  __GST_NAN_BYTES
};
#define NAN    (__gst_nan_union.__d)
#endif

G_END_DECLS

#endif /* !__GI_SCANNER__ */

#endif /* __GST_MATH_COMPAT_H__ */
