 /*
 * GStreamer
 * Copyright (C) 2012 Matthew Waters <ystreet00@gmail.com>
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

#ifndef __GST_GL_COMPAT_H__
#define __GST_GL_COMPAT_H__

#include <gst/gl/gstglconfig.h>

/* undefined typedefs */
#if !GST_GL_HAVE_GLEGLIMAGEOES
typedef gpointer GLeglImageOES;
#endif
#if !GST_GL_HAVE_GLCHAR
typedef gchar GLchar;
#endif
#if !GST_GL_HAVE_GLSIZEIPTR
typedef ptrdiff_t GLsizeiptr;
#endif
#if !GST_GL_HAVE_GLINTPTR
typedef ptrdiff_t GLintptr;
#endif
#if !GST_GL_HAVE_GLSYNC
typedef gpointer GLsync;
#endif
#if !GST_GL_HAVE_GLUINT64
typedef guint64 GLuint64;
#endif
#if !GST_GL_HAVE_GLINT64
typedef gint64 GLint64;
#endif

#if !defined(GST_GL_DEBUG_PROC)
#if defined(GLDEBUGPROC)
#define GST_GL_DEBUG_PROC GLDEBUGPROC
#elif defined(GLDEBUGPROCARB)
#define GST_GL_DEBUG_PROC GLDEBUGPROCARB
#elif defined(GLDEBUGPROCKHR)
#define GST_GL_DEBUG_PROC GLDEBUGPROCKHR
#else
typedef void (GSTGLAPI *GST_GL_DEBUG_PROC) (GLenum source,
                                            GLenum type,
                                            GLuint id,
                                            GLenum severity,
                                            GLsizei length,
                                            const gchar* message,
                                            gpointer user_data);
#endif
#endif

#endif
