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

#ifndef __GST_GL_FUNCS_H__
#define __GST_GL_FUNCS_H__

#include <gst/gl/gstglconfig.h>

/* This mimic GCC behaviour with system headers files even if GL headers may
 * not be in the system header path. */
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wredundant-decls"
#endif

/* OpenGL 2.0 for Embedded Systems */
#if GST_GL_HAVE_GLES2
# if GST_GL_HAVE_PLATFORM_EAGL
#  include <OpenGLES/ES2/gl.h>
#  include <OpenGLES/ES2/glext.h>
# else
#  if GST_GL_HAVE_GLES3
#   include <GLES3/gl3.h>
#   if GST_GL_HAVE_GLES3EXT3_H
#    include <GLES3/gl3ext.h>
#   endif
#   include <GLES2/gl2ext.h>
#  else
#   include <GLES2/gl2.h>
#   include <GLES2/gl2ext.h>
#  endif
# endif
# if !GST_GL_HAVE_OPENGL
#  include <gst/gl/glprototypes/gstgl_gles2compat.h>
# endif
#endif

/* OpenGL for desktop systems */
#if GST_GL_HAVE_OPENGL
# ifdef __APPLE__
#  if !defined(MAC_OS_X_VERSION_MAX_ALLOWED) || MAC_OS_X_VERSION_MAX_ALLOWED >= 1014
#   define GL_SILENCE_DEPRECATION
#  endif
#  include <OpenGL/OpenGL.h>
#  include <OpenGL/gl.h>
#  if MAC_OS_X_VERSION_MAX_ALLOWED >= 1070
#   define GL_DO_NOT_WARN_IF_MULTI_GL_VERSION_HEADERS_INCLUDED
#   include <OpenGL/gl3.h>
#  endif
# else
#  if defined(_MSC_VER)
#   include <windows.h>
#  endif
#  include <GL/gl.h>
#  if defined(__WIN32__) || defined(_WIN32)
#   include <GL/glext.h>
#  endif
# endif
#endif

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#if defined(WINAPI)
#define GSTGLAPI WINAPI
#else
#define GSTGLAPI
#endif
#include <gst/gl/glprototypes/gstgl_compat.h>

#include <gst/gst.h>

G_BEGIN_DECLS

#define GST_GL_EXT_BEGIN(name, gl_availability, min_gl, maj_gl, gles_maj, \
    gles_min, ext_suf, ext_name)
#define GST_GL_EXT_FUNCTION(ret, name, args) \
  ret (GSTGLAPI *name) args;
#define GST_GL_EXT_END()

struct _GstGLFuncs
{
#include <gst/gl/glprototypes/all_functions.h>
  gpointer padding[GST_PADDING_LARGE*6-2];
};

#undef GST_GL_EXT_BEGIN
#undef GST_GL_EXT_FUNCTION
#undef GST_GL_EXT_END

G_END_DECLS

#endif /* __GST_GL_API_H__ */
