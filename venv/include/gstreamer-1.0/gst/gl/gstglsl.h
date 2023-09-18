/*
 * GStreamer
 * Copyright (C) 2015 Matthew Waters <matthew@centricular.com>
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

#ifndef __GST_GLSL_H__
#define __GST_GLSL_H__

#include <gst/gl/gstgl_fwd.h>

G_BEGIN_DECLS

GST_GL_API
GQuark gst_glsl_error_quark (void);

/**
 * GST_GLSL_ERROR:
 *
 * Error domain for GStreamer's GLSL module. Errors in this domain will be
 * from the #GstGLSLError enumeration
 */
#define GST_GLSL_ERROR (gst_glsl_error_quark ())

/**
 * GstGLSLError:
 * @GST_GLSL_ERROR_COMPILE: Compilation error occurred
 * @GST_GLSL_ERROR_LINK: Link error occurred
 * @GST_GLSL_ERROR_PROGRAM: General program error occurred
 *
 * Compilation stage that caused an error
 *
 * Since: 1.8
 */
typedef enum {
  GST_GLSL_ERROR_COMPILE,
  GST_GLSL_ERROR_LINK,
  GST_GLSL_ERROR_PROGRAM,
} GstGLSLError;

/**
 * GstGLSLVersion:
 * @GST_GLSL_VERSION_NONE: no version
 * @GST_GLSL_VERSION_100: version 100 (only valid for ES)
 * @GST_GLSL_VERSION_110: version 110 (only valid for compatibility desktop GL)
 * @GST_GLSL_VERSION_120: version 120 (only valid for compatibility desktop GL)
 * @GST_GLSL_VERSION_130: version 130 (only valid for compatibility desktop GL)
 * @GST_GLSL_VERSION_140: version 140 (only valid for compatibility desktop GL)
 * @GST_GLSL_VERSION_150: version 150 (valid for compatibility/core desktop GL)
 * @GST_GLSL_VERSION_300: version 300 (only valid for ES)
 * @GST_GLSL_VERSION_310: version 310 (only valid for ES)
 * @GST_GLSL_VERSION_320: version 320 (only valid for ES)
 * @GST_GLSL_VERSION_330: version 330 (valid for compatibility/core desktop GL)
 * @GST_GLSL_VERSION_400: version 400 (valid for compatibility/core desktop GL)
 * @GST_GLSL_VERSION_410: version 410 (valid for compatibility/core desktop GL)
 * @GST_GLSL_VERSION_420: version 420 (valid for compatibility/core desktop GL)
 * @GST_GLSL_VERSION_430: version 430 (valid for compatibility/core desktop GL)
 * @GST_GLSL_VERSION_440: version 440 (valid for compatibility/core desktop GL)
 * @GST_GLSL_VERSION_450: version 450 (valid for compatibility/core desktop GL)
 *
 * GLSL version list
 *
 * Since: 1.8
 */
typedef enum
{
  GST_GLSL_VERSION_NONE = 0,

  GST_GLSL_VERSION_100 = 100, /* ES */
  GST_GLSL_VERSION_110 = 110, /* GL */
  GST_GLSL_VERSION_120 = 120, /* GL */
  GST_GLSL_VERSION_130 = 130, /* GL */
  GST_GLSL_VERSION_140 = 140, /* GL */
  GST_GLSL_VERSION_150 = 150, /* GL */
  GST_GLSL_VERSION_300 = 300, /* ES */
  GST_GLSL_VERSION_310 = 310, /* ES */
  GST_GLSL_VERSION_320 = 320, /* ES */
  GST_GLSL_VERSION_330 = 330, /* GL */
  GST_GLSL_VERSION_400 = 400, /* GL */
  GST_GLSL_VERSION_410 = 410, /* GL */
  GST_GLSL_VERSION_420 = 420, /* GL */
  GST_GLSL_VERSION_430 = 430, /* GL */
  GST_GLSL_VERSION_440 = 440, /* GL */
  GST_GLSL_VERSION_450 = 450, /* GL */
} GstGLSLVersion;

/**
 * GstGLSLProfile:
 * @GST_GLSL_PROFILE_NONE: no profile supported/available
 * @GST_GLSL_PROFILE_ES: OpenGL|ES profile
 * @GST_GLSL_PROFILE_CORE: OpenGL core profile
 * @GST_GLSL_PROFILE_COMPATIBILITY: OpenGL compatibility profile
 * @GST_GLSL_PROFILE_ANY: any OpenGL/OpenGL|ES profile
 *
 * GLSL profiles
 *
 * Since: 1.8
 */
/* FIXME: For GST_GLSL_PROFILE_ANY ~0 -> 0xffffffff see
 * https://bugzilla.gnome.org/show_bug.cgi?id=732633
*/
typedef enum
{
  /* XXX: maybe make GstGLAPI instead */
  GST_GLSL_PROFILE_NONE = 0,

  GST_GLSL_PROFILE_ES = (1 << 0),
  GST_GLSL_PROFILE_CORE = (1 << 1),
  GST_GLSL_PROFILE_COMPATIBILITY = (1 << 2),

  GST_GLSL_PROFILE_ANY = (gint) (0xffffffff),
} GstGLSLProfile;

GST_GL_API
GstGLSLVersion gst_glsl_version_from_string         (const gchar * string);
GST_GL_API
const gchar *  gst_glsl_version_to_string           (GstGLSLVersion version);

GST_GL_API
GstGLSLProfile gst_glsl_profile_from_string         (const gchar * string);
GST_GL_API
const gchar *  gst_glsl_profile_to_string           (GstGLSLProfile profile);

GST_GL_API
gchar *        gst_glsl_version_profile_to_string   (GstGLSLVersion version,
                                                     GstGLSLProfile profile);
GST_GL_API
gboolean       gst_glsl_version_profile_from_string (const gchar * string,
                                                     GstGLSLVersion * version_ret,
                                                     GstGLSLProfile * profile_ret);

GST_GL_API
gboolean       gst_glsl_string_get_version_profile  (const gchar *s,
                                                     GstGLSLVersion * version,
                                                     GstGLSLProfile * profile);

GST_GL_API
GstGLSLVersion gst_gl_version_to_glsl_version       (GstGLAPI gl_api, gint maj, gint min);
GST_GL_API
gboolean       gst_gl_context_supports_glsl_profile_version (GstGLContext * context,
                                                             GstGLSLVersion version,
                                                             GstGLSLProfile profile);

GST_GL_API
gboolean       gst_gl_context_supports_precision        (GstGLContext * context,
                                                         GstGLSLVersion version,
                                                         GstGLSLProfile profile);
GST_GL_API
gboolean       gst_gl_context_supports_precision_highp  (GstGLContext * context,
                                                         GstGLSLVersion version,
                                                         GstGLSLProfile profile);

G_END_DECLS

#endif /* __GST_GLSL_H__ */
