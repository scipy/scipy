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

#ifndef __GST_GL_SHADER_STRINGS_H__
#define __GST_GL_SHADER_STRINGS_H__

#include <gst/gst.h>
#include <gst/gl/gl-prelude.h>
#include <gst/gl/gl.h>

G_BEGIN_DECLS

GST_GL_API
const gchar *gst_gl_shader_string_fragment_highp_precision;
GST_GL_API
const gchar *gst_gl_shader_string_fragment_mediump_precision;

GST_GL_API
const gchar *gst_gl_shader_string_vertex_default;
GST_GL_DEPRECATED_FOR(gst_gl_shader_string_fragment_get_default)
const gchar *gst_gl_shader_string_fragment_default;

GST_GL_API
const gchar *gst_gl_shader_string_vertex_mat4_texture_transform;
GST_GL_API
const gchar *gst_gl_shader_string_vertex_mat4_vertex_transform;
GST_GL_DEPRECATED_FOR(gst_gl_shader_string_fragment_external_oes_get_default)
const gchar *gst_gl_shader_string_fragment_external_oes_default;

GST_GL_API
const gchar * gst_gl_shader_string_get_highest_precision        (GstGLContext * context,
                                                                 GstGLSLVersion version,
                                                                 GstGLSLProfile profile);

GST_GL_API
gchar * gst_gl_shader_string_fragment_get_default               (GstGLContext * context,
                                                                 GstGLSLVersion version,
                                                                 GstGLSLProfile profile);
GST_GL_API
gchar * gst_gl_shader_string_fragment_external_oes_get_default  (GstGLContext * context,
                                                                 GstGLSLVersion version,
                                                                 GstGLSLProfile profile);

G_END_DECLS

#endif /* __GST_GL_SHADER_STRINGS_H__ */
