/*
 * GStreamer
 * Copyright (C) 2013 Matthew Waters <ystreet00@gmail.com>
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

#ifndef __GST_GL_UTILS_H__
#define __GST_GL_UTILS_H__

#include <gst/video/video.h>
#include <gst/video/gstvideoaffinetransformationmeta.h>

#include <gst/gl/gstgl_fwd.h>

G_BEGIN_DECLS

GST_GL_API
gboolean gst_gl_ensure_element_data (gpointer element,
    GstGLDisplay **display_ptr, GstGLContext ** other_context_ptr);
GST_GL_API
gboolean gst_gl_handle_set_context (GstElement * element, GstContext * context,
    GstGLDisplay ** display, GstGLContext ** other_context);
GST_GL_API
gboolean gst_gl_handle_context_query (GstElement * element, GstQuery * query,
    GstGLDisplay * display, GstGLContext * context, GstGLContext * other_context);
GST_GL_API
gboolean gst_gl_query_local_gl_context (GstElement * element, GstPadDirection direction,
    GstGLContext ** context_ptr);

GST_GL_API
void gst_gl_element_propagate_display_context (GstElement * element, GstGLDisplay * display);

GST_GL_API
gsize gst_gl_get_plane_data_size (const GstVideoInfo * info, const GstVideoAlignment * align,
    guint plane);
GST_GL_API
gsize gst_gl_get_plane_start (const GstVideoInfo * info, const GstVideoAlignment * valign,
    guint plane);

GST_GL_API
gboolean gst_gl_value_set_texture_target_from_mask (GValue * value,
    GstGLTextureTarget target_mask);
GST_GL_API
gboolean gst_gl_value_set_texture_target (GValue * value, GstGLTextureTarget target);
GST_GL_API
GstGLTextureTarget gst_gl_value_get_texture_target_mask (const GValue * value);

GST_GL_API
void gst_gl_get_affine_transformation_meta_as_ndc (GstVideoAffineTransformationMeta * meta, gfloat * matrix);
GST_GL_API
void gst_gl_set_affine_transformation_meta_from_ndc (GstVideoAffineTransformationMeta * meta, const gfloat * matrix);

GST_GL_API
void gst_gl_multiply_matrix4 (const gfloat * a, const gfloat * b, gfloat * result);


G_END_DECLS

#endif /* __GST_GL_UTILS_H__ */
