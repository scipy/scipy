/*
 * GStreamer
 * Copyright (C) 2007 David Schleef <ds@schleef.org>
 * Copyright (C) 2008 Julien Isorce <julien.isorce@gmail.com>
 * Copyright (C) 2008 Filippo Argiolas <filippo.argiolas@gmail.com>
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

#ifndef _GST_GL_FILTER_H_
#define _GST_GL_FILTER_H_

#include <gst/gst.h>
#include <gst/video/video.h>

#include <gst/gl/gl.h>

G_BEGIN_DECLS

GST_GL_API
GType gst_gl_filter_get_type(void);
#define GST_TYPE_GL_FILTER            (gst_gl_filter_get_type())
#define GST_GL_FILTER(obj)            (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_GL_FILTER,GstGLFilter))
#define GST_IS_GL_FILTER(obj)         (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_GL_FILTER))
#define GST_GL_FILTER_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST((klass) ,GST_TYPE_GL_FILTER,GstGLFilterClass))
#define GST_IS_GL_FILTER_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass) ,GST_TYPE_GL_FILTER))
#define GST_GL_FILTER_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS((obj) ,GST_TYPE_GL_FILTER,GstGLFilterClass))

/**
 * GstGLFilterRenderFunc:
 * @filter: the #GstGLFilter
 * @in_tex: the input #GstGLMemory to render
 * @user_data: user data
 *
 * Returns: whether the render succeeded
 *
 * Since: 1.10
 */
typedef gboolean (*GstGLFilterRenderFunc) (GstGLFilter * filter, GstGLMemory * in_tex, gpointer user_data);

/**
 * GstGLFilter:
 * @in_info: the video info for input buffers
 * @out_info: the video info for output buffers
 * @in_texture_target: The texture target of the input buffers (usually 2D)
 * @out_texture_target: The texture target of the output buffers (usually 2D)
 * @out_caps: the output #GstCaps
 * @fbo: #GstGLFramebuffer object used for transformations (only for subclass usage)
 */
struct _GstGLFilter
{
  GstGLBaseFilter    parent;

  /*< public >*/
  GstVideoInfo       in_info;
  GstVideoInfo       out_info;
  GstGLTextureTarget in_texture_target;
  GstGLTextureTarget out_texture_target;

  GstCaps           *out_caps;

  /* protected */
  GstGLFramebuffer  *fbo;

  /*< private >*/
  gboolean           gl_result;
  GstBuffer         *inbuf;
  GstBuffer         *outbuf;

  GstGLShader       *default_shader;
  gboolean           valid_attributes;

  guint              vao;
  guint              vbo_indices;
  guint              vertex_buffer;
  gint               draw_attr_position_loc;
  gint               draw_attr_texture_loc;

  gpointer          _padding[GST_PADDING];
};

/**
 * GstGLFilterClass:
 * @set_caps: mirror from #GstBaseTransform
 * @filter: perform operations on the input and output buffers.  In general,
 *          you should avoid using this method if at all possible. One valid
 *          use-case for using this is keeping previous buffers for future calculations.
 *          Note: If @filter exists, then @filter_texture is not run
 * @filter_texture: given @in_tex, transform it into @out_tex.  Not used
 *                  if @filter exists
 * @init_fbo: perform initialization when the Framebuffer object is created
 * @transform_internal_caps: Perform sub-class specific modifications of the
 *   caps to be processed between upload on input and before download for output.
 */
struct _GstGLFilterClass
{
  GstGLBaseFilterClass parent_class;

  /*< public >*/
  gboolean (*set_caps)          (GstGLFilter* filter, GstCaps* incaps, GstCaps* outcaps);
  gboolean (*filter)            (GstGLFilter *filter, GstBuffer *inbuf, GstBuffer *outbuf);
  gboolean (*filter_texture)    (GstGLFilter *filter, GstGLMemory *input, GstGLMemory *output);
  gboolean (*init_fbo)          (GstGLFilter *filter);

  GstCaps *(*transform_internal_caps) (GstGLFilter *filter,
    GstPadDirection direction, GstCaps * caps, GstCaps * filter_caps);

  /*< private >*/
  gpointer                      _padding[GST_PADDING];
};

GST_GL_API
void gst_gl_filter_add_rgba_pad_templates (GstGLFilterClass *klass);

GST_GL_API
gboolean gst_gl_filter_filter_texture (GstGLFilter * filter, GstBuffer * input,
                                       GstBuffer * output);

GST_GL_API
gboolean gst_gl_filter_render_to_target             (GstGLFilter *filter,
                                                     GstGLMemory * input,
                                                     GstGLMemory * output,
                                                     GstGLFilterRenderFunc func,
                                                     gpointer data);

GST_GL_API
void gst_gl_filter_draw_fullscreen_quad             (GstGLFilter *filter);
GST_GL_API
void gst_gl_filter_render_to_target_with_shader     (GstGLFilter * filter,
                                                     GstGLMemory * input,
                                                     GstGLMemory * output,
                                                     GstGLShader *shader);

G_END_DECLS

#endif /* _GST_GL_FILTER_H_ */
