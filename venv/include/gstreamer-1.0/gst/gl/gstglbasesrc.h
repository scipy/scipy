/* 
 * GStreamer
 * Copyright (C) <1999> Erik Walthinsen <omega@cse.ogi.edu>
 * Copyright (C) 2002,2007 David A. Schleef <ds@schleef.org>
 * Copyright (C) 2008 Julien Isorce <julien.isorce@gmail.com>
 * Copyright (C) 2019 Philippe Normand <philn@igalia.com>
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

#ifndef __GST_GL_BASE_SRC_H__
#define __GST_GL_BASE_SRC_H__

#include <gst/base/gstpushsrc.h>
#include <gst/gl/gstgl_fwd.h>

G_BEGIN_DECLS

GST_GL_API
GType gst_gl_base_src_get_type(void);

#define GST_TYPE_GL_BASE_SRC            (gst_gl_base_src_get_type())
#define GST_GL_BASE_SRC(obj)            (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_GL_BASE_SRC,GstGLBaseSrc))
#define GST_GL_BASE_SRC_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_GL_BASE_SRC,GstGLBaseSrcClass))
#define GST_IS_GL_BASE_SRC(obj)         (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_GL_BASE_SRC))
#define GST_IS_GL_BASE_SRC_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_GL_BASE_SRC))
#define GST_GL_BASE_SRC_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS((obj) ,GST_TYPE_GL_BASE_SRC,GstGLBaseSrcClass))

/**
 * GstGLBaseSrc:
 * @display: the currently configured #GstGLDisplay
 * @context: the currently configured #GstGLContext
 * @out_caps: the currently configured output #GstCaps
 * @out_info: the currently configured output #GstVideoInfo
 * @running_time: the total running time
 *
 * The parent instance type of a base GStreamer GL Video source.
 *
 * Since: 1.18
 */
struct _GstGLBaseSrc {
  GstPushSrc parent;

  /*< public >*/
  GstGLDisplay *display;
  GstGLContext *context;

  /* video state */
  GstVideoInfo out_info;
  GstCaps *out_caps;

  /* total running time */
  GstClockTime running_time;

  /*< private >*/
  gpointer           _padding[GST_PADDING];

  GstGLBaseSrcPrivate *priv;
};

/**
 * GstGLBaseSrcClass:
 * @supported_gl_api: the logical-OR of #GstGLAPI's supported by this element
 * @gl_start: called in the GL thread to setup the element GL state.
 * @gl_stop: called in the GL thread to setup the element GL state.
 * @fill_gl_memory: called in the GL thread to fill the current video texture.
 *
 * The base class for GStreamer GL Video sources.
 *
 * Since: 1.18
 */
struct _GstGLBaseSrcClass {
  GstPushSrcClass parent_class;

  /*< public >*/
  GstGLAPI supported_gl_api;
  gboolean (*gl_start)          (GstGLBaseSrc *src);
  void     (*gl_stop)           (GstGLBaseSrc *src);
  gboolean (*fill_gl_memory)    (GstGLBaseSrc *src, GstGLMemory *mem);

  /*< private >*/
  gpointer _padding[GST_PADDING];
};

G_END_DECLS

#endif /* __GST_GL_BASE_SRC_H__ */
