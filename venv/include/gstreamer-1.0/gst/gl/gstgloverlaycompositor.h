/*
 * GStreamer
 * Copyright (C) 2015 Lubosz Sarnecki <lubosz.sarnecki@collabora.co.uk>
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

#ifndef __GST_GL_OVERLAY_COMPOSITOR_H__
#define __GST_GL_OVERLAY_COMPOSITOR_H__

#include <gst/video/video.h>
#include <gst/gl/gstgl_fwd.h>

#define GST_TYPE_GL_OVERLAY_COMPOSITOR (gst_gl_overlay_compositor_get_type())
#define GST_GL_OVERLAY_COMPOSITOR(obj) (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_GL_OVERLAY_COMPOSITOR,GstGLOverlayCompositor))
#define GST_GL_OVERLAY_COMPOSITOR_CLASS(klass) (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_GL_OVERLAY_COMPOSITOR,GstGLOverlayCompositorClass))
#define GST_IS_GL_OVERLAY_COMPOSITOR(obj) (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_GL_OVERLAY_COMPOSITOR))
#define GST_IS_GL_OVERLAY_COMPOSITOR_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_GL_OVERLAY_COMPOSITOR))
#define GST_GL_OVERLAY_COMPOSITOR_CAST(obj) ((GstGLOverlayCompositor*)(obj))

G_BEGIN_DECLS

GST_GL_API
GType gst_gl_overlay_compositor_get_type (void);

/**
 * GstGLOverlayCompositor:
 *
 * Opaque #GstGLOverlayCompositor object
 */
struct _GstGLOverlayCompositor
{
  /*< private >*/
  GstObject parent;

  GstGLContext *context;

  guint last_window_width;
  guint last_window_height;

  GList * overlays;
 
  GstGLShader *shader;
  gint  position_attrib;
  gint  texcoord_attrib;

  gpointer _padding[GST_PADDING];
};

/**
 * GstGLOverlayCompositorClass:
 *
 */
struct _GstGLOverlayCompositorClass
{
  GstObjectClass object_class;

  /*< private >*/
  gpointer _padding[GST_PADDING];
};

GST_GL_API
GstGLOverlayCompositor *gst_gl_overlay_compositor_new (GstGLContext * context);

GST_GL_API
void gst_gl_overlay_compositor_free_overlays (GstGLOverlayCompositor * compositor);

GST_GL_API
void gst_gl_overlay_compositor_upload_overlays (GstGLOverlayCompositor * compositor,
        GstBuffer * buf);

GST_GL_API
void gst_gl_overlay_compositor_draw_overlays (GstGLOverlayCompositor * compositor);

GST_GL_API
GstCaps * gst_gl_overlay_compositor_add_caps(GstCaps * caps);

G_END_DECLS
#endif /* __GST_GL_OVERLAY_COMPOSITOR_H__ */
