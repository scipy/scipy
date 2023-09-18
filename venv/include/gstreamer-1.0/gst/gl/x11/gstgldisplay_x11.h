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

#ifndef __GST_GL_DISPLAY_X11_H__
#define __GST_GL_DISPLAY_X11_H__

#include <gst/gst.h>

#include <X11/Xlib-xcb.h>

#include <gst/gl/gstgldisplay.h>

G_BEGIN_DECLS

GST_GL_API
GType gst_gl_display_x11_get_type (void);

#define GST_TYPE_GL_DISPLAY_X11             (gst_gl_display_x11_get_type())
#define GST_GL_DISPLAY_X11(obj)             (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_GL_DISPLAY_X11,GstGLDisplayX11))
#define GST_GL_DISPLAY_X11_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST((klass), GST_TYPE_GL_DISPLAY_X11,GstGLDisplayX11Class))
#define GST_IS_GL_DISPLAY_X11(obj)          (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_GL_DISPLAY_X11))
#define GST_IS_GL_DISPLAY_X11_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE((klass), GST_TYPE_GL_DISPLAY_X11))
#define GST_GL_DISPLAY_X11_CAST(obj)        ((GstGLDisplayX11*)(obj))

typedef struct _GstGLDisplayX11 GstGLDisplayX11;
typedef struct _GstGLDisplayX11Class GstGLDisplayX11Class;

/**
 * GstGLDisplayX11:
 *
 * the contents of a #GstGLDisplayX11 are private and should only be accessed
 * through the provided API
 */
struct _GstGLDisplayX11
{
  /*< private >*/
  GstGLDisplay          parent;

  gchar *name;
  Display *display;
  xcb_connection_t *xcb_connection;
  gboolean foreign_display;

  gpointer _padding[GST_PADDING];
};

struct _GstGLDisplayX11Class
{
  GstGLDisplayClass object_class;

  gpointer _padding[GST_PADDING];
};

GST_GL_API
GstGLDisplayX11 *gst_gl_display_x11_new (const gchar * name);

GST_GL_API
GstGLDisplayX11 *gst_gl_display_x11_new_with_display (Display *display);

G_END_DECLS

#endif /* __GST_GL_DISPLAY_X11_H__ */
