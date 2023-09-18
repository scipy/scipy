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

#ifndef __GST_GL_CONTEXT_H__
#define __GST_GL_CONTEXT_H__

#include <gst/gst.h>

#include <gst/gl/gstgl_fwd.h>

G_BEGIN_DECLS

GST_GL_API
GType gst_gl_context_get_type       (void);
#define GST_TYPE_GL_CONTEXT         (gst_gl_context_get_type())

#define GST_GL_CONTEXT(o)           (G_TYPE_CHECK_INSTANCE_CAST((o), GST_TYPE_GL_CONTEXT, GstGLContext))
#define GST_GL_CONTEXT_CLASS(k)     (G_TYPE_CHECK_CLASS_CAST((k), GST_TYPE_GL_CONTEXT, GstGLContextClass))
#define GST_IS_GL_CONTEXT(o)        (G_TYPE_CHECK_INSTANCE_TYPE((o), GST_TYPE_GL_CONTEXT))
#define GST_IS_GL_CONTEXT_CLASS(k)  (G_TYPE_CHECK_CLASS_TYPE((k), GST_TYPE_GL_CONTEXT))
#define GST_GL_CONTEXT_GET_CLASS(o) (G_TYPE_INSTANCE_GET_CLASS((o), GST_TYPE_GL_CONTEXT, GstGLContextClass))

GST_GL_API
GQuark gst_gl_context_error_quark (void);

/**
 * GST_GL_CONTEXT_ERROR:
 *
 * Error domain for GStreamer's GL context module. Errors in this domain will
 * be from the #GstGLContextError enumeration
 */
#define GST_GL_CONTEXT_ERROR (gst_gl_context_error_quark ())

/**
 * GstGLContextThreadFunc:
 * @context: a #GstGLContext
 * @data: user data
 *
 * Represents a function to run in the GL thread with @context and @data
 */
typedef void (*GstGLContextThreadFunc) (GstGLContext * context, gpointer data);

#define GST_GL_CONTEXT_TYPE_CGL "gst.gl.context.CGL"
#define GST_GL_CONTEXT_TYPE_GLX "gst.gl.context.GLX"
#define GST_GL_CONTEXT_TYPE_EGL "gst.gl.context.EGL"
#define GST_GL_CONTEXT_TYPE_WGL "gst.gl.context.WGL"
#define GST_GL_CONTEXT_TYPE_EAGL "gst.gl.context.EAGL"

/**
 * GstGLContextError:
 * @GST_GL_CONTEXT_ERROR_FAILED: Failed for an unspecified reason
 * @GST_GL_CONTEXT_ERROR_WRONG_CONFIG: The configuration requested is not correct
 * @GST_GL_CONTEXT_ERROR_WRONG_API: The OpenGL API requested is not correct
 * @GST_GL_CONTEXT_ERROR_OLD_LIBS: The OpenGL libraries are too old
 * @GST_GL_CONTEXT_ERROR_CREATE_CONTEXT: glXCreateContext (or similar) failed
 * @GST_GL_CONTEXT_ERROR_RESOURCE_UNAVAILABLE: A resource is not available
 *
 * OpenGL context errors.
 */
typedef enum
{
  GST_GL_CONTEXT_ERROR_FAILED,
  GST_GL_CONTEXT_ERROR_WRONG_CONFIG,
  GST_GL_CONTEXT_ERROR_WRONG_API,
  GST_GL_CONTEXT_ERROR_OLD_LIBS,
  GST_GL_CONTEXT_ERROR_CREATE_CONTEXT,
  GST_GL_CONTEXT_ERROR_RESOURCE_UNAVAILABLE,
} GstGLContextError;

/**
 * GstGLContext:
 * @gl_vtable: a list of OpenGL function pointers
 *
 * Opaque #GstGLContext object
 */
struct _GstGLContext {
  /*< private >*/
  GstObject parent;

  GstGLDisplay *display;
  GstGLWindow  *window;

  /*< public >*/
  GstGLFuncs *gl_vtable;

  /*< private >*/
  GstGLContextPrivate *priv;

  gpointer _reserved[GST_PADDING];
};

/**
 * GstGLContextClass:
 * @get_gl_context: get the backing platform specific OpenGL context
 * @get_gl_api: get the available OpenGL api's that this context can work with
 * @get_proc_address: get an function pointer to an OpenGL function
 * @activate: call eglMakeCurrent or similar
 * @choose_format: choose a format for the framebuffer
 * @create_context: create the OpenGL context
 * @destroy_context: destroy the OpenGL context
 * @swap_buffers: swap the default framebuffer's front/back buffers
 */
/**
 * GstGLContextClass::get_config:
 * @context: the #GstGLContext
 *
 * Retrieve the configuration in use by this context.  See also
 * gst_gl_context_get_config().
 *
 * Returns: (transfer full) (nullable): the configuration chosen for this
 *          #GstGLContext
 *
 * Since: 1.20
 */
/**
 * GstGLContextClass::request_config:
 * @context: the #GstGLContext
 * @gl_config: (nullable) (transfer full): a configuration structure for
 *             configuring on @context
 *
 * Request a configuration for this @context to use.
 *
 * Unknown fields within @gl_config should be ignored by subclasses.
 *
 * See also gst_gl_context_request_config().
 *
 * Returns: Whether @gl_config could be successfull set on @context.
 *
 * Since: 1.20
 */
struct _GstGLContextClass {
  GstObjectClass parent_class;

  guintptr      (*get_current_context) (void);
  guintptr      (*get_gl_context)     (GstGLContext *context);
  GstGLAPI      (*get_gl_api)         (GstGLContext *context);
  GstGLPlatform (*get_gl_platform)    (GstGLContext *context);
  gpointer      (*get_proc_address)   (GstGLAPI gl_api, const gchar *name);
  gboolean      (*activate)           (GstGLContext *context, gboolean activate);
  gboolean      (*choose_format)      (GstGLContext *context, GError **error);
  gboolean      (*create_context)     (GstGLContext *context, GstGLAPI gl_api,
                                       GstGLContext *other_context, GError ** error);
  void          (*destroy_context)    (GstGLContext *context);
  void          (*swap_buffers)       (GstGLContext *context);
  gboolean      (*check_feature)      (GstGLContext *context, const gchar *feature);
  void          (*get_gl_platform_version) (GstGLContext *context, gint *major, gint *minor);
  GstStructure *(*get_config)         (GstGLContext * context);
  gboolean      (*request_config)     (GstGLContext * context, GstStructure * gl_config);

  /*< private >*/
  gpointer _reserved[GST_PADDING-2];
};

/* methods */

GST_GL_API
GstGLContext * gst_gl_context_new  (GstGLDisplay *display);
GST_GL_API
GstGLContext * gst_gl_context_new_wrapped (GstGLDisplay *display,
                                           guintptr handle,
                                           GstGLPlatform context_type,
                                           GstGLAPI available_apis);

GST_GL_API
GstStructure * gst_gl_context_get_config      (GstGLContext * context);
GST_GL_API
gboolean      gst_gl_context_request_config   (GstGLContext * context, GstStructure * gl_config);

GST_GL_API
gboolean      gst_gl_context_activate         (GstGLContext *context, gboolean activate);
GST_GL_API
GThread *     gst_gl_context_get_thread       (GstGLContext *context);
GST_GL_API
GstGLContext * gst_gl_context_get_current     (void);

GST_GL_API
GstGLDisplay * gst_gl_context_get_display (GstGLContext *context);
GST_GL_API
gpointer      gst_gl_context_get_proc_address (GstGLContext *context, const gchar *name);
GST_GL_API
GstGLPlatform gst_gl_context_get_gl_platform  (GstGLContext *context);
GST_GL_API
GstGLAPI      gst_gl_context_get_gl_api       (GstGLContext *context);
GST_GL_API
guintptr      gst_gl_context_get_gl_context   (GstGLContext *context);
GST_GL_API
gboolean      gst_gl_context_can_share        (GstGLContext * context, GstGLContext *other_context);
GST_GL_API
void          gst_gl_context_swap_buffers     (GstGLContext * context);

GST_GL_API
gboolean      gst_gl_context_create           (GstGLContext *context, GstGLContext *other_context, GError ** error);
GST_GL_API
void          gst_gl_context_destroy          (GstGLContext *context);

GST_GL_API
gpointer      gst_gl_context_default_get_proc_address (GstGLAPI gl_api, const gchar *name);
GST_GL_API
gpointer      gst_gl_context_get_proc_address_with_platform (GstGLPlatform context_type, GstGLAPI gl_api, const gchar *name);

GST_GL_API
gboolean      gst_gl_context_set_window (GstGLContext *context, GstGLWindow *window);
GST_GL_API
GstGLWindow * gst_gl_context_get_window (GstGLContext *context);

GST_GL_API
void          gst_gl_context_get_gl_version (GstGLContext *context, gint *maj, gint *min);
GST_GL_API
gboolean      gst_gl_context_check_gl_version (GstGLContext * context, GstGLAPI api, gint maj, gint min);
GST_GL_API
gboolean      gst_gl_context_check_feature (GstGLContext *context, const gchar *feature);
GST_GL_API
void          gst_gl_context_get_gl_platform_version (GstGLContext * context, gint * major, gint * minor);

GST_GL_API
guintptr      gst_gl_context_get_current_gl_context     (GstGLPlatform context_type);
GST_GL_API
GstGLAPI      gst_gl_context_get_current_gl_api         (GstGLPlatform platform, guint *major, guint *minor);

GST_GL_API
gboolean      gst_gl_context_is_shared                  (GstGLContext * context);
GST_GL_API
void          gst_gl_context_set_shared_with            (GstGLContext * context, GstGLContext * share);

GST_GL_API
gboolean gst_gl_context_fill_info (GstGLContext * context, GError ** error);

/* FIXME: remove */
GST_GL_API
void gst_gl_context_thread_add (GstGLContext * context,
    GstGLContextThreadFunc func, gpointer data);

G_END_DECLS

#endif /* __GST_GL_CONTEXT_H__ */
