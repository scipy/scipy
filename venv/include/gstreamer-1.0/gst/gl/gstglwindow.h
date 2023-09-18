/*
 * GStreamer
 * Copyright (C) 2008 Julien Isorce <julien.isorce@gmail.com>
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

#ifndef __GST_GL_WINDOW_H__
#define __GST_GL_WINDOW_H__

#include <gst/gst.h>

#include <gst/gl/gstgl_fwd.h>
#include <gst/gl/gstglcontext.h>
#include <gst/gl/gstgldisplay.h>

G_BEGIN_DECLS

GST_GL_API
GType gst_gl_window_get_type       (void);
#define GST_TYPE_GL_WINDOW         (gst_gl_window_get_type())

#define GST_GL_WINDOW(o)           (G_TYPE_CHECK_INSTANCE_CAST((o), GST_TYPE_GL_WINDOW, GstGLWindow))
#define GST_GL_WINDOW_CLASS(k)     (G_TYPE_CHECK_CLASS_CAST((k), GST_TYPE_GL_WINDOW, GstGLWindowClass))
#define GST_IS_GL_WINDOW(o)        (G_TYPE_CHECK_INSTANCE_TYPE((o), GST_TYPE_GL_WINDOW))
#define GST_IS_GL_WINDOW_CLASS(k)  (G_TYPE_CHECK_CLASS_TYPE((k), GST_TYPE_GL_WINDOW))
#define GST_GL_WINDOW_GET_CLASS(o) (G_TYPE_INSTANCE_GET_CLASS((o), GST_TYPE_GL_WINDOW, GstGLWindowClass))

#define GST_GL_WINDOW_LOCK(w) g_mutex_lock(&GST_GL_WINDOW(w)->lock)
#define GST_GL_WINDOW_UNLOCK(w) g_mutex_unlock(&GST_GL_WINDOW(w)->lock)
#define GST_GL_WINDOW_GET_LOCK(w) (&GST_GL_WINDOW(w)->lock)

GST_GL_API
GQuark gst_gl_window_error_quark (void);
/**
 * GST_GL_WINDOW_ERROR:
 *
 * Error domain for GStreamer's GL window module. Errors in this domain will be
 * from the #GstGLWindowError enumeration
 */
#define GST_GL_WINDOW_ERROR (gst_gl_window_error_quark ())

/**
 * GstGLWindowError:
 * @GST_GL_WINDOW_ERROR_FAILED: failed for a unspecified reason
 * @GST_GL_WINDOW_ERROR_OLD_LIBS: the implementation is too old
 * @GST_GL_WINDOW_ERROR_RESOURCE_UNAVAILABLE: no such resource was found
 */
typedef enum
{
  GST_GL_WINDOW_ERROR_FAILED,
  GST_GL_WINDOW_ERROR_OLD_LIBS,
  GST_GL_WINDOW_ERROR_RESOURCE_UNAVAILABLE,
} GstGLWindowError;

typedef void (*GstGLWindowCB) (gpointer data);
typedef void (*GstGLWindowResizeCB) (gpointer data, guint width, guint height);

/**
 * GST_GL_WINDOW_CB:
 * @f: the function to cast
 *
 * Cast to the current function type for generic window callbacks
 */
#define	GST_GL_WINDOW_CB(f)			 ((GstGLWindowCB) (f))

/**
 * GST_GL_WINDOW_RESIZE_CB:
 * @f: the function to cast
 *
 * Cast to the current function type for window resize callbacks
 */
#define	GST_GL_WINDOW_RESIZE_CB(f)		 ((GstGLWindowResizeCB) (f))

/**
 * GstGLWindow:
 *
 * #GstGLWindow is an opaque struct and should only be accessed through the
 * provided api.
 */
struct _GstGLWindow {
  /*< private >*/
  GstObject parent;

  GMutex        lock;

  GstGLDisplay *display;
  GWeakRef      context_ref;

  /*< protected >*/
  gboolean      is_drawing;

  GstGLWindowCB         draw;
  gpointer              draw_data;
  GDestroyNotify        draw_notify;
  GstGLWindowCB         close;
  gpointer              close_data;
  GDestroyNotify        close_notify;
  GstGLWindowResizeCB   resize;
  gpointer              resize_data;
  GDestroyNotify        resize_notify;

  gboolean              queue_resize;

  GMainContext         *main_context; /* default main_context */

  /*< private >*/
  GstGLWindowPrivate *priv;

  gpointer _reserved[GST_PADDING];
};

/**
 * GstGLWindowClass:
 * @parent_class: Parent class
 * @get_display: Gets the current windowing system display connection
 * @set_window_handle: Set a window handle to render into
 * @get_window_handle: Gets the current window handle that this #GstGLWindow is
 *                     rendering into.  This may return a different value to
 *                     what is passed into @set_window_handle
 * @draw: redraw the window with the specified dimensions
 * @run: run the mainloop
 * @quit: send a quit to the mainloop
 * @send_message: invoke a function on the window thread.  Required to be reentrant.
 * @send_message_async: invoke a function on the window thread. @run may or may
 *                      not have been called.  Required to be reentrant.
 * @open: open the connection to the display
 * @close: close the connection to the display
 * @handle_events: whether to handle 'extra' events from the windowing system.
 *                 Basic events like surface moves and resizes are still valid
 *                 things to listen for.
 * @set_preferred_size: request that the window change surface size.  The
 *                      implementation is free to ignore this information.
 * @show: request that the window be shown to the user
 * @set_render_rectangle: request a rectangle to render into.  See #GstVideoOverlay
 * @queue_resize: request a resize to occur when possible
 * @controls_viewport: Whether the window takes care of glViewport setup.
 *                     and the user does not need to deal with viewports
 * @has_output_surface: Whether the window has output surface or not. (Since: 1.18)
 */
struct _GstGLWindowClass {
  GstObjectClass parent_class;

  guintptr (*get_display)        (GstGLWindow *window);
  void     (*set_window_handle)  (GstGLWindow *window, guintptr handle);
  guintptr (*get_window_handle)  (GstGLWindow *window);
  void     (*draw)               (GstGLWindow *window);
  void     (*run)                (GstGLWindow *window);
  void     (*quit)               (GstGLWindow *window);
  void     (*send_message)       (GstGLWindow *window, GstGLWindowCB callback, gpointer data);
  void     (*send_message_async) (GstGLWindow *window, GstGLWindowCB callback, gpointer data, GDestroyNotify destroy);

  gboolean (*open)               (GstGLWindow *window, GError **error);
  void     (*close)              (GstGLWindow *window);
  void     (*handle_events)      (GstGLWindow *window, gboolean handle_events);
  void     (*set_preferred_size) (GstGLWindow *window, gint width, gint height);
  void     (*show)               (GstGLWindow *window);
  gboolean (*set_render_rectangle)(GstGLWindow *window, gint x, gint y, gint width, gint height);
  void     (*queue_resize)       (GstGLWindow *window);
  gboolean (*controls_viewport)  (GstGLWindow *window);
  gboolean (*has_output_surface) (GstGLWindow *window);

  /*< private >*/
  gpointer _reserved[GST_PADDING-2];
};

GST_GL_API
GstGLWindow * gst_gl_window_new  (GstGLDisplay *display);

/* callbacks */
GST_GL_API
void     gst_gl_window_set_draw_callback    (GstGLWindow *window,
                                             GstGLWindowCB callback,
                                             gpointer data,
                                             GDestroyNotify destroy_notify);
GST_GL_API
void     gst_gl_window_set_resize_callback  (GstGLWindow *window,
                                             GstGLWindowResizeCB callback,
                                             gpointer data,
                                             GDestroyNotify destroy_notify);
GST_GL_API
void     gst_gl_window_set_close_callback   (GstGLWindow *window,
                                             GstGLWindowCB callback,
                                             gpointer data,
                                             GDestroyNotify destroy_notify);

GST_GL_API
void     gst_gl_window_set_window_handle    (GstGLWindow *window, guintptr handle);
GST_GL_API
guintptr gst_gl_window_get_window_handle    (GstGLWindow *window);

/* loop/events */
GST_GL_API
void     gst_gl_window_run                  (GstGLWindow *window);
GST_GL_API
void     gst_gl_window_quit                 (GstGLWindow *window);
GST_GL_API
void     gst_gl_window_send_message         (GstGLWindow *window,
                                             GstGLWindowCB callback,
                                             gpointer data);
GST_GL_API
void     gst_gl_window_send_message_async   (GstGLWindow *window,
                                             GstGLWindowCB callback,
                                             gpointer data,
                                             GDestroyNotify destroy);

/* navigation */
GST_GL_API
void     gst_gl_window_handle_events        (GstGLWindow * window,
                                             gboolean handle_events);

GST_GL_API
void     gst_gl_window_send_key_event       (GstGLWindow * window,
                                             const char * event_type,
                                             const char * key_str);
GST_GL_API
void     gst_gl_window_send_mouse_event     (GstGLWindow * window,
                                             const char * event_type,
                                             int button,
                                             double posx,
                                             double posy);

GST_GL_API
void     gst_gl_window_send_scroll_event    (GstGLWindow * window,
                                             double posx,
                                             double posy,
                                             double delta_x,
                                             double delta_y);

/* surfaces/rendering */
GST_GL_API
void     gst_gl_window_queue_resize         (GstGLWindow *window);
GST_GL_API
void     gst_gl_window_draw                 (GstGLWindow *window);
GST_GL_API
void     gst_gl_window_show                 (GstGLWindow *window);
GST_GL_API
void     gst_gl_window_set_preferred_size   (GstGLWindow * window,
                                             gint width,
                                             gint height);
GST_GL_API
void     gst_gl_window_get_surface_dimensions (GstGLWindow * window,
                                               guint * width,
                                               guint * height);
GST_GL_API
gboolean gst_gl_window_set_render_rectangle   (GstGLWindow * window,
                                               gint x,
                                               gint y,
                                               gint width,
                                               gint height);
GST_GL_API
gboolean gst_gl_window_controls_viewport      (GstGLWindow * window);

/* subclass usage only */
GST_GL_API
void     gst_gl_window_resize               (GstGLWindow *window, guint width, guint height);

GST_GL_API
GstGLContext * gst_gl_window_get_context    (GstGLWindow *window);
GST_GL_API
guintptr       gst_gl_window_get_display    (GstGLWindow *window);

GST_GL_API
gboolean       gst_gl_window_has_output_surface (GstGLWindow *window);

G_END_DECLS

#endif /* __GST_GL_WINDOW_H__ */
