/*
 * GStreamer
 * Copyright (C) 2016 Matthew Waters <matthew@centricular.com>
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

#ifndef __GST_GL_QUERY_H__
#define __GST_GL_QUERY_H__

#include <gst/gl/gstgl_fwd.h>
#include <gst/gl/gstgldebug.h>

G_BEGIN_DECLS

/**
 * GstGLQueryType:
 * @GST_GL_QUERY_NONE: no query
 * @GST_GL_QUERY_TIME_ELAPSED: query the time elapsed
 * @GST_GL_QUERY_TIMESTAMP: query the current time
 */
typedef enum
{
  GST_GL_QUERY_NONE,
  GST_GL_QUERY_TIME_ELAPSED,
  GST_GL_QUERY_TIMESTAMP,
} GstGLQueryType;

/**
 * GstGLQuery:
 *
 * Opaque #GstGLQuery struct
 */
struct _GstGLQuery
{
  /*< private >*/
  GstGLContext *    context;
  guint             query_type;
  guint             query_id;
  gboolean          supported;

  gboolean          start_called;
  GstGLAsyncDebug   debug;

  gpointer          _padding[GST_PADDING];
};

GST_GL_API
void                gst_gl_query_init               (GstGLQuery * query,
                                                     GstGLContext * context,
                                                     GstGLQueryType query_type);
GST_GL_API
void                gst_gl_query_unset              (GstGLQuery * query);
GST_GL_API
GstGLQuery *        gst_gl_query_new                (GstGLContext * context,
                                                     GstGLQueryType query_type);
GST_GL_API
void                gst_gl_query_free               (GstGLQuery * query);

GST_GL_API
void                gst_gl_query_start              (GstGLQuery * query);
GST_GL_API
void                gst_gl_query_end                (GstGLQuery * query);
GST_GL_API
void                gst_gl_query_counter            (GstGLQuery * query);
GST_GL_API
guint64             gst_gl_query_result             (GstGLQuery * query);

#define gst_gl_query_start_log_valist(query,cat,level,object,format,varargs) \
  G_STMT_START {    \
    GST_GL_ASYNC_CAT_LEVEL_LOG_valist (&(query)->debug, cat, level, object, format, varargs); \
    gst_gl_async_debug_freeze (&(query)->debug); \
    gst_gl_query_start (query); \
    gst_gl_async_debug_thaw (&(query)->debug); \
  } G_STMT_END

#define gst_gl_query_counter_log_valist(query,cat,level,object,format,varargs) \
  G_STMT_START {    \
    GST_GL_ASYNC_CAT_LEVEL_LOG_valist (&(query)->debug, cat, level, object, format, varargs); \
    gst_gl_async_debug_freeze (&(query)->debug); \
    gst_gl_query_counter (query); \
    gst_gl_async_debug_thaw (&(query)->debug); \
  } G_STMT_END

#ifdef G_HAVE_ISO_VARARGS

#define gst_gl_query_start_log(query,cat,level,object,format,...) \
  G_STMT_START {    \
    GST_GL_ASYNC_CAT_LEVEL_LOG (&(query)->debug, cat, level, object, format, __VA_ARGS__); \
    gst_gl_async_debug_freeze (&(query)->debug); \
    gst_gl_query_start (query); \
    gst_gl_async_debug_thaw (&(query)->debug); \
  } G_STMT_END
#define gst_gl_query_counter_log(query,cat,level,object,format,...) \
  G_STMT_START {    \
    GST_GL_ASYNC_CAT_LEVEL_LOG (&(query)->debug, cat, level, object, format, __VA_ARGS__); \
    gst_gl_async_debug_freeze (&(query)->debug); \
    gst_gl_query_counter (query); \
    gst_gl_async_debug_thaw (&(query)->debug); \
  } G_STMT_END

#else /* G_HAVE_ISO_VARARGS */
#if G_HAVE_GNUC_VARARGS

#define gst_gl_query_start_log(query,cat,level,object,format,args...) \
  G_STMT_START {    \
    GST_GL_ASYNC_CAT_LEVEL_LOG (&(query)->debug, cat, level, object, format, ##args); \
    gst_gl_async_debug_freeze (&(query)->debug); \
    gst_gl_query_start (query); \
    gst_gl_async_debug_thaw (&(query)->debug); \
  } G_STMT_END
#define gst_gl_query_counter_log(query,cat,level,object,format,args...) \
  G_STMT_START {    \
    GST_GL_ASYNC_CAT_LEVEL_LOG (&(query)->debug, cat, level, object, format, ##args); \
    gst_gl_async_debug_freeze (&(query)->debug); \
    gst_gl_query_counter (query); \
    gst_gl_async_debug_thaw (&(query)->debug); \
  } G_STMT_END

#else /* G_HAVE_GNUC_VARARGS */

static inline void
gst_gl_query_start_log(GstGLQuery * query, GstDebugCategory * cat,
    GstDebugLevel level, GObject * object, const gchar * format, ...)
{
  va_list varargs;

  va_start (varargs, format);
  gst_gl_query_start_log_valist (query, cat, level, object, format, varargs);
  va_end (varargs);
}

static inline void
gst_gl_query_counter_log(GstGLQuery * query, GstDebugCategory * cat,
    GstDebugLevel level, GObject * object, const gchar * format, ...)
{
  va_list varargs;

  va_start (varargs, format);
  gst_gl_query_counter_log_valist (query, cat, level, object, format, varargs);
  va_end (varargs);
}

#endif /* G_HAVE_GNUC_VARARGS */
#endif /* G_HAVE_ISO_VARARGS */

G_END_DECLS

#endif /* __GST_GL_QUERY_H__ */
