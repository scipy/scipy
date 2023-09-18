/* GStreamer
 * Copyright (C) 2007 David Schleef <ds@schleef.org>
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

#ifndef _GST_APP_SRC_H_
#define _GST_APP_SRC_H_

#include <gst/gst.h>
#include <gst/base/gstpushsrc.h>
#include <gst/app/app-prelude.h>
#include <gst/app/app-enumtypes.h>

G_BEGIN_DECLS

#define GST_TYPE_APP_SRC \
  (gst_app_src_get_type())
#define GST_APP_SRC(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_APP_SRC,GstAppSrc))
#define GST_APP_SRC_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_APP_SRC,GstAppSrcClass))
#define GST_IS_APP_SRC(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_APP_SRC))
#define GST_IS_APP_SRC_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_APP_SRC))
#define GST_APP_SRC_CAST(obj) \
  ((GstAppSrc*)(obj))

typedef struct _GstAppSrc GstAppSrc;
typedef struct _GstAppSrcClass GstAppSrcClass;
typedef struct _GstAppSrcPrivate GstAppSrcPrivate;

/* FIXME 2.0: Make the instance/class struct private */

/**
 * GstAppSrcCallbacks: (skip)
 * @need_data: Called when the appsrc needs more data. A buffer or EOS should be
 *    pushed to appsrc from this thread or another thread. @length is just a hint
 *    and when it is set to -1, any number of bytes can be pushed into @appsrc.
 * @enough_data: Called when appsrc has enough data. It is recommended that the
 *    application stops calling push-buffer until the need_data callback is
 *    emitted again to avoid excessive buffer queueing.
 * @seek_data: Called when a seek should be performed to the offset.
 *    The next push-buffer should produce buffers from the new @offset.
 *    This callback is only called for seekable stream types.
 *
 * A set of callbacks that can be installed on the appsrc with
 * gst_app_src_set_callbacks().
 */
typedef struct {
  void      (*need_data)    (GstAppSrc *src, guint length, gpointer user_data);
  void      (*enough_data)  (GstAppSrc *src, gpointer user_data);
  gboolean  (*seek_data)    (GstAppSrc *src, guint64 offset, gpointer user_data);

  /*< private >*/
  gpointer     _gst_reserved[GST_PADDING];
} GstAppSrcCallbacks;

/**
 * GstAppStreamType:
 * @GST_APP_STREAM_TYPE_STREAM: No seeking is supported in the stream, such as a
 * live stream.
 * @GST_APP_STREAM_TYPE_SEEKABLE: The stream is seekable but seeking might not
 * be very fast, such as data from a webserver.
 * @GST_APP_STREAM_TYPE_RANDOM_ACCESS: The stream is seekable and seeking is fast,
 * such as in a local file.
 *
 * The stream type.
 */
typedef enum
{
  GST_APP_STREAM_TYPE_STREAM,
  GST_APP_STREAM_TYPE_SEEKABLE,
  GST_APP_STREAM_TYPE_RANDOM_ACCESS
} GstAppStreamType;

/**
 * GstAppLeakyType:
 * @GST_APP_LEAKY_TYPE_NONE: Not Leaky
 * @GST_APP_LEAKY_TYPE_UPSTREAM: Leaky on upstream (new buffers)
 * @GST_APP_LEAKY_TYPE_DOWNSTREAM: Leaky on downstream (old buffers)
 *
 * Buffer dropping scheme to avoid the element's internal queue to block when
 * full.
 *
 * Since: 1.20
 */
typedef enum {
  GST_APP_LEAKY_TYPE_NONE,
  GST_APP_LEAKY_TYPE_UPSTREAM,
  GST_APP_LEAKY_TYPE_DOWNSTREAM
} GstAppLeakyType;

struct _GstAppSrc
{
  GstBaseSrc basesrc;

  /*< private >*/
  GstAppSrcPrivate *priv;

  /*< private >*/
  gpointer     _gst_reserved[GST_PADDING];
};

struct _GstAppSrcClass
{
  GstBaseSrcClass basesrc_class;

  /* signals */
  void          (*need_data)       (GstAppSrc *appsrc, guint length);
  void          (*enough_data)     (GstAppSrc *appsrc);
  gboolean      (*seek_data)       (GstAppSrc *appsrc, guint64 offset);

  /* actions */
  GstFlowReturn (*push_buffer)     (GstAppSrc *appsrc, GstBuffer *buffer);
  GstFlowReturn (*end_of_stream)   (GstAppSrc *appsrc);
  GstFlowReturn (*push_sample)     (GstAppSrc *appsrc, GstSample *sample);
  GstFlowReturn (*push_buffer_list) (GstAppSrc *appsrc, GstBufferList *buffer_list);

  /*< private >*/
  gpointer     _gst_reserved[GST_PADDING-2];
};

GST_APP_API
GType            gst_app_src_get_type                (void);

GST_APP_API
void             gst_app_src_set_caps                (GstAppSrc *appsrc, const GstCaps *caps);

GST_APP_API
GstCaps*         gst_app_src_get_caps                (GstAppSrc *appsrc);

GST_APP_API
void             gst_app_src_set_size                (GstAppSrc *appsrc, gint64 size);

GST_APP_API
gint64           gst_app_src_get_size                (GstAppSrc *appsrc);

GST_APP_API
void             gst_app_src_set_duration            (GstAppSrc *appsrc, GstClockTime duration);

GST_APP_API
GstClockTime     gst_app_src_get_duration            (GstAppSrc *appsrc);

GST_APP_API
void             gst_app_src_set_stream_type         (GstAppSrc *appsrc, GstAppStreamType type);

GST_APP_API
GstAppStreamType gst_app_src_get_stream_type         (GstAppSrc *appsrc);

GST_APP_API
void             gst_app_src_set_max_bytes           (GstAppSrc *appsrc, guint64 max);

GST_APP_API
guint64          gst_app_src_get_max_bytes           (GstAppSrc *appsrc);

GST_APP_API
guint64          gst_app_src_get_current_level_bytes (GstAppSrc *appsrc);

GST_APP_API
void             gst_app_src_set_max_buffers           (GstAppSrc *appsrc, guint64 max);

GST_APP_API
guint64          gst_app_src_get_max_buffers           (GstAppSrc *appsrc);

GST_APP_API
guint64          gst_app_src_get_current_level_buffers (GstAppSrc *appsrc);

GST_APP_API
void             gst_app_src_set_max_time            (GstAppSrc *appsrc, GstClockTime max);

GST_APP_API
GstClockTime     gst_app_src_get_max_time            (GstAppSrc *appsrc);

GST_APP_API
GstClockTime     gst_app_src_get_current_level_time  (GstAppSrc *appsrc);

GST_APP_API
void             gst_app_src_set_leaky_type          (GstAppSrc *appsrc, GstAppLeakyType leaky);

GST_APP_API
GstAppLeakyType  gst_app_src_get_leaky_type          (GstAppSrc *appsrc);

GST_APP_API
void             gst_app_src_set_latency             (GstAppSrc *appsrc, guint64 min, guint64 max);

GST_APP_API
void             gst_app_src_get_latency             (GstAppSrc *appsrc, guint64 *min, guint64 *max);

GST_APP_API
void             gst_app_src_set_emit_signals        (GstAppSrc *appsrc, gboolean emit);

GST_APP_API
gboolean         gst_app_src_get_emit_signals        (GstAppSrc *appsrc);

GST_APP_API
GstFlowReturn    gst_app_src_push_buffer             (GstAppSrc *appsrc, GstBuffer *buffer);

GST_APP_API
GstFlowReturn    gst_app_src_push_buffer_list        (GstAppSrc * appsrc, GstBufferList * buffer_list);

GST_APP_API
GstFlowReturn    gst_app_src_end_of_stream           (GstAppSrc *appsrc);

GST_APP_API
GstFlowReturn    gst_app_src_push_sample             (GstAppSrc *appsrc, GstSample *sample);

GST_APP_API
void             gst_app_src_set_callbacks           (GstAppSrc * appsrc,
                                                      GstAppSrcCallbacks *callbacks,
                                                      gpointer user_data,
                                                      GDestroyNotify notify);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstAppSrc, gst_object_unref)

G_END_DECLS

#endif
