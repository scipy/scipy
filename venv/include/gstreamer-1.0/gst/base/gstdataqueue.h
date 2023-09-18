/* GStreamer
 * Copyright (C) 2006 Edward Hervey <edward@fluendo.com>
 *
 * gstdataqueue.h:
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


#ifndef __GST_DATA_QUEUE_H__
#define __GST_DATA_QUEUE_H__

#include <gst/gst.h>
#include <gst/base/base-prelude.h>

G_BEGIN_DECLS
#define GST_TYPE_DATA_QUEUE \
  (gst_data_queue_get_type())
#define GST_DATA_QUEUE(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_DATA_QUEUE,GstDataQueue))
#define GST_DATA_QUEUE_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_DATA_QUEUE,GstDataQueueClass))
#define GST_IS_DATA_QUEUE(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_DATA_QUEUE))
#define GST_IS_DATA_QUEUE_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_DATA_QUEUE))

typedef struct _GstDataQueue GstDataQueue;
typedef struct _GstDataQueueClass GstDataQueueClass;
typedef struct _GstDataQueueSize GstDataQueueSize;
typedef struct _GstDataQueueItem GstDataQueueItem;
typedef struct _GstDataQueuePrivate GstDataQueuePrivate;

/**
 * GstDataQueueItem: (skip)
 * @object: the #GstMiniObject to queue.
 * @size: the size in bytes of the miniobject.
 * @duration: the duration in #GstClockTime of the miniobject. Can not be
 * %GST_CLOCK_TIME_NONE.
 * @visible: %TRUE if @object should be considered as a visible object.
 * @destroy: The #GDestroyNotify function to use to free the #GstDataQueueItem.
 * This function should also drop the reference to @object the owner of the
 * #GstDataQueueItem is assumed to hold.
 *
 * Structure used by #GstDataQueue. You can supply a different structure, as
 * long as the top of the structure is identical to this structure.
 */

struct _GstDataQueueItem
{
  GstMiniObject *object;
  guint size;
  guint64 duration;
  gboolean visible;

  /* user supplied destroy function */
  GDestroyNotify destroy;

  /* < private > */
  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstDataQueueSize: (skip)
 * @visible: number of buffers
 * @bytes: number of bytes
 * @time: amount of time
 *
 * Structure describing the size of a queue.
 */
struct _GstDataQueueSize
{
  guint visible;
  guint bytes;
  guint64 time;
};

/**
 * GstDataQueueCheckFullFunction: (skip)
 * @queue: a #GstDataQueue.
 * @visible: The number of visible items currently in the queue.
 * @bytes: The amount of bytes currently in the queue.
 * @time: The accumulated duration of the items currently in the queue.
 * @checkdata: The #gpointer registered when the #GstDataQueue was created.
 *
 * The prototype of the function used to inform the queue that it should be
 * considered as full.
 *
 * Returns: %TRUE if the queue should be considered full.
 */
typedef gboolean (*GstDataQueueCheckFullFunction) (GstDataQueue * queue,
    guint visible, guint bytes, guint64 time, gpointer checkdata);

typedef void (*GstDataQueueFullCallback) (GstDataQueue * queue, gpointer checkdata);
typedef void (*GstDataQueueEmptyCallback) (GstDataQueue * queue, gpointer checkdata);

/**
 * GstDataQueue:
 * @object: the parent structure
 *
 * Opaque #GstDataQueue structure.
 */
struct _GstDataQueue
{
  GObject object;

  /*< private >*/
  GstDataQueuePrivate *priv;
  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstDataQueueClass:
 */
struct _GstDataQueueClass
{
  GObjectClass parent_class;

  /* signals */
  void (*empty) (GstDataQueue * queue);
  void (*full) (GstDataQueue * queue);

  gpointer _gst_reserved[GST_PADDING];
};

GST_BASE_API
GType          gst_data_queue_get_type (void);

GST_BASE_API
GstDataQueue * gst_data_queue_new            (GstDataQueueCheckFullFunction checkfull,
					      GstDataQueueFullCallback fullcallback,
					      GstDataQueueEmptyCallback emptycallback,
					      gpointer checkdata) G_GNUC_MALLOC;
GST_BASE_API
gboolean       gst_data_queue_push           (GstDataQueue * queue, GstDataQueueItem * item);

GST_BASE_API
gboolean       gst_data_queue_push_force     (GstDataQueue * queue, GstDataQueueItem * item);

GST_BASE_API
gboolean       gst_data_queue_pop            (GstDataQueue * queue, GstDataQueueItem ** item);

GST_BASE_API
gboolean       gst_data_queue_peek           (GstDataQueue * queue, GstDataQueueItem ** item);

GST_BASE_API
void           gst_data_queue_flush          (GstDataQueue * queue);

GST_BASE_API
void           gst_data_queue_set_flushing   (GstDataQueue * queue, gboolean flushing);

GST_BASE_API
gboolean       gst_data_queue_drop_head      (GstDataQueue * queue, GType type);

GST_BASE_API
gboolean       gst_data_queue_is_full        (GstDataQueue * queue);

GST_BASE_API
gboolean       gst_data_queue_is_empty       (GstDataQueue * queue);

GST_BASE_API
void           gst_data_queue_get_level      (GstDataQueue * queue, GstDataQueueSize *level);

GST_BASE_API
void           gst_data_queue_limits_changed (GstDataQueue * queue);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstDataQueue, gst_object_unref)

G_END_DECLS

#endif /* __GST_DATA_QUEUE_H__ */
