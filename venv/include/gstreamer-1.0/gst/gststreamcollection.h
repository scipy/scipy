/* GStreamer
 * Copyright (C) 2015 Centricular Ltd
 *  @author: Edward Hervey <edward@centricular.com>
 *  @author: Jan Schmidt <jan@centricular.com>
 *
 * gststreams.h : Header for GstStreamCollection subsystem
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


#ifndef __GST_STREAM_COLLECTION_H__
#define __GST_STREAM_COLLECTION_H__

#include <gst/gstobject.h>

G_BEGIN_DECLS

#define GST_TYPE_STREAM_COLLECTION             (gst_stream_collection_get_type ())
#define GST_IS_STREAM_COLLECTION(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_STREAM_COLLECTION))
#define GST_IS_STREAM_COLLECTION_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), GST_TYPE_STREAM_COLLECTION))
#define GST_STREAM_COLLECTION_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_STREAM_COLLECTION, GstStreamCollectionClass))
#define GST_STREAM_COLLECTION(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_STREAM_COLLECTION, GstStreamCollection))
#define GST_STREAM_COLLECTION_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), GST_TYPE_STREAM_COLLECTION, GstStreamCollectionClass))
#define GST_STREAM_COLLECTION_CAST(obj)        ((GstStreamCollection*)(obj))

typedef struct _GstStreamCollection GstStreamCollection;
typedef struct _GstStreamCollectionClass GstStreamCollectionClass;
typedef struct _GstStreamCollectionPrivate GstStreamCollectionPrivate;

#include <gst/gststreamcollection.h>
#include <gst/gststreams.h>

/**
 * GstStreamCollection:
 *
 * A collection of #GstStream that are available.
 *
 * A #GstStreamCollection will be provided by elements that can make those
 * streams available. Applications can use the collection to show the user
 * what streams are available by using %gst_stream_collection_get_stream()
 *
 * Once posted, a #GstStreamCollection is immutable. Updates are made by sending
 * a new #GstStreamCollection message, which may or may not share some of
 * the #GstStream objects from the collection it replaces. The receiver can check
 * the sender of a stream collection message to know which collection is
 * obsoleted.
 *
 * Several elements in a pipeline can provide #GstStreamCollection.
 *
 * Applications can activate streams from a collection by using the
 * #GST_EVENT_SELECT_STREAMS event on a pipeline, bin or element.
 *
 * Since: 1.10
 */
struct _GstStreamCollection {
  /*< private >*/
  GstObject object;

  gchar *upstream_id;
  GstStreamCollectionPrivate *priv;

  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstStreamCollectionClass:
 * @parent_class: the parent class structure
 * @stream_notify: default signal handler for the stream-notify signal
 *
 * GstStreamCollection class structure
 */
struct _GstStreamCollectionClass {
  GstObjectClass parent_class;

  /* signals */
  void  (*stream_notify)      (GstStreamCollection *collection, GstStream *stream, GParamSpec * pspec);

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_API
GType gst_stream_collection_get_type (void);

GST_API
GstStreamCollection *gst_stream_collection_new (const gchar *upstream_id);

GST_API
const gchar *gst_stream_collection_get_upstream_id (GstStreamCollection *collection);

GST_API
guint gst_stream_collection_get_size (GstStreamCollection *collection);

GST_API
GstStream *gst_stream_collection_get_stream (GstStreamCollection *collection, guint index);

GST_API
gboolean gst_stream_collection_add_stream (GstStreamCollection *collection,
                                           GstStream *stream);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstStreamCollection, gst_object_unref)

G_END_DECLS

#endif /* __GST_STREAM_COLLECTION_H__ */
