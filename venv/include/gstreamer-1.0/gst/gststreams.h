/* GStreamer
 * Copyright (C) 2015 Centricular Ltd
 *  @author: Edward Hervey <edward@centricular.com>
 *  @author: Jan Schmidt <jan@centricular.com>
 *
 * gststreams.h : Header for GstStream subsystem
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


#ifndef __GST_STREAMS_H__
#define __GST_STREAMS_H__

#include <gst/gstobject.h>

G_BEGIN_DECLS

#define GST_TYPE_STREAM             (gst_stream_get_type ())
#define GST_IS_STREAM(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_STREAM))
#define GST_IS_STREAM_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), GST_TYPE_STREAM))
#define GST_STREAM_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_STREAM, GstStreamClass))
#define GST_STREAM(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_STREAM, GstStream))
#define GST_STREAM_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), GST_TYPE_STREAM, GstStreamClass))
#define GST_STREAM_CAST(obj)        ((GstStream*)(obj))

/**
 * GstStreamType:
 * @GST_STREAM_TYPE_UNKNOWN: The stream is of unknown (unclassified) type.
 * @GST_STREAM_TYPE_AUDIO: The stream is of audio data
 * @GST_STREAM_TYPE_VIDEO: The stream carries video data
 * @GST_STREAM_TYPE_CONTAINER: The stream is a muxed container type
 * @GST_STREAM_TYPE_TEXT: The stream contains subtitle / subpicture data.
 *
 * #GstStreamType describes a high level classification set for
 * flows of data in #GstStream objects.
 *
 * Note that this is a flag, and therefore users should not assume it
 * will be a single value. Do not use the equality operator for checking
 * whether a stream is of a certain type.
 *
 * Since: 1.10
 */
typedef enum {
  GST_STREAM_TYPE_UNKNOWN   = 1 << 0,
  GST_STREAM_TYPE_AUDIO     = 1 << 1,
  GST_STREAM_TYPE_VIDEO     = 1 << 2,
  GST_STREAM_TYPE_CONTAINER = 1 << 3,
  GST_STREAM_TYPE_TEXT      = 1 << 4
} GstStreamType;


typedef struct _GstStream GstStream;
typedef struct _GstStreamClass GstStreamClass;
typedef struct _GstStreamPrivate GstStreamPrivate;

/**
 * GstStream:
 * @stream_id: The Stream Identifier for this #GstStream
 *
 * A high-level object representing a single stream. It might be backed, or
 * not, by an actual flow of data in a pipeline (#GstPad).
 *
 * A #GstStream does not care about data changes (such as decoding, encoding,
 * parsing,...) as long as the underlying data flow corresponds to the same
 * high-level flow (ex: a certain audio track).
 *
 * A #GstStream contains all the information pertinent to a stream, such as
 * stream-id, tags, caps, type, ...
 *
 * Elements can subclass a #GstStream for internal usage (to contain information
 * pertinent to streams of data).
 *
 * Since: 1.10
 */
struct _GstStream {
  /*< private >*/
  GstObject object;

  /*< public >*/
  const gchar *stream_id;

  /*< private >*/
  GstStreamPrivate *priv;

  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstStreamClass:
 * @parent_class: the parent class structure
 *
 * GstStream class structure
 */
struct _GstStreamClass {
  GstObjectClass parent_class;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_API
GType     gst_stream_get_type (void);

#include <gst/gstevent.h>

GST_API
GstStream *gst_stream_new            (const gchar *stream_id,
				      GstCaps *caps,
				      GstStreamType type,
				      GstStreamFlags flags);
GST_API
const gchar *  gst_stream_get_stream_id (GstStream *stream);

GST_API
void           gst_stream_set_stream_flags (GstStream *stream, GstStreamFlags flags);

GST_API
GstStreamFlags gst_stream_get_stream_flags (GstStream *stream);

GST_API
void           gst_stream_set_stream_type (GstStream *stream, GstStreamType stream_type);

GST_API
GstStreamType  gst_stream_get_stream_type (GstStream *stream);

GST_API
void           gst_stream_set_tags (GstStream *stream, GstTagList *tags);

GST_API
GstTagList *   gst_stream_get_tags (GstStream *stream);

GST_API
void           gst_stream_set_caps (GstStream *stream, GstCaps *caps);

GST_API
GstCaps *      gst_stream_get_caps (GstStream *stream);

GST_API
const gchar *  gst_stream_type_get_name (GstStreamType stype);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstStream, gst_object_unref)

G_END_DECLS

#endif /* __GST_STREAMS_H__ */
