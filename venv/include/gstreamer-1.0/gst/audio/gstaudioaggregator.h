/* GStreamer
 * Copyright (C) 2014 Collabora
 *   Author: Olivier Crete <olivier.crete@collabora.com>
 *
 * gstaudioaggregator.h:
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

#ifndef __GST_AUDIO_AGGREGATOR_H__
#define __GST_AUDIO_AGGREGATOR_H__

#include <gst/gst.h>
#include <gst/base/gstaggregator.h>
#include <gst/audio/audio.h>

G_BEGIN_DECLS

/*******************************
 * GstAudioAggregator Structs  *
 *******************************/

typedef struct _GstAudioAggregator GstAudioAggregator;
typedef struct _GstAudioAggregatorPrivate GstAudioAggregatorPrivate;
typedef struct _GstAudioAggregatorClass GstAudioAggregatorClass;


/************************
 * GstAudioAggregatorPad API *
 ***********************/

#define GST_TYPE_AUDIO_AGGREGATOR_PAD            (gst_audio_aggregator_pad_get_type())
#define GST_AUDIO_AGGREGATOR_PAD(obj)            (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_AUDIO_AGGREGATOR_PAD, GstAudioAggregatorPad))
#define GST_AUDIO_AGGREGATOR_PAD_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_AUDIO_AGGREGATOR_PAD, GstAudioAggregatorPadClass))
#define GST_AUDIO_AGGREGATOR_PAD_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS ((obj),GST_TYPE_AUDIO_AGGREGATOR_PAD, GstAudioAggregatorPadClass))
#define GST_IS_AUDIO_AGGREGATOR_PAD(obj)         (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_AUDIO_AGGREGATOR_PAD))
#define GST_IS_AUDIO_AGGREGATOR_PAD_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_AUDIO_AGGREGATOR_PAD))

/****************************
 * GstAudioAggregatorPad Structs *
 ***************************/

typedef struct _GstAudioAggregatorPad GstAudioAggregatorPad;
typedef struct _GstAudioAggregatorPadClass GstAudioAggregatorPadClass;
typedef struct _GstAudioAggregatorPadPrivate GstAudioAggregatorPadPrivate;

/**
 * GstAudioAggregatorPad:
 * @info: The audio info for this pad set from the incoming caps
 *
 * The default implementation of GstPad used with #GstAudioAggregator
 *
 * Since: 1.14
 */
struct _GstAudioAggregatorPad
{
  GstAggregatorPad                  parent;

  /*< public >*/
  /* read-only, with OBJECT_LOCK */
  GstAudioInfo                      info;

  /*< private >*/
  GstAudioAggregatorPadPrivate     *priv;

  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstAudioAggregatorPadClass:
 * @convert_buffer: Convert a buffer from one format to another.
 * @update_conversion_info: Called when either the input or output
 *  formats have changed.
 *
 * Since: 1.14
 */
struct _GstAudioAggregatorPadClass
  {
  GstAggregatorPadClass   parent_class;

  GstBuffer * (* convert_buffer) (GstAudioAggregatorPad * pad,
                                  GstAudioInfo *in_info,
                                  GstAudioInfo *out_info,
                                  GstBuffer * buffer);

  void        (* update_conversion_info) (GstAudioAggregatorPad *pad);

  /*< private >*/
  gpointer      _gst_reserved[GST_PADDING_LARGE];
};

GST_AUDIO_API
GType gst_audio_aggregator_pad_get_type           (void);

#define GST_TYPE_AUDIO_AGGREGATOR_CONVERT_PAD            (gst_audio_aggregator_convert_pad_get_type())
#define GST_AUDIO_AGGREGATOR_CONVERT_PAD(obj)            (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_AUDIO_AGGREGATOR_CONVERT_PAD, GstAudioAggregatorConvertPad))
#define GST_AUDIO_AGGREGATOR_CONVERT_PAD_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_AUDIO_AGGREGATOR_CONVERT_PAD, GstAudioAggregatorConvertPadClass))
#define GST_AUDIO_AGGREGATOR_CONVERT_PAD_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS ((obj),GST_TYPE_AUDIO_AGGREGATOR_CONVERT_PAD, GstAudioAggregatorConvertPadClass))
#define GST_IS_AUDIO_AGGREGATOR_CONVERT_PAD(obj)         (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_AUDIO_AGGREGATOR_CONVERT_PAD))
#define GST_IS_AUDIO_AGGREGATOR_CONVERT_PAD_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_AUDIO_AGGREGATOR_CONVERT_PAD))

/****************************
 * GstAudioAggregatorPad Structs *
 ***************************/

typedef struct _GstAudioAggregatorConvertPad GstAudioAggregatorConvertPad;
typedef struct _GstAudioAggregatorConvertPadClass GstAudioAggregatorConvertPadClass;
typedef struct _GstAudioAggregatorConvertPadPrivate GstAudioAggregatorConvertPadPrivate;

/**
 * GstAudioAggregatorConvertPad:
 *
 * An implementation of GstPad that can be used with #GstAudioAggregator.
 *
 * See #GstAudioAggregator for more details.
 *
 * Since: 1.14
 */
struct _GstAudioAggregatorConvertPad
{
  /*< private >*/
  GstAudioAggregatorPad                  parent;

  GstAudioAggregatorConvertPadPrivate   *priv;

  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstAudioAggregatorConvertPadClass:
 *
 * Since: 1.14
 */
struct _GstAudioAggregatorConvertPadClass
{
  GstAudioAggregatorPadClass   parent_class;

  /*< private >*/
  gpointer      _gst_reserved[GST_PADDING];
};

GST_AUDIO_API
GType gst_audio_aggregator_convert_pad_get_type           (void);

/**************************
 * GstAudioAggregator API *
 **************************/

#define GST_TYPE_AUDIO_AGGREGATOR            (gst_audio_aggregator_get_type())
#define GST_AUDIO_AGGREGATOR(obj)            (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_AUDIO_AGGREGATOR,GstAudioAggregator))
#define GST_AUDIO_AGGREGATOR_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_AUDIO_AGGREGATOR,GstAudioAggregatorClass))
#define GST_AUDIO_AGGREGATOR_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS ((obj),GST_TYPE_AUDIO_AGGREGATOR,GstAudioAggregatorClass))
#define GST_IS_AUDIO_AGGREGATOR(obj)         (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_AUDIO_AGGREGATOR))
#define GST_IS_AUDIO_AGGREGATOR_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_AUDIO_AGGREGATOR))

/**
 * GstAudioAggregator:
 * @current_caps: The caps set by the subclass
 *
 * GstAudioAggregator object
 *
 * Since: 1.14
 */
struct _GstAudioAggregator
{
  GstAggregator              parent;

  /*< public >*/
  GstCaps                   *current_caps;

  /*< private >*/
  GstAudioAggregatorPrivate *priv;

  gpointer                  _gst_reserved[GST_PADDING];
};

/**
 * GstAudioAggregatorClass:
 * @create_output_buffer: Create a new output buffer contains num_frames frames.
 * @aggregate_one_buffer: Aggregates one input buffer to the output
 *  buffer.  The in_offset and out_offset are in "frames", which is
 *  the size of a sample times the number of channels. Returns TRUE if
 *  any non-silence was added to the buffer
 *
 * Since: 1.14
 */
struct _GstAudioAggregatorClass {
  GstAggregatorClass   parent_class;

  /*< public >*/
  GstBuffer * (* create_output_buffer) (GstAudioAggregator * aagg,
      guint num_frames);
  gboolean (* aggregate_one_buffer) (GstAudioAggregator * aagg,
      GstAudioAggregatorPad * pad, GstBuffer * inbuf, guint in_offset,
      GstBuffer * outbuf, guint out_offset, guint num_frames);

  /*< private >*/
  gpointer          _gst_reserved[GST_PADDING_LARGE];
};

/*************************
 * GstAggregator methods *
 ************************/

GST_AUDIO_API
GType gst_audio_aggregator_get_type(void);

GST_AUDIO_API
void  gst_audio_aggregator_set_sink_caps (GstAudioAggregator    * aagg,
                                          GstAudioAggregatorPad * pad,
                                          GstCaps               * caps);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstAudioAggregator, gst_object_unref)
G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstAudioAggregatorPad, gst_object_unref)
G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstAudioAggregatorConvertPad, gst_object_unref)

G_END_DECLS

#endif /* __GST_AUDIO_AGGREGATOR_H__ */
