/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2005 Wim Taymans <wim@fluendo.com>
 *
 * gstaudiobasesrc.h:
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

/* a base class for audio sources.
 */

#ifndef __GST_AUDIO_AUDIO_H__
#include <gst/audio/audio.h>
#endif

#ifndef __GST_AUDIO_BASE_SRC_H__
#define __GST_AUDIO_BASE_SRC_H__

#include <gst/gst.h>
#include <gst/base/gstpushsrc.h>

G_BEGIN_DECLS

#define GST_TYPE_AUDIO_BASE_SRC                 (gst_audio_base_src_get_type())
#define GST_AUDIO_BASE_SRC(obj)                 (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_AUDIO_BASE_SRC,GstAudioBaseSrc))
#define GST_AUDIO_BASE_SRC_CAST(obj)            ((GstAudioBaseSrc*)obj)
#define GST_AUDIO_BASE_SRC_CLASS(klass)         (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_AUDIO_BASE_SRC,GstAudioBaseSrcClass))
#define GST_AUDIO_BASE_SRC_GET_CLASS(obj)       (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_AUDIO_BASE_SRC, GstAudioBaseSrcClass))
#define GST_IS_AUDIO_BASE_SRC(obj)              (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_AUDIO_BASE_SRC))
#define GST_IS_AUDIO_BASE_SRC_CLASS(klass)      (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_AUDIO_BASE_SRC))

/**
 * GST_AUDIO_BASE_SRC_CLOCK:
 * @obj: a #GstAudioBaseSrc
 *
 * Get the #GstClock of @obj.
 */
#define GST_AUDIO_BASE_SRC_CLOCK(obj)    (GST_AUDIO_BASE_SRC (obj)->clock)
/**
 * GST_AUDIO_BASE_SRC_PAD:
 * @obj: a #GstAudioBaseSrc
 *
 * Get the source #GstPad of @obj.
 */
#define GST_AUDIO_BASE_SRC_PAD(obj)      (GST_BASE_SRC (obj)->srcpad)

typedef struct _GstAudioBaseSrc GstAudioBaseSrc;
typedef struct _GstAudioBaseSrcClass GstAudioBaseSrcClass;
typedef struct _GstAudioBaseSrcPrivate GstAudioBaseSrcPrivate;

/* FIXME 2.0: Should be "retimestamp" not "re-timestamp" */

/**
 * GstAudioBaseSrcSlaveMethod:
 * @GST_AUDIO_BASE_SRC_SLAVE_RESAMPLE: Resample to match the master clock.
 * @GST_AUDIO_BASE_SRC_SLAVE_RE_TIMESTAMP: Retimestamp output buffers with master
 * clock time.
 * @GST_AUDIO_BASE_SRC_SLAVE_SKEW: Adjust capture pointer when master clock
 * drifts too much.
 * @GST_AUDIO_BASE_SRC_SLAVE_NONE: No adjustment is done.
 *
 * Different possible clock slaving algorithms when the internal audio clock was
 * not selected as the pipeline clock.
 */
typedef enum
{
  GST_AUDIO_BASE_SRC_SLAVE_RESAMPLE,
  GST_AUDIO_BASE_SRC_SLAVE_RE_TIMESTAMP,
  GST_AUDIO_BASE_SRC_SLAVE_SKEW,
  GST_AUDIO_BASE_SRC_SLAVE_NONE
} GstAudioBaseSrcSlaveMethod;

#define GST_AUDIO_BASE_SRC_SLAVE_RETIMESTAMP GST_AUDIO_BASE_SRC_SLAVE_RE_TIMESTAMP

/**
 * GstAudioBaseSrc:
 *
 * Opaque #GstAudioBaseSrc.
 */
struct _GstAudioBaseSrc {
  GstPushSrc          element;

  /*< protected >*/ /* with LOCK */
  /* our ringbuffer */
  GstAudioRingBuffer *ringbuffer;

  /* required buffer and latency */
  GstClockTime        buffer_time;
  GstClockTime        latency_time;

  /* the next sample to write */
  guint64             next_sample;

  /* clock */
  GstClock           *clock;

  /*< private >*/
  GstAudioBaseSrcPrivate *priv;

  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstAudioBaseSrcClass:
 * @parent_class: the parent class.
 * @create_ringbuffer: create and return a #GstAudioRingBuffer to read from.
 *
 * #GstAudioBaseSrc class. Override the vmethod to implement
 * functionality.
 */
struct _GstAudioBaseSrcClass {
  GstPushSrcClass      parent_class;

  /* subclass ringbuffer allocation */
  GstAudioRingBuffer* (*create_ringbuffer)  (GstAudioBaseSrc *src);

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_AUDIO_API
GType      gst_audio_base_src_get_type(void);

GST_AUDIO_API
GstAudioRingBuffer *
           gst_audio_base_src_create_ringbuffer        (GstAudioBaseSrc *src);

GST_AUDIO_API
void       gst_audio_base_src_set_provide_clock        (GstAudioBaseSrc *src, gboolean provide);

GST_AUDIO_API
gboolean   gst_audio_base_src_get_provide_clock        (GstAudioBaseSrc *src);

GST_AUDIO_API
void       gst_audio_base_src_set_slave_method         (GstAudioBaseSrc *src,
                                                        GstAudioBaseSrcSlaveMethod method);
GST_AUDIO_API
GstAudioBaseSrcSlaveMethod
           gst_audio_base_src_get_slave_method         (GstAudioBaseSrc *src);


G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstAudioBaseSrc, gst_object_unref)

G_END_DECLS

#endif /* __GST_AUDIO_BASE_SRC_H__ */
