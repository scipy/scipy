/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2005 Wim Taymans <wim@fluendo.com>
 *
 * gstaudiobasesink.h:
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

/* a base class for audio sinks.
 *
 * It uses a ringbuffer to schedule playback of samples. This makes
 * it very easy to drop or insert samples to align incoming
 * buffers to the exact playback timestamp.
 *
 * Subclasses must provide a ringbuffer pointing to either DMA
 * memory or regular memory. A subclass should also call a callback
 * function when it has played N segments in the buffer. The subclass
 * is free to use a thread to signal this callback, use EIO or any
 * other mechanism.
 *
 * The base class is able to operate in push or pull mode. The chain
 * mode will queue the samples in the ringbuffer as much as possible.
 * The available space is calculated in the callback function.
 *
 * The pull mode will pull_range() a new buffer of N samples with a
 * configurable latency. This allows for high-end real time
 * audio processing pipelines driven by the audiosink. The callback
 * function will be used to perform a pull_range() on the sinkpad.
 * The thread scheduling the callback can be a real-time thread.
 *
 * Subclasses must implement a GstAudioRingBuffer in addition to overriding
 * the methods in GstBaseSink and this class.
 */

#ifndef __GST_AUDIO_AUDIO_H__
#include <gst/audio/audio.h>
#endif

#ifndef __GST_AUDIO_BASE_SINK_H__
#define __GST_AUDIO_BASE_SINK_H__

#include <gst/base/gstbasesink.h>

G_BEGIN_DECLS

#define GST_TYPE_AUDIO_BASE_SINK                (gst_audio_base_sink_get_type())
#define GST_AUDIO_BASE_SINK(obj)                (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_AUDIO_BASE_SINK,GstAudioBaseSink))
#define GST_AUDIO_BASE_SINK_CAST(obj)           ((GstAudioBaseSink*)obj)
#define GST_AUDIO_BASE_SINK_CLASS(klass)        (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_AUDIO_BASE_SINK,GstAudioBaseSinkClass))
#define GST_AUDIO_BASE_SINK_GET_CLASS(obj)      (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_AUDIO_BASE_SINK, GstAudioBaseSinkClass))
#define GST_IS_AUDIO_BASE_SINK(obj)             (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_AUDIO_BASE_SINK))
#define GST_IS_AUDIO_BASE_SINK_CLASS(klass)     (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_AUDIO_BASE_SINK))

/**
 * GST_AUDIO_BASE_SINK_CLOCK:
 * @obj: a #GstAudioBaseSink
 *
 * Get the #GstClock of @obj.
 */
#define GST_AUDIO_BASE_SINK_CLOCK(obj)   (GST_AUDIO_BASE_SINK (obj)->clock)
/**
 * GST_AUDIO_BASE_SINK_PAD:
 * @obj: a #GstAudioBaseSink
 *
 * Get the sink #GstPad of @obj.
 */
#define GST_AUDIO_BASE_SINK_PAD(obj)     (GST_BASE_SINK (obj)->sinkpad)

/**
 * GstAudioBaseSinkSlaveMethod:
 * @GST_AUDIO_BASE_SINK_SLAVE_RESAMPLE: Resample to match the master clock
 * @GST_AUDIO_BASE_SINK_SLAVE_SKEW: Adjust playout pointer when master clock
 * drifts too much.
 * @GST_AUDIO_BASE_SINK_SLAVE_NONE: No adjustment is done.
 * @GST_AUDIO_BASE_SINK_SLAVE_CUSTOM: Use custom clock slaving algorithm (Since: 1.6)
 *
 * Different possible clock slaving algorithms used when the internal audio
 * clock is not selected as the pipeline master clock.
 */
typedef enum
{
  GST_AUDIO_BASE_SINK_SLAVE_RESAMPLE,
  GST_AUDIO_BASE_SINK_SLAVE_SKEW,
  GST_AUDIO_BASE_SINK_SLAVE_NONE,
  GST_AUDIO_BASE_SINK_SLAVE_CUSTOM
} GstAudioBaseSinkSlaveMethod;

typedef struct _GstAudioBaseSink GstAudioBaseSink;
typedef struct _GstAudioBaseSinkClass GstAudioBaseSinkClass;
typedef struct _GstAudioBaseSinkPrivate GstAudioBaseSinkPrivate;

/**
 * GstAudioBaseSinkDiscontReason:
 * @GST_AUDIO_BASE_SINK_DISCONT_REASON_NO_DISCONT: No discontinuity occurred
 * @GST_AUDIO_BASE_SINK_DISCONT_REASON_NEW_CAPS: New caps are set, causing renegotiotion
 * @GST_AUDIO_BASE_SINK_DISCONT_REASON_FLUSH: Samples have been flushed
 * @GST_AUDIO_BASE_SINK_DISCONT_REASON_SYNC_LATENCY: Sink was synchronized to the estimated latency (occurs during initialization)
 * @GST_AUDIO_BASE_SINK_DISCONT_REASON_ALIGNMENT: Aligning buffers failed because the timestamps are too discontinuous
 * @GST_AUDIO_BASE_SINK_DISCONT_REASON_DEVICE_FAILURE: Audio output device experienced and recovered from an error but introduced latency in the process (see also gst_audio_base_sink_report_device_failure())
 *
 * Different possible reasons for discontinuities. This enum is useful for the custom
 * slave method.
 *
 * Since: 1.6
 */
typedef enum
{
  GST_AUDIO_BASE_SINK_DISCONT_REASON_NO_DISCONT,
  GST_AUDIO_BASE_SINK_DISCONT_REASON_NEW_CAPS,
  GST_AUDIO_BASE_SINK_DISCONT_REASON_FLUSH,
  GST_AUDIO_BASE_SINK_DISCONT_REASON_SYNC_LATENCY,
  GST_AUDIO_BASE_SINK_DISCONT_REASON_ALIGNMENT,
  GST_AUDIO_BASE_SINK_DISCONT_REASON_DEVICE_FAILURE
} GstAudioBaseSinkDiscontReason;

/**
 * GstAudioBaseSinkCustomSlavingCallback:
 * @sink: a #GstAudioBaseSink
 * @etime: external clock time
 * @itime: internal clock time
 * @requested_skew: skew amount requested by the callback
 * @discont_reason: reason for discontinuity (if any)
 * @user_data: user data
 *
 * This function is set with gst_audio_base_sink_set_custom_slaving_callback()
 * and is called during playback. It receives the current time of external and
 * internal clocks, which the callback can then use to apply any custom
 * slaving/synchronization schemes.
 *
 * The external clock is the sink's element clock, the internal one is the
 * internal audio clock. The internal audio clock's calibration is applied to
 * the timestamps before they are passed to the callback. The difference between
 * etime and itime is the skew; how much internal and external clock lie apart
 * from each other. A skew of 0 means both clocks are perfectly in sync.
 * itime > etime means the external clock is going slower, while itime < etime
 * means it is going faster than the internal clock. etime and itime are always
 * valid timestamps, except for when a discontinuity happens.
 *
 * requested_skew is an output value the callback can write to. It informs the
 * sink of whether or not it should move the playout pointer, and if so, by how
 * much. This pointer is only NULL if a discontinuity occurs; otherwise, it is
 * safe to write to *requested_skew. The default skew is 0.
 *
 * The sink may experience discontinuities. If one happens, discont is TRUE,
 * itime, etime are set to GST_CLOCK_TIME_NONE, and requested_skew is NULL.
 * This makes it possible to reset custom clock slaving algorithms when a
 * discontinuity happens.
 *
 * Since: 1.6
 */
typedef void (*GstAudioBaseSinkCustomSlavingCallback) (GstAudioBaseSink *sink, GstClockTime etime, GstClockTime itime, GstClockTimeDiff *requested_skew, GstAudioBaseSinkDiscontReason discont_reason, gpointer user_data);

/**
 * GstAudioBaseSink:
 *
 * Opaque #GstAudioBaseSink.
 */
struct _GstAudioBaseSink {
  GstBaseSink         element;

  /*< protected >*/ /* with LOCK */
  /* our ringbuffer */
  GstAudioRingBuffer *ringbuffer;

  /* required buffer and latency in microseconds */
  guint64             buffer_time;
  guint64             latency_time;

  /* the next sample to write */
  guint64             next_sample;

  /* clock */
  GstClock           *provided_clock;

  /* with g_atomic_; currently rendering eos */
  gboolean            eos_rendering;

  /*< private >*/
  GstAudioBaseSinkPrivate *priv;

  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstAudioBaseSinkClass:
 * @parent_class: the parent class.
 * @create_ringbuffer: create and return a #GstAudioRingBuffer to write to.
 * @payload: payload data in a format suitable to write to the sink. If no
 *           payloading is required, returns a reffed copy of the original
 *           buffer, else returns the payloaded buffer with all other metadata
 *           copied.
 *
 * #GstAudioBaseSink class. Override the vmethod to implement
 * functionality.
 */
struct _GstAudioBaseSinkClass {
  GstBaseSinkClass     parent_class;

  /* subclass ringbuffer allocation */
  GstAudioRingBuffer* (*create_ringbuffer)  (GstAudioBaseSink *sink);

  /* subclass payloader */
  GstBuffer*          (*payload)            (GstAudioBaseSink *sink,
                                             GstBuffer        *buffer);
  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_AUDIO_API
GType gst_audio_base_sink_get_type(void);

GST_AUDIO_API
GstAudioRingBuffer *
           gst_audio_base_sink_create_ringbuffer       (GstAudioBaseSink *sink);

GST_AUDIO_API
void       gst_audio_base_sink_set_provide_clock       (GstAudioBaseSink *sink, gboolean provide);

GST_AUDIO_API
gboolean   gst_audio_base_sink_get_provide_clock       (GstAudioBaseSink *sink);

GST_AUDIO_API
void       gst_audio_base_sink_set_slave_method        (GstAudioBaseSink *sink,
                                                        GstAudioBaseSinkSlaveMethod method);
GST_AUDIO_API
GstAudioBaseSinkSlaveMethod
           gst_audio_base_sink_get_slave_method        (GstAudioBaseSink *sink);

GST_AUDIO_API
void       gst_audio_base_sink_set_drift_tolerance     (GstAudioBaseSink *sink,
                                                        gint64 drift_tolerance);
GST_AUDIO_API
gint64     gst_audio_base_sink_get_drift_tolerance     (GstAudioBaseSink *sink);

GST_AUDIO_API
void       gst_audio_base_sink_set_alignment_threshold (GstAudioBaseSink * sink,
                                                        GstClockTime alignment_threshold);
GST_AUDIO_API
GstClockTime
           gst_audio_base_sink_get_alignment_threshold (GstAudioBaseSink * sink);

GST_AUDIO_API
void       gst_audio_base_sink_set_discont_wait        (GstAudioBaseSink * sink,
                                                        GstClockTime discont_wait);
GST_AUDIO_API
GstClockTime
           gst_audio_base_sink_get_discont_wait        (GstAudioBaseSink * sink);

GST_AUDIO_API
void
gst_audio_base_sink_set_custom_slaving_callback        (GstAudioBaseSink * sink,
                                                        GstAudioBaseSinkCustomSlavingCallback callback,
                                                        gpointer user_data,
                                                        GDestroyNotify notify);

GST_AUDIO_API
void gst_audio_base_sink_report_device_failure         (GstAudioBaseSink * sink);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstAudioBaseSink, gst_object_unref)

G_END_DECLS

#endif /* __GST_AUDIO_BASE_SINK_H__ */
