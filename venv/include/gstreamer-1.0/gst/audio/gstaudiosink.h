/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2005 Wim Taymans <wim@fluendo.com>
 *
 * gstaudiosink.h:
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

#ifndef __GST_AUDIO_AUDIO_H__
#include <gst/audio/audio.h>
#endif

#ifndef __GST_AUDIO_SINK_H__
#define __GST_AUDIO_SINK_H__

#include <gst/gst.h>
#include <gst/audio/gstaudiobasesink.h>

G_BEGIN_DECLS

#define GST_TYPE_AUDIO_SINK             (gst_audio_sink_get_type())
#define GST_AUDIO_SINK(obj)             (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_AUDIO_SINK,GstAudioSink))
#define GST_AUDIO_SINK_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_AUDIO_SINK,GstAudioSinkClass))
#define GST_AUDIO_SINK_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj),GST_TYPE_AUDIO_SINK,GstAudioSinkClass))
#define GST_IS_AUDIO_SINK(obj)          (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_AUDIO_SINK))
#define GST_IS_AUDIO_SINK_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_AUDIO_SINK))

typedef struct _GstAudioSink GstAudioSink;
typedef struct _GstAudioSinkClass GstAudioSinkClass;
typedef struct _GstAudioSinkClassExtension GstAudioSinkClassExtension;

/**
 * GstAudioSink:
 *
 * Opaque #GstAudioSink.
 */
struct _GstAudioSink {
  GstAudioBaseSink       element;

  /*< private >*/ /* with LOCK */
  GThread   *thread;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstAudioSinkClass:
 * @parent_class: the parent class structure.
 * @open: Open the device. No configuration needs to be done at this point.
 *        This function is also used to check if the device is available.
 * @prepare: Prepare the device to operate with the specified parameters.
 * @unprepare: Undo operations done in prepare.
 * @close: Close the device.
 * @write: Write data to the device.
 *         This vmethod is allowed to block until all the data is written.
 *         If such is the case then it is expected that pause, stop and
 *         reset will unblock the write when called.
 * @delay: Return how many frames are still in the device. Participates in
 *         computing the time for audio clocks and drives the synchronisation.
 * @reset: Returns as quickly as possible from a write and flush any pending
 *         samples from the device.
 *         This vmethod is deprecated. Please provide pause and stop instead.
 * @pause: Pause the device and unblock write as fast as possible.
 *         For retro compatibility, the audio sink will fallback
 *         to calling reset if this vmethod is not provided. Since: 1.18
 * @resume: Resume the device. Since: 1.18
 * @stop: Stop the device and unblock write as fast as possible.
 *        Pending samples are flushed from the device.
 *        For retro compatibility, the audio sink will fallback
 *        to calling reset if this vmethod is not provided. Since: 1.18
 * @extension: class extension structure. Since: 1.18
 */
struct _GstAudioSinkClass {
  GstAudioBaseSinkClass parent_class;

  /* vtable */

  /* open the device with given specs */
  gboolean (*open)      (GstAudioSink *sink);
  /* prepare resources and state to operate with the given specs */
  gboolean (*prepare)   (GstAudioSink *sink, GstAudioRingBufferSpec *spec);
  /* undo anything that was done in prepare() */
  gboolean (*unprepare) (GstAudioSink *sink);
  /* close the device */
  gboolean (*close)     (GstAudioSink *sink);
  /**
   * GstAudioSinkClass::write:
   * @data: (type guint8) (array length=length): the sample data
   *
   * Write samples to the device.
   */
  gint     (*write)     (GstAudioSink *sink, gpointer data, guint length);
  /* get number of frames queued in the device */
  guint    (*delay)     (GstAudioSink *sink);
  /* deprecated: reset the audio device, unblock from a write */
  void     (*reset)     (GstAudioSink *sink);
  /* pause the audio device, unblock from a write */
  void     (*pause)     (GstAudioSink *sink);
  /* resume the audio device */
  void     (*resume)    (GstAudioSink *sink);
  /* stop the audio device, unblock from a write */
  void     (*stop)      (GstAudioSink *sink);

  GstAudioSinkClassExtension *extension;
};

/**
 * GstAudioSinkClassExtension:
 * @clear-all: Clear the device. Since: 1.18
 */
struct _GstAudioSinkClassExtension
{
  /* clear the audio device */
  void     (*clear_all) (GstAudioSink *sink);

  /* no padding needed  */
};

GST_AUDIO_API
GType gst_audio_sink_get_type(void);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstAudioSink, gst_object_unref)

G_END_DECLS

#endif /* __GST_AUDIO_SINK_H__ */
