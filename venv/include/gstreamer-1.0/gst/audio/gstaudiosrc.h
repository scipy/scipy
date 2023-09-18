/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2005 Wim Taymans <wim@fluendo.com>
 *
 * gstaudiosrc.h:
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

#ifndef __GST_AUDIO_SRC_H__
#define __GST_AUDIO_SRC_H__

#include <gst/gst.h>
#include <gst/audio/gstaudiobasesrc.h>

G_BEGIN_DECLS

#define GST_TYPE_AUDIO_SRC              (gst_audio_src_get_type())
#define GST_AUDIO_SRC(obj)              (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_AUDIO_SRC,GstAudioSrc))
#define GST_AUDIO_SRC_CLASS(klass)      (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_AUDIO_SRC,GstAudioSrcClass))
#define GST_AUDIO_SRC_GET_CLASS(obj)    (G_TYPE_INSTANCE_GET_CLASS ((obj),GST_TYPE_AUDIO_SRC,GstAudioSrcClass))
#define GST_IS_AUDIO_SRC(obj)           (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_AUDIO_SRC))
#define GST_IS_AUDIO_SRC_CLASS(klass)   (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_AUDIO_SRC))

typedef struct _GstAudioSrc GstAudioSrc;
typedef struct _GstAudioSrcClass GstAudioSrcClass;

/**
 * GstAudioSrc:
 *
 * Base class for simple audio sources.
 */
struct _GstAudioSrc {
  GstAudioBaseSrc        element;

  /*< private >*/ /* with LOCK */
  GThread   *thread;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstAudioSrcClass:
 * @parent_class: the parent class.
 * @open: open the device with the specified caps
 * @prepare: configure device with format
 * @unprepare: undo the configuration
 * @close: close the device
 * @read: read samples from the audio device
 * @delay: the number of frames queued in the device
 * @reset: unblock a read to the device and reset.
 *
 * #GstAudioSrc class. Override the vmethod to implement
 * functionality.
 */
struct _GstAudioSrcClass {
  GstAudioBaseSrcClass parent_class;

  /* vtable */

  /* open the device with given specs */
  gboolean (*open)      (GstAudioSrc *src);
  /* prepare resources and state to operate with the given specs */
  gboolean (*prepare)   (GstAudioSrc *src, GstAudioRingBufferSpec *spec);
  /* undo anything that was done in prepare() */
  gboolean (*unprepare) (GstAudioSrc *src);
  /* close the device */
  gboolean (*close)     (GstAudioSrc *src);
  /**
   * GstAudioSrcClass::read:
   * @data: (type guint8) (array length=length): the sample data
   * @timestamp: (out): a #GstClockTime
   *
   * Read samples from the device.
   */
  guint    (*read)      (GstAudioSrc *src, gpointer data, guint length,
      GstClockTime *timestamp);
  /* get number of frames queued in the device */
  guint    (*delay)     (GstAudioSrc *src);
  /* reset the audio device, unblock from a write */
  void     (*reset)     (GstAudioSrc *src);

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_AUDIO_API
GType gst_audio_src_get_type(void);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstAudioSrc, gst_object_unref)

G_END_DECLS

#endif /* __GST_AUDIO_SRC_H__ */
