/* GStreamer
 * Copyright (C) <1999> Erik Walthinsen <omega@cse.ogi.edu>
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

#ifndef __GST_AUDIO_FILTER_H__
#define __GST_AUDIO_FILTER_H__

#include <gst/gst.h>
#include <gst/base/gstbasetransform.h>

G_BEGIN_DECLS

typedef struct _GstAudioFilter GstAudioFilter;
typedef struct _GstAudioFilterClass GstAudioFilterClass;

#define GST_TYPE_AUDIO_FILTER \
  (gst_audio_filter_get_type())
#define GST_AUDIO_FILTER(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_AUDIO_FILTER,GstAudioFilter))
#define GST_AUDIO_FILTER_CAST(obj) \
  ((GstAudioFilter *) (obj))
#define GST_AUDIO_FILTER_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_AUDIO_FILTER,GstAudioFilterClass))
#define GST_AUDIO_FILTER_CLASS_CAST(klass) \
  ((GstAudioFilterClass *) (klass))
#define GST_AUDIO_FILTER_GET_CLASS(obj) \
  (G_TYPE_INSTANCE_GET_CLASS((obj),GST_TYPE_AUDIO_FILTER,GstAudioFilterClass))
#define GST_IS_AUDIO_FILTER(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_AUDIO_FILTER))
#define GST_IS_AUDIO_FILTER_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_AUDIO_FILTER))

/**
 * GstAudioFilter:
 *
 * Base class for audio filters with the same format for input and output.
 */
struct _GstAudioFilter {
  GstBaseTransform basetransform;

  /*< protected >*/
  GstAudioInfo info;   /* currently configured format */

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

#define GST_AUDIO_FILTER_INFO(filter)     (&GST_AUDIO_FILTER_CAST(filter)->info)

#define GST_AUDIO_FILTER_FORMAT(filter)   (GST_AUDIO_INFO_FORMAT(GST_AUDIO_FILTER_INFO(filter)))
#define GST_AUDIO_FILTER_RATE(filter)     (GST_AUDIO_INFO_RATE(GST_AUDIO_FILTER_INFO(filter)))
#define GST_AUDIO_FILTER_CHANNELS(filter) (GST_AUDIO_INFO_CHANNELS(GST_AUDIO_FILTER_INFO(filter)))
#define GST_AUDIO_FILTER_BPF(filter)      (GST_AUDIO_INFO_BPF(GST_AUDIO_FILTER_INFO(filter)))
#define GST_AUDIO_FILTER_BPS(filter)      (GST_AUDIO_INFO_BPS(GST_AUDIO_FILTER_INFO(filter)))

/**
 * GstAudioFilterClass:
 * @basetransformclass: parent class
 * @setup: virtual function called whenever the format changes
 *
 * In addition to the @setup virtual function, you should also override the
 * GstBaseTransform::transform and/or GstBaseTransform::transform_ip virtual
 * function.
 */

struct _GstAudioFilterClass {
  GstBaseTransformClass basetransformclass;

  /* virtual function, called whenever the format changes */
  gboolean  (*setup) (GstAudioFilter * filter, const GstAudioInfo * info);

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_AUDIO_API
GType   gst_audio_filter_get_type (void);

GST_AUDIO_API
void    gst_audio_filter_class_add_pad_templates (GstAudioFilterClass * klass,
                                                  GstCaps             * allowed_caps);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstAudioFilter, gst_object_unref)

G_END_DECLS

#endif /* __GST_AUDIO_FILTER_H__ */

