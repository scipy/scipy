/* GStreamer StreamVolume
 * Copyright (C) 2009 Sebastian Dröge <sebastian.droege@collabora.co.uk>
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

#ifndef __GST_STREAM_VOLUME_H__
#define __GST_STREAM_VOLUME_H__

#include <gst/gst.h>
#include <gst/audio/audio-prelude.h>

G_BEGIN_DECLS

#define GST_TYPE_STREAM_VOLUME (gst_stream_volume_get_type ())
GST_AUDIO_API
G_DECLARE_INTERFACE (GstStreamVolume, gst_stream_volume, GST, STREAM_VOLUME,
    GObject)

#define GST_STREAM_VOLUME_GET_INTERFACE(obj) GST_STREAM_VOLUME_GET_IFACE(obj)

struct _GstStreamVolumeInterface {
  GTypeInterface iface;
};

/**
 * GstStreamVolumeFormat:
 * @GST_STREAM_VOLUME_FORMAT_LINEAR: Linear scale factor, 1.0 = 100%
 * @GST_STREAM_VOLUME_FORMAT_CUBIC: Cubic volume scale
 * @GST_STREAM_VOLUME_FORMAT_DB: Logarithmic volume scale (dB, amplitude not power)
 *
 * Different representations of a stream volume. gst_stream_volume_convert_volume()
 * allows to convert between the different representations.
 *
 * Formulas to convert from a linear to a cubic or dB volume are
 * cbrt(val) and 20 * log10 (val).
 */
typedef enum {
  GST_STREAM_VOLUME_FORMAT_LINEAR = 0,
  GST_STREAM_VOLUME_FORMAT_CUBIC,
  GST_STREAM_VOLUME_FORMAT_DB
} GstStreamVolumeFormat;

GST_AUDIO_API
void            gst_stream_volume_set_volume      (GstStreamVolume *volume,
                                                   GstStreamVolumeFormat format,
                                                   gdouble val);

GST_AUDIO_API
gdouble         gst_stream_volume_get_volume      (GstStreamVolume *volume,
                                                   GstStreamVolumeFormat format);

GST_AUDIO_API
void            gst_stream_volume_set_mute        (GstStreamVolume *volume,
                                                   gboolean mute);

GST_AUDIO_API
gboolean        gst_stream_volume_get_mute        (GstStreamVolume *volume);

GST_AUDIO_API
gdouble         gst_stream_volume_convert_volume  (GstStreamVolumeFormat from,
                                                   GstStreamVolumeFormat to,
                                                   gdouble val) G_GNUC_CONST;

G_END_DECLS

#endif /* __GST_STREAM_VOLUME_H__ */
