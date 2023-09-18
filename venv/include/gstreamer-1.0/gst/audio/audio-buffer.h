/* GStreamer
 * Copyright (C) <2018> Collabora Ltd.
 *   @author George Kiagiadakis <george.kiagiadakis@collabora.com>
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

#ifndef __GST_AUDIO_BUFFER_H__
#define __GST_AUDIO_BUFFER_H__

G_BEGIN_DECLS

/**
 * GstAudioBuffer:
 * @info: a #GstAudioInfo describing the audio properties of this buffer
 * @n_samples: the size of the buffer in samples
 * @n_planes: the number of planes available
 * @planes: an array of @n_planes pointers pointing to the start of each
 *   plane in the mapped buffer
 * @buffer: the mapped buffer
 *
 * A structure containing the result of an audio buffer map operation,
 * which is executed with gst_audio_buffer_map(). For non-interleaved (planar)
 * buffers, the beginning of each channel in the buffer has its own pointer in
 * the @planes array. For interleaved buffers, the @planes array only contains
 * one item, which is the pointer to the beginning of the buffer, and @n_planes
 * equals 1.
 *
 * The different channels in @planes are always in the GStreamer channel order.
 *
 * Since: 1.16
 */
typedef struct {
  GstAudioInfo info;

  gsize        n_samples;
  gint         n_planes;
  gpointer     *planes;

  GstBuffer    *buffer;

  /*< private >*/
  GstMapInfo   *map_infos;
  gpointer     priv_planes_arr[8];
  GstMapInfo   priv_map_infos_arr[8];

  gpointer     _gst_reserved[GST_PADDING];
} GstAudioBuffer;


GST_AUDIO_API
gboolean gst_audio_buffer_map (GstAudioBuffer *buffer, const GstAudioInfo *info,
                               GstBuffer *gstbuffer, GstMapFlags flags);

GST_AUDIO_API
void gst_audio_buffer_unmap (GstAudioBuffer *buffer);


#define GST_AUDIO_BUFFER_FORMAT(b)          (GST_AUDIO_INFO_FORMAT(&(b)->info))
#define GST_AUDIO_BUFFER_CHANNELS(b)        (GST_AUDIO_INFO_CHANNELS(&(b)->info))
#define GST_AUDIO_BUFFER_LAYOUT(b)          (GST_AUDIO_INFO_LAYOUT(&(b)->info))
#define GST_AUDIO_BUFFER_RATE(b)            (GST_AUDIO_INFO_RATE(&(b)->info))

#define GST_AUDIO_BUFFER_WIDTH(b)           (GST_AUDIO_INFO_WIDTH(&(b)->info))
#define GST_AUDIO_BUFFER_DEPTH(b)           (GST_AUDIO_INFO_DEPTH(&(b)->info))
#define GST_AUDIO_BUFFER_SAMPLE_STRIDE(b)   (GST_AUDIO_INFO_WIDTH(&(b)->info) >> 3)
#define GST_AUDIO_BUFFER_BPS(b)             (GST_AUDIO_INFO_DEPTH(&(b)->info) >> 3)
#define GST_AUDIO_BUFFER_BPF(b)             (GST_AUDIO_INFO_BPF(&(b)->info))

#define GST_AUDIO_BUFFER_N_SAMPLES(b)       ((b)->n_samples)
#define GST_AUDIO_BUFFER_N_PLANES(b)        ((b)->n_planes)
#define GST_AUDIO_BUFFER_PLANE_DATA(b,p)    ((b)->planes[p])

/* the size of each plane in bytes */
#define GST_AUDIO_BUFFER_PLANE_SIZE(b)      \
    (GST_AUDIO_BUFFER_N_SAMPLES(b) * GST_AUDIO_BUFFER_SAMPLE_STRIDE(b) * \
     GST_AUDIO_BUFFER_CHANNELS(b) / GST_AUDIO_BUFFER_N_PLANES(b))

G_END_DECLS

#endif /* __GST_AUDIO_BUFFER_H__ */
