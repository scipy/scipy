/* GStreamer
 * Copyright (C) <2011> Wim Taymans <wim.taymans@gmail.com>
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

#ifndef __GST_AUDIO_META_H__
#define __GST_AUDIO_META_H__

#include <gst/audio/audio.h>

G_BEGIN_DECLS

#define GST_AUDIO_DOWNMIX_META_API_TYPE (gst_audio_downmix_meta_api_get_type())
#define GST_AUDIO_DOWNMIX_META_INFO  (gst_audio_downmix_meta_get_info())

typedef struct _GstAudioDownmixMeta GstAudioDownmixMeta;

/**
 * GstAudioDownmixMeta:
 * @meta: parent #GstMeta
 * @from_position: the channel positions of the source
 * @to_position: the channel positions of the destination
 * @from_channels: the number of channels of the source
 * @to_channels: the number of channels of the destination
 * @matrix: the matrix coefficients.
 *
 * Extra buffer metadata describing audio downmixing matrix. This metadata is
 * attached to audio buffers and contains a matrix to downmix the buffer number
 * of channels to @channels.
 *
 * @matrix is an two-dimensional array of @to_channels times @from_channels
 * coefficients, i.e. the i-th output channels is constructed by multiplicating
 * the input channels with the coefficients in @matrix[i] and taking the sum
 * of the results.
 */
struct _GstAudioDownmixMeta {
  GstMeta      meta;

  GstAudioChannelPosition *from_position;
  GstAudioChannelPosition *to_position;
  gint        from_channels, to_channels;
  gfloat       **matrix;
};

GST_AUDIO_API
GType gst_audio_downmix_meta_api_get_type (void);

GST_AUDIO_API
const GstMetaInfo * gst_audio_downmix_meta_get_info (void);

#define gst_buffer_get_audio_downmix_meta(b) ((GstAudioDownmixMeta*)gst_buffer_get_meta((b), GST_AUDIO_DOWNMIX_META_API_TYPE))
GST_AUDIO_API
GstAudioDownmixMeta * gst_buffer_get_audio_downmix_meta_for_channels    (GstBuffer *buffer,
                                                                         const GstAudioChannelPosition *to_position,
                                                                         gint                           to_channels);

GST_AUDIO_API
GstAudioDownmixMeta * gst_buffer_add_audio_downmix_meta (GstBuffer    *buffer,
                                                         const GstAudioChannelPosition *from_position,
                                                         gint                           from_channels,
                                                         const GstAudioChannelPosition *to_position,
                                                         gint                           to_channels,
                                                         const gfloat                 **matrix);


#define GST_AUDIO_CLIPPING_META_API_TYPE (gst_audio_clipping_meta_api_get_type())
#define GST_AUDIO_CLIPPING_META_INFO  (gst_audio_clipping_meta_get_info())

typedef struct _GstAudioClippingMeta GstAudioClippingMeta;

/**
 * GstAudioClippingMeta:
 * @meta: parent #GstMeta
 * @format: GstFormat of @start and @stop, GST_FORMAT_DEFAULT is samples
 * @start: Amount of audio to clip from start of buffer
 * @end: Amount of  to clip from end of buffer
 *
 * Extra buffer metadata describing how much audio has to be clipped from
 * the start or end of a buffer. This is used for compressed formats, where
 * the first frame usually has some additional samples due to encoder and
 * decoder delays, and the last frame usually has some additional samples to
 * be able to fill the complete last frame.
 *
 * This is used to ensure that decoded data in the end has the same amount of
 * samples, and multiply decoded streams can be gaplessly concatenated.
 *
 * Note: If clipping of the start is done by adjusting the segment, this meta
 * has to be dropped from buffers as otherwise clipping could happen twice.
 *
 * Since: 1.8
 */
struct _GstAudioClippingMeta {
  GstMeta   meta;

  GstFormat format;
  guint64   start;
  guint64   end;
};

GST_AUDIO_API
GType gst_audio_clipping_meta_api_get_type (void);

GST_AUDIO_API
const GstMetaInfo * gst_audio_clipping_meta_get_info (void);

#define gst_buffer_get_audio_clipping_meta(b) ((GstAudioClippingMeta*)gst_buffer_get_meta((b), GST_AUDIO_CLIPPING_META_API_TYPE))

GST_AUDIO_API
GstAudioClippingMeta * gst_buffer_add_audio_clipping_meta (GstBuffer *buffer,
                                                           GstFormat  format,
                                                           guint64    start,
                                                           guint64    end);


#define GST_AUDIO_META_API_TYPE (gst_audio_meta_api_get_type())
#define GST_AUDIO_META_INFO  (gst_audio_meta_get_info())

typedef struct _GstAudioMeta GstAudioMeta;

/**
 * GstAudioMeta:
 * @meta: parent #GstMeta
 * @info: the audio properties of the buffer
 * @samples: the number of valid samples in the buffer
 * @offsets: the offsets (in bytes) where each channel plane starts in the
 *   buffer or %NULL if the buffer has interleaved layout; if not %NULL, this
 *   is guaranteed to be an array of @info.channels elements
 *
 * Buffer metadata describing how data is laid out inside the buffer. This
 * is useful for non-interleaved (planar) buffers, where it is necessary to
 * have a place to store where each plane starts and how long each plane is.
 *
 * It is a requirement for non-interleaved buffers to have this metadata
 * attached and to be mapped with gst_audio_buffer_map() in order to ensure
 * correct handling of clipping and channel reordering.
 *
 * The different channels in @offsets are always in the GStreamer channel order.
 * Zero-copy channel reordering can be implemented by swapping the values in
 * @offsets.
 *
 * It is not allowed for channels to overlap in memory,
 * i.e. for each i in [0, channels), the range
 * [@offsets[i], @offsets[i] + @samples * sample_stride) must not overlap
 * with any other such range.
 *
 * It is, however, allowed to have parts of the buffer memory unused,
 * by using @offsets and @samples in such a way that leave gaps on it.
 * This is used to implement zero-copy clipping in non-interleaved buffers.
 *
 * Obviously, due to the above, it is not safe to infer the
 * number of valid samples from the size of the buffer. You should always
 * use the @samples variable of this metadata.
 *
 * Note that for interleaved audio it is not a requirement to have this
 * metadata attached and at the moment of writing, there is actually no use
 * case to do so. It is, however, allowed to attach it, for some potential
 * future use case.
 *
 * Since: 1.16
 */
struct _GstAudioMeta {
  GstMeta      meta;

  GstAudioInfo info;
  gsize        samples;
  gsize        *offsets;

  /*< private >*/
  gsize        priv_offsets_arr[8];
  gpointer     _gst_reserved[GST_PADDING];
};

GST_AUDIO_API
GType gst_audio_meta_api_get_type (void);

GST_AUDIO_API
const GstMetaInfo * gst_audio_meta_get_info (void);

#define gst_buffer_get_audio_meta(b) \
    ((GstAudioMeta*)gst_buffer_get_meta((b), GST_AUDIO_META_API_TYPE))

GST_AUDIO_API
GstAudioMeta * gst_buffer_add_audio_meta (GstBuffer *buffer,
                                          const GstAudioInfo *info,
                                          gsize samples, gsize offsets[]);

/**
 * GST_AUDIO_LEVEL_META_API_TYPE:
 *
 * The #GType associated with #GstAudioLevelMeta.
 *
 * Since: 1.20
 */
#define GST_AUDIO_LEVEL_META_API_TYPE  (gst_audio_level_meta_api_get_type())
/**
 * GST_AUDIO_LEVEL_META_INFO:
 *
 * The #GstMetaInfo associated with #GstAudioLevelMeta.
 *
 * Since: 1.20
 */
#define GST_AUDIO_LEVEL_META_INFO  (gst_audio_level_meta_get_info())
typedef struct _GstAudioLevelMeta GstAudioLevelMeta;

/**
 * GstAudioLevelMeta:
 * @meta: parent #GstMeta
 * @level: the -dBov from 0-127 (127 is silence).
 * @voice_activity: whether the buffer contains voice activity
 *
 * Meta containing Audio Level Indication: https://tools.ietf.org/html/rfc6464
 *
 * Since: 1.20
 */
struct _GstAudioLevelMeta
{
  GstMeta meta;

  guint8 level;
  gboolean voice_activity;
};

GST_AUDIO_API
GType                  gst_audio_level_meta_api_get_type                (void);

GST_AUDIO_API
const GstMetaInfo *    gst_audio_level_meta_get_info                    (void);

GST_AUDIO_API
GstAudioLevelMeta * gst_buffer_add_audio_level_meta                     (GstBuffer * buffer,
                                                                         guint8 level,
                                                                         gboolean voice_activity);
GST_AUDIO_API
GstAudioLevelMeta * gst_buffer_get_audio_level_meta                     (GstBuffer * buffer);

G_END_DECLS

#endif /* __GST_AUDIO_META_H__ */
