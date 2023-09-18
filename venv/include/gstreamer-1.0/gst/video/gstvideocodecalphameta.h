/* GStreamer
 * Copyright (C) 2021 Collabora Ltd.
 *   Author: Nicolas Dufresne <nicolas.dufresne@collabora.com>
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

#ifndef __GST_VIDEO_CODEC_ALPHA_META_H__
#define __GST_VIDEO_CODEC_ALPHA_META_H__

#include <gst/gst.h>
#include <gst/video/video.h>

G_BEGIN_DECLS

/**
 * GST_VIDEO_CODEC_ALPHA_META_API_TYPE:
 *
 * Since: 1.20
 */
#define GST_VIDEO_CODEC_ALPHA_META_API_TYPE (gst_video_codec_alpha_meta_api_get_type())

/**
 * GST_VIDEO_CODEC_ALPHA_META_INFO:
 *
 * Since: 1.20
 */
#define GST_VIDEO_CODEC_ALPHA_META_INFO  (gst_video_codec_alpha_meta_get_info())

typedef struct _GstVideoCodecAlphaMeta GstVideoCodecAlphaMeta;

/**
 * GstVideoCodecAlphaMeta:
 * @meta: parent #GstMeta
 * @buffer: the encoded alpha frame
 *
 * Encapsulate an extra frame containing the encoded alpha channel for the
 * currently negotiated CODEC. The streams must be of the same dimention as
 * the original one.
 *
 * Since: 1.20
 */
struct _GstVideoCodecAlphaMeta
{
  GstMeta meta;

  GstBuffer *buffer;
};

GST_VIDEO_API
GType gst_video_codec_alpha_meta_api_get_type          (void);

GST_VIDEO_API
const GstMetaInfo *gst_video_codec_alpha_meta_get_info (void);

/**
 * gst_buffer_get_video_codec_alpha_meta:
 * @b: A #GstBuffer pointer, must be writable.
 *
 * Helper macro to get #GstVideoCodecAlphaMeta from an existing #GstBuffer.
 *
 * Returns: (nullable): the #GstVideoCodecAlphaMeta pointer, or %NULL if none.
 *
 * Since: 1.20
 */
#define gst_buffer_get_video_codec_alpha_meta(b) \
    ((GstVideoCodecAlphaMeta *)gst_buffer_get_meta((b),GST_VIDEO_CODEC_ALPHA_META_API_TYPE))

GST_VIDEO_API
GstVideoCodecAlphaMeta *gst_buffer_add_video_codec_alpha_meta (GstBuffer * buffer,
                                                               GstBuffer * alpha_buffer);

G_END_DECLS

#endif /* __GST_VIDEO_CODEC_ALPHA_META_H__ */
