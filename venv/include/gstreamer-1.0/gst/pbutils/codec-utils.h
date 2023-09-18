/* GStreamer base utils library codec-specific utility functions
 * Copyright (C) 2010 Arun Raghavan <arun.raghavan@collabora.co.uk>
 *               2010 Collabora Multimedia
 *               2010 Nokia Corporation
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

#ifndef __GST_PB_UTILS_CODEC_UTILS_H__
#define __GST_PB_UTILS_CODEC_UTILS_H__

#include <gst/gst.h>
#include <gst/pbutils/pbutils-prelude.h>

G_BEGIN_DECLS

/* AAC */

GST_PBUTILS_API
guint         gst_codec_utils_aac_get_sample_rate_from_index (guint sr_idx);

GST_PBUTILS_API
gint          gst_codec_utils_aac_get_index_from_sample_rate (guint rate);

GST_PBUTILS_API
const gchar * gst_codec_utils_aac_get_profile (const guint8 * audio_config, guint len);

GST_PBUTILS_API
const gchar * gst_codec_utils_aac_get_level   (const guint8 * audio_config, guint len);

GST_PBUTILS_API
guint         gst_codec_utils_aac_get_sample_rate (const guint8 * audio_config, guint len);

GST_PBUTILS_API
guint         gst_codec_utils_aac_get_channels (const guint8 * audio_config, guint len);

GST_PBUTILS_API
gboolean      gst_codec_utils_aac_caps_set_level_and_profile (GstCaps      * caps,
                                                              const guint8 * audio_config,
                                                              guint          len);

/* H.264 */

GST_PBUTILS_API
const gchar * gst_codec_utils_h264_get_profile (const guint8 * sps, guint len);

GST_PBUTILS_API
const gchar * gst_codec_utils_h264_get_level   (const guint8 * sps, guint len);

GST_PBUTILS_API
guint8        gst_codec_utils_h264_get_level_idc (const gchar * level);

GST_PBUTILS_API
gboolean      gst_codec_utils_h264_caps_set_level_and_profile (GstCaps      * caps,
                                                               const guint8 * sps,
                                                               guint          len);

GST_PBUTILS_API
gboolean      gst_codec_utils_h264_get_profile_flags_level (const guint8 * codec_data,
                                                            guint len,
                                                            guint8 * profile,
                                                            guint8 * flags,
                                                            guint8 * level);

/* H.265 */

GST_PBUTILS_API
const gchar * gst_codec_utils_h265_get_profile                     (const guint8 * profile_tier_level,
                                                                    guint len);

GST_PBUTILS_API
const gchar * gst_codec_utils_h265_get_tier                        (const guint8 * profile_tier_level,
                                                                    guint len);

GST_PBUTILS_API
const gchar * gst_codec_utils_h265_get_level                       (const guint8 * profile_tier_level,
                                                                    guint len);

GST_PBUTILS_API
guint8        gst_codec_utils_h265_get_level_idc                   (const gchar  * level);

GST_PBUTILS_API
gboolean      gst_codec_utils_h265_caps_set_level_tier_and_profile (GstCaps      * caps,
                                                                    const guint8 * profile_tier_level,
                                                                    guint          len);
/* MPEG-4 part 2 */

GST_PBUTILS_API
const gchar * gst_codec_utils_mpeg4video_get_profile (const guint8 * vis_obj_seq, guint len);

GST_PBUTILS_API
const gchar * gst_codec_utils_mpeg4video_get_level   (const guint8 * vis_obj_seq, guint len);

GST_PBUTILS_API
gboolean      gst_codec_utils_mpeg4video_caps_set_level_and_profile (GstCaps      * caps,
                                                                     const guint8 * vis_obj_seq,
                                                                     guint          len);

/* Opus */

GST_PBUTILS_API
gboolean      gst_codec_utils_opus_parse_caps (GstCaps   * caps,
                                               guint32   * rate,
                                               guint8    * channels,
                                               guint8    * channel_mapping_family,
                                               guint8    * stream_count,
                                               guint8    * coupled_count,
                                               guint8      channel_mapping[256]);

GST_PBUTILS_API
GstCaps *     gst_codec_utils_opus_create_caps (guint32        rate,
                                                guint8         channels,
                                                guint8         channel_mapping_family,
                                                guint8         stream_count,
                                                guint8         coupled_count,
                                                const guint8 * channel_mapping);

GST_PBUTILS_API
GstCaps *    gst_codec_utils_opus_create_caps_from_header (GstBuffer * header, GstBuffer * comments);

GST_PBUTILS_API
GstBuffer *  gst_codec_utils_opus_create_header (guint32        rate,
                                                 guint8         channels,
                                                 guint8         channel_mapping_family,
                                                 guint8         stream_count,
                                                 guint8         coupled_count,
                                                 const guint8 * channel_mapping,
                                                 guint16        pre_skip,
                                                 gint16         output_gain);

GST_PBUTILS_API
gboolean  gst_codec_utils_opus_parse_header (GstBuffer * header,
                                             guint32   * rate,
                                             guint8    * channels,
                                             guint8    * channel_mapping_family,
                                             guint8    * stream_count,
                                             guint8    * coupled_count,
                                             guint8      channel_mapping[256],
                                             guint16   * pre_skip,
                                             gint16    * output_gain);

/* General */
GST_PBUTILS_API
gchar * gst_codec_utils_caps_get_mime_codec (GstCaps * caps);

GST_PBUTILS_API
GstCaps * gst_codec_utils_caps_from_mime_codec (const gchar *codecs_field);

G_END_DECLS

#endif /* __GST_PB_UTILS_CODEC_UTILS_H__ */
