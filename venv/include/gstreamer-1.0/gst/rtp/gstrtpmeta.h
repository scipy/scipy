/* GStreamer
 * Copyright (C) <2016> Stian Selnes <stian@pexip.com>
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

#ifndef __GST_RTP_META_H__
#define __GST_RTP_META_H__

#include <gst/gst.h>
#include <gst/rtp/rtp-prelude.h>

G_BEGIN_DECLS

#define GST_RTP_SOURCE_META_API_TYPE  (gst_rtp_source_meta_api_get_type())
#define GST_RTP_SOURCE_META_INFO  (gst_rtp_source_meta_get_info())
typedef struct _GstRTPSourceMeta GstRTPSourceMeta;

#define GST_RTP_SOURCE_META_MAX_CSRC_COUNT 15

/**
 * GstRTPSourceMeta:
 * @meta: parent #GstMeta
 * @ssrc: the SSRC
 * @ssrc_valid: whether @ssrc is set and valid
 * @csrc: (allow-none): pointer to the CSRCs
 * @csrc_count: number of elements in @csrc
 *
 * Meta describing the source(s) of the buffer.
 *
 * Since: 1.16
 */
struct _GstRTPSourceMeta
{
  GstMeta meta;

  guint32 ssrc;
  gboolean ssrc_valid;
  guint32 csrc[GST_RTP_SOURCE_META_MAX_CSRC_COUNT];
  guint csrc_count;
};

GST_RTP_API
GType               gst_rtp_source_meta_api_get_type     (void);

GST_RTP_API
GstRTPSourceMeta *  gst_buffer_add_rtp_source_meta       (GstBuffer * buffer, const guint32 * ssrc,
                                                          const guint32 * csrc, guint csrc_count);
GST_RTP_API
GstRTPSourceMeta *  gst_buffer_get_rtp_source_meta       (GstBuffer * buffer);

GST_RTP_API
guint               gst_rtp_source_meta_get_source_count (const GstRTPSourceMeta * meta);

GST_RTP_API
gboolean            gst_rtp_source_meta_set_ssrc         (GstRTPSourceMeta * meta, guint32 * ssrc);

GST_RTP_API
gboolean            gst_rtp_source_meta_append_csrc      (GstRTPSourceMeta * meta,
                                                          const guint32 * csrc, guint csrc_count);
GST_RTP_API
const GstMetaInfo * gst_rtp_source_meta_get_info         (void);

G_END_DECLS

#endif /* __GST_RTP_META_H__ */
