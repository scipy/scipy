/* GStreamer
 * Copyright (C) <2005> Philippe Khalaf <burger@speedy.org>
 *               <2005> Wim Taymans <wim@fluendo.com>
 *
 * gstrtpbuffer.h: various helper functions to manipulate buffers
 *     with RTP payload.
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

#ifndef __GST_RTPBUFFER_H__
#define __GST_RTPBUFFER_H__

#include <gst/gst.h>
#include <gst/rtp/gstrtppayloads.h>

G_BEGIN_DECLS

/**
 * GST_RTP_VERSION:
 *
 * The supported RTP version 2.
 */
#define GST_RTP_VERSION 2


typedef struct _GstRTPBuffer GstRTPBuffer;

/**
 * GstRTPBuffer:
 * @buffer: pointer to RTP buffer
 * @state: internal state
 * @data: array of data
 * @size: array of size
 * @map: array of #GstMapInfo
 *
 * Data structure that points to an RTP packet.
 * The size of the structure is made public to allow stack allocations.
 */
struct _GstRTPBuffer
{
  GstBuffer   *buffer;
  guint        state;
  gpointer     data[4];
  gsize        size[4];
  GstMapInfo   map[4];
};

#define GST_RTP_BUFFER_INIT { NULL, 0, { NULL, NULL, NULL, NULL}, { 0, 0, 0, 0 }, \
  { GST_MAP_INFO_INIT, GST_MAP_INFO_INIT, GST_MAP_INFO_INIT, GST_MAP_INFO_INIT} }

/* creating buffers */

GST_RTP_API
void            gst_rtp_buffer_allocate_data         (GstBuffer *buffer, guint payload_len,
                                                      guint8 pad_len, guint8 csrc_count);

GST_RTP_API
GstBuffer*      gst_rtp_buffer_new_take_data         (gpointer data, gsize len);

GST_RTP_API
GstBuffer*      gst_rtp_buffer_new_copy_data         (gconstpointer data, gsize len);

GST_RTP_API
GstBuffer*      gst_rtp_buffer_new_allocate          (guint payload_len, guint8 pad_len, guint8 csrc_count);

GST_RTP_API
GstBuffer*      gst_rtp_buffer_new_allocate_len      (guint packet_len, guint8 pad_len, guint8 csrc_count);

GST_RTP_API
guint           gst_rtp_buffer_calc_header_len       (guint8 csrc_count);

GST_RTP_API
guint           gst_rtp_buffer_calc_packet_len       (guint payload_len, guint8 pad_len, guint8 csrc_count);

GST_RTP_API
guint           gst_rtp_buffer_calc_payload_len      (guint packet_len, guint8 pad_len, guint8 csrc_count);

GST_RTP_API
gboolean        gst_rtp_buffer_map                   (GstBuffer *buffer, GstMapFlags flags, GstRTPBuffer *rtp);

GST_RTP_API
void            gst_rtp_buffer_unmap                 (GstRTPBuffer *rtp);

GST_RTP_API
void            gst_rtp_buffer_set_packet_len        (GstRTPBuffer *rtp, guint len);

GST_RTP_API
guint           gst_rtp_buffer_get_packet_len        (GstRTPBuffer *rtp);

GST_RTP_API
guint           gst_rtp_buffer_get_header_len        (GstRTPBuffer *rtp);

GST_RTP_API
guint8          gst_rtp_buffer_get_version           (GstRTPBuffer *rtp);

GST_RTP_API
void            gst_rtp_buffer_set_version           (GstRTPBuffer *rtp, guint8 version);

GST_RTP_API
gboolean        gst_rtp_buffer_get_padding           (GstRTPBuffer *rtp);

GST_RTP_API
void            gst_rtp_buffer_set_padding           (GstRTPBuffer *rtp, gboolean padding);

GST_RTP_API
void            gst_rtp_buffer_pad_to                (GstRTPBuffer *rtp, guint len);

GST_RTP_API
gboolean        gst_rtp_buffer_get_extension         (GstRTPBuffer *rtp);

GST_RTP_API
void            gst_rtp_buffer_set_extension         (GstRTPBuffer *rtp, gboolean extension);

GST_RTP_API
gboolean        gst_rtp_buffer_get_extension_data    (GstRTPBuffer *rtp, guint16 *bits,
                                                      gpointer *data, guint *wordlen);

GST_RTP_API
GBytes*         gst_rtp_buffer_get_extension_bytes   (GstRTPBuffer *rtp, guint16 *bits);

GST_RTP_API
gboolean        gst_rtp_buffer_set_extension_data    (GstRTPBuffer *rtp, guint16 bits, guint16 length);

GST_RTP_API
void            gst_rtp_buffer_remove_extension_data (GstRTPBuffer *rtp);

GST_RTP_API
guint32         gst_rtp_buffer_get_ssrc              (GstRTPBuffer *rtp);

GST_RTP_API
void            gst_rtp_buffer_set_ssrc              (GstRTPBuffer *rtp, guint32 ssrc);

GST_RTP_API
guint8          gst_rtp_buffer_get_csrc_count        (GstRTPBuffer *rtp);

GST_RTP_API
guint32         gst_rtp_buffer_get_csrc              (GstRTPBuffer *rtp, guint8 idx);

GST_RTP_API
void            gst_rtp_buffer_set_csrc              (GstRTPBuffer *rtp, guint8 idx, guint32 csrc);

GST_RTP_API
gboolean        gst_rtp_buffer_get_marker            (GstRTPBuffer *rtp);

GST_RTP_API
void            gst_rtp_buffer_set_marker            (GstRTPBuffer *rtp, gboolean marker);

GST_RTP_API
guint8          gst_rtp_buffer_get_payload_type      (GstRTPBuffer *rtp);

GST_RTP_API
void            gst_rtp_buffer_set_payload_type      (GstRTPBuffer *rtp, guint8 payload_type);

GST_RTP_API
guint16         gst_rtp_buffer_get_seq               (GstRTPBuffer *rtp);

GST_RTP_API
void            gst_rtp_buffer_set_seq               (GstRTPBuffer *rtp, guint16 seq);

GST_RTP_API
guint32         gst_rtp_buffer_get_timestamp         (GstRTPBuffer *rtp);

GST_RTP_API
void            gst_rtp_buffer_set_timestamp         (GstRTPBuffer *rtp, guint32 timestamp);

GST_RTP_API
GstBuffer*      gst_rtp_buffer_get_payload_buffer    (GstRTPBuffer *rtp);

GST_RTP_API
GstBuffer*      gst_rtp_buffer_get_payload_subbuffer (GstRTPBuffer *rtp, guint offset, guint len);

GST_RTP_API
guint           gst_rtp_buffer_get_payload_len       (GstRTPBuffer *rtp);

GST_RTP_API
gpointer        gst_rtp_buffer_get_payload           (GstRTPBuffer *rtp);

GST_RTP_API
GBytes*         gst_rtp_buffer_get_payload_bytes     (GstRTPBuffer *rtp);

/* some helpers */

GST_RTP_API
guint32         gst_rtp_buffer_default_clock_rate    (guint8 payload_type);

GST_RTP_API
gint            gst_rtp_buffer_compare_seqnum        (guint16 seqnum1, guint16 seqnum2);

GST_RTP_API
guint64         gst_rtp_buffer_ext_timestamp         (guint64 *exttimestamp, guint32 timestamp);

GST_RTP_API
gboolean        gst_rtp_buffer_get_extension_onebyte_header  (GstRTPBuffer *rtp,
                                                              guint8 id,
                                                              guint nth,
                                                              gpointer * data,
                                                              guint * size);

GST_RTP_API
gboolean        gst_rtp_buffer_get_extension_twobytes_header (GstRTPBuffer *rtp,
                                                              guint8 * appbits,
                                                              guint8 id,
                                                              guint nth,
                                                              gpointer * data,
                                                              guint * size);

GST_RTP_API
gboolean       gst_rtp_buffer_add_extension_onebyte_header  (GstRTPBuffer *rtp,
                                                             guint8 id,
                                                             gconstpointer data,
                                                             guint size);

GST_RTP_API
gboolean       gst_rtp_buffer_add_extension_twobytes_header (GstRTPBuffer *rtp,
                                                             guint8 appbits,
                                                             guint8 id,
                                                             gconstpointer data,
                                                             guint size);

GST_RTP_API
gboolean gst_rtp_buffer_get_extension_onebyte_header_from_bytes (GBytes * bytes,
                                                                 guint16 bit_pattern,
                                                                 guint8 id,
                                                                 guint nth,
                                                                 gpointer * data,
                                                                 guint * size);

/**
 * GstRTPBufferFlags:
 * @GST_RTP_BUFFER_FLAG_RETRANSMISSION: The #GstBuffer was once wrapped
 *           in a retransmitted packet as specified by RFC 4588.
 * @GST_RTP_BUFFER_FLAG_REDUNDANT:      The packet represents redundant RTP packet.
 *           The flag is used in gstrtpstorage to be able to hold the packetback
 *           and use it only for recovery from packet loss.
 *           Since: 1.14
 * @GST_RTP_BUFFER_FLAG_LAST:           Offset to define more flags.
 *
 * Additional RTP buffer flags. These flags can potentially be used on any
 * buffers carrying RTP packets.
 *
 * Note that these are only valid for #GstCaps of type: application/x-rtp (x-rtcp).
 * They can conflict with other extended buffer flags.
 *
 * Since: 1.10
 */
typedef enum {
  GST_RTP_BUFFER_FLAG_RETRANSMISSION = (GST_BUFFER_FLAG_LAST << 0),
  GST_RTP_BUFFER_FLAG_REDUNDANT      = (GST_BUFFER_FLAG_LAST << 1),
  GST_RTP_BUFFER_FLAG_LAST           = (GST_BUFFER_FLAG_LAST << 8)
} GstRTPBufferFlags;

/**
 * GstRTPBufferMapFlags:
 * @GST_RTP_BUFFER_MAP_FLAG_SKIP_PADDING: Skip mapping and validation of RTP
 *           padding and RTP pad count when present. Useful for buffers where
 *           the padding may be encrypted.
 * @GST_RTP_BUFFER_MAP_FLAG_LAST: Offset to define more flags
 *
 * Additional mapping flags for gst_rtp_buffer_map().
 *
 * Since: 1.6.1
 */
typedef enum {
  GST_RTP_BUFFER_MAP_FLAG_SKIP_PADDING = (GST_MAP_FLAG_LAST << 0),
  GST_RTP_BUFFER_MAP_FLAG_LAST         = (GST_MAP_FLAG_LAST << 8)
  /* 8 more flags possible afterwards */
} GstRTPBufferMapFlags;

G_END_DECLS

#endif /* __GST_RTPBUFFER_H__ */

