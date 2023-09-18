/* GStreamer
 * Copyright (C) <2007> Wim Taymans <wim@fluendo.com>
 *
 * gstrtcpbuffer.h: various helper functions to manipulate buffers
 *     with RTCP payload.
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

#ifndef __GST_RTCPBUFFER_H__
#define __GST_RTCPBUFFER_H__

#include <gst/gst.h>
#include <gst/rtp/rtp-prelude.h>

G_BEGIN_DECLS

/**
 * GST_RTCP_VERSION:
 *
 * The supported RTCP version 2.
 */
#define GST_RTCP_VERSION 2

/**
 * GstRTCPType:
 * @GST_RTCP_TYPE_INVALID: Invalid type
 * @GST_RTCP_TYPE_SR: Sender report
 * @GST_RTCP_TYPE_RR: Receiver report
 * @GST_RTCP_TYPE_SDES: Source description
 * @GST_RTCP_TYPE_BYE: Goodbye
 * @GST_RTCP_TYPE_APP: Application defined
 * @GST_RTCP_TYPE_RTPFB: Transport layer feedback.
 * @GST_RTCP_TYPE_PSFB: Payload-specific feedback.
 * @GST_RTCP_TYPE_XR: Extended report.
 *
 * Different RTCP packet types.
 */
typedef enum
{
  GST_RTCP_TYPE_INVALID = 0,
  GST_RTCP_TYPE_SR      = 200,
  GST_RTCP_TYPE_RR      = 201,
  GST_RTCP_TYPE_SDES    = 202,
  GST_RTCP_TYPE_BYE     = 203,
  GST_RTCP_TYPE_APP     = 204,
  GST_RTCP_TYPE_RTPFB   = 205,
  GST_RTCP_TYPE_PSFB    = 206,
  GST_RTCP_TYPE_XR      = 207
} GstRTCPType;

/* FIXME 2.0: backwards compatibility define for enum typo */
#define GST_RTCP_RTPFB_TYPE_RCTP_SR_REQ GST_RTCP_RTPFB_TYPE_RTCP_SR_REQ

/**
 * GstRTCPFBType:
 * @GST_RTCP_FB_TYPE_INVALID: Invalid type
 * @GST_RTCP_RTPFB_TYPE_NACK: Generic NACK
 * @GST_RTCP_RTPFB_TYPE_TMMBR: Temporary Maximum Media Stream Bit Rate Request
 * @GST_RTCP_RTPFB_TYPE_TMMBN: Temporary Maximum Media Stream Bit Rate
 *    Notification
 * @GST_RTCP_RTPFB_TYPE_RTCP_SR_REQ: Request an SR packet for early
 *    synchronization
 * @GST_RTCP_PSFB_TYPE_PLI: Picture Loss Indication
 * @GST_RTCP_PSFB_TYPE_SLI: Slice Loss Indication
 * @GST_RTCP_PSFB_TYPE_RPSI: Reference Picture Selection Indication
 * @GST_RTCP_PSFB_TYPE_AFB: Application layer Feedback
 * @GST_RTCP_PSFB_TYPE_FIR: Full Intra Request Command
 * @GST_RTCP_PSFB_TYPE_TSTR: Temporal-Spatial Trade-off Request
 * @GST_RTCP_PSFB_TYPE_TSTN: Temporal-Spatial Trade-off Notification
 * @GST_RTCP_PSFB_TYPE_VBCN: Video Back Channel Message
 *
 * Different types of feedback messages.
 */
typedef enum
{
  /* generic */
  GST_RTCP_FB_TYPE_INVALID        = 0,
  /* RTPFB types */
  GST_RTCP_RTPFB_TYPE_NACK        = 1,
  /* RTPFB types assigned in RFC 5104 */
  GST_RTCP_RTPFB_TYPE_TMMBR       = 3,
  GST_RTCP_RTPFB_TYPE_TMMBN       = 4,
  /* RTPFB types assigned in RFC 6051 */
  GST_RTCP_RTPFB_TYPE_RTCP_SR_REQ = 5,
  /* draft-holmer-rmcat-transport-wide-cc-extensions-01 */
  GST_RTCP_RTPFB_TYPE_TWCC         = 15,

  /* PSFB types */
  GST_RTCP_PSFB_TYPE_PLI          = 1,
  GST_RTCP_PSFB_TYPE_SLI          = 2,
  GST_RTCP_PSFB_TYPE_RPSI         = 3,
  GST_RTCP_PSFB_TYPE_AFB          = 15,
  /* PSFB types assigned in RFC 5104 */
  GST_RTCP_PSFB_TYPE_FIR          = 4,
  GST_RTCP_PSFB_TYPE_TSTR         = 5,
  GST_RTCP_PSFB_TYPE_TSTN         = 6,
  GST_RTCP_PSFB_TYPE_VBCN         = 7,
} GstRTCPFBType;

/**
 * GstRTCPSDESType:
 * @GST_RTCP_SDES_INVALID: Invalid SDES entry
 * @GST_RTCP_SDES_END: End of SDES list
 * @GST_RTCP_SDES_CNAME: Canonical name
 * @GST_RTCP_SDES_NAME: User name
 * @GST_RTCP_SDES_EMAIL: User's electronic mail address
 * @GST_RTCP_SDES_PHONE: User's phone number
 * @GST_RTCP_SDES_LOC: Geographic user location
 * @GST_RTCP_SDES_TOOL: Name of application or tool
 * @GST_RTCP_SDES_NOTE: Notice about the source
 * @GST_RTCP_SDES_PRIV: Private extensions
 *
 * Different types of SDES content.
 */
/**
 * GST_RTCP_SDES_H323_CADDR:
 *
 * H.323 callable address
 *
 * Since: 1.20:
 */
/**
 * GST_RTCP_SDES_APSI:
 *
 * Application Specific Identifier (RFC6776)
 *
 * Since: 1.20:
 */
/**
 * GST_RTCP_SDES_RGRP:
 *
 * Reporting Group Identifier (RFC8861)
 *
 * Since: 1.20:
 */
/**
 * GST_RTCP_SDES_RTP_STREAM_ID:
 *
 * RtpStreamId SDES item (RFC8852).
 *
 * Since: 1.20:
 */
/**
 * GST_RTCP_SDES_REPAIRED_RTP_STREAM_ID:
 *
 * RepairedRtpStreamId SDES item (RFC8852).
 *
 * Since: 1.20:
 */
/**
 * GST_RTCP_SDES_CCID:
 *
 * CLUE CaptId (RFC8849)
 *
 * Since: 1.20:
 */
/**
 * GST_RTCP_SDES_MID:
 *
 * MID SDES item (RFC8843).
 *
 * Since: 1.20:
 */
typedef enum
{
  GST_RTCP_SDES_INVALID                 = -1,
  GST_RTCP_SDES_END                     = 0,
  GST_RTCP_SDES_CNAME                   = 1,
  GST_RTCP_SDES_NAME                    = 2,
  GST_RTCP_SDES_EMAIL                   = 3,
  GST_RTCP_SDES_PHONE                   = 4,
  GST_RTCP_SDES_LOC                     = 5,
  GST_RTCP_SDES_TOOL                    = 6,
  GST_RTCP_SDES_NOTE                    = 7,
  GST_RTCP_SDES_PRIV                    = 8,
  GST_RTCP_SDES_H323_CADDR              = 9,
  GST_RTCP_SDES_APSI                    = 10,
  GST_RTCP_SDES_RGRP                    = 11,
  GST_RTCP_SDES_RTP_STREAM_ID           = 12,
  GST_RTCP_SDES_REPAIRED_RTP_STREAM_ID  = 13,
  GST_RTCP_SDES_CCID                    = 14,
  GST_RTCP_SDES_MID                     = 15,
} GstRTCPSDESType;

/**
 * GstRTCPXRType:
 * @GST_RTCP_XR_TYPE_INVALID: Invalid XR Report Block
 * @GST_RTCP_XR_TYPE_LRLE: Loss RLE Report Block
 * @GST_RTCP_XR_TYPE_DRLE: Duplicate RLE Report Block
 * @GST_RTCP_XR_TYPE_PRT: Packet Receipt Times Report Block
 * @GST_RTCP_XR_TYPE_RRT: Receiver Reference Time Report Block
 * @GST_RTCP_XR_TYPE_DLRR: Delay since the last Receiver Report
 * @GST_RTCP_XR_TYPE_SSUMM: Statistics Summary Report Block
 * @GST_RTCP_XR_TYPE_VOIP_METRICS: VoIP Metrics Report Block
 *
 * Types of RTCP Extended Reports, those are defined in RFC 3611 and other RFCs
 * according to the [IANA registry](https://www.iana.org/assignments/rtcp-xr-block-types/rtcp-xr-block-types.xhtml).
 *
 * Since: 1.16
 */
typedef enum
{
  GST_RTCP_XR_TYPE_INVALID      = -1,
  GST_RTCP_XR_TYPE_LRLE         = 1,
  GST_RTCP_XR_TYPE_DRLE         = 2,
  GST_RTCP_XR_TYPE_PRT          = 3,
  GST_RTCP_XR_TYPE_RRT          = 4,
  GST_RTCP_XR_TYPE_DLRR         = 5,
  GST_RTCP_XR_TYPE_SSUMM        = 6,
  GST_RTCP_XR_TYPE_VOIP_METRICS = 7
} GstRTCPXRType;

/**
 * GST_RTCP_MAX_SDES:
 *
 * The maximum text length for an SDES item.
 */
#define GST_RTCP_MAX_SDES 255

/**
 * GST_RTCP_MAX_RB_COUNT:
 *
 * The maximum amount of Receiver report blocks in RR and SR messages.
 */
#define GST_RTCP_MAX_RB_COUNT   31

/**
 * GST_RTCP_MAX_SDES_ITEM_COUNT:
 *
 * The maximum amount of SDES items.
 */
#define GST_RTCP_MAX_SDES_ITEM_COUNT   31

/**
 * GST_RTCP_MAX_BYE_SSRC_COUNT:
 *
 * The maximum amount of SSRCs in a BYE packet.
 */
#define GST_RTCP_MAX_BYE_SSRC_COUNT   31

/**
 * GST_RTCP_VALID_MASK:
 *
 * Mask for version, padding bit and packet type pair
 */
#define GST_RTCP_VALID_MASK (0xc000 | 0x2000 | 0xfe)

/**
 * GST_RTCP_REDUCED_SIZE_VALID_MASK:
 *
 * Mask for version and packet type pair allowing reduced size
 * packets, basically it accepts other types than RR and SR
 */
#define GST_RTCP_REDUCED_SIZE_VALID_MASK (0xc000 | 0xf8)

/**
 * GST_RTCP_VALID_VALUE:
 *
 * Valid value for the first two bytes of an RTCP packet after applying
 * #GST_RTCP_VALID_MASK to them.
 */
#define GST_RTCP_VALID_VALUE ((GST_RTCP_VERSION << 14) | GST_RTCP_TYPE_SR)

typedef struct _GstRTCPBuffer GstRTCPBuffer;
typedef struct _GstRTCPPacket GstRTCPPacket;

struct _GstRTCPBuffer
{
  GstBuffer   *buffer;
  GstMapInfo   map;
};

#define GST_RTCP_BUFFER_INIT { NULL, GST_MAP_INFO_INIT }

/**
 * GstRTCPPacket:
 * @rtcp: pointer to RTCP buffer
 * @offset: offset of packet in buffer data
 *
 * Data structure that points to a packet at @offset in @buffer.
 * The size of the structure is made public to allow stack allocations.
 */
struct _GstRTCPPacket
{
  /*< public >*/
  GstRTCPBuffer *rtcp;
  guint          offset;

  /*< private >*/
  gboolean       padding;      /* padding field of current packet */
  guint8         count;        /* count field of current packet */
  GstRTCPType    type;         /* type of current packet */
  guint16        length;       /* length of current packet in 32-bits words minus one, this is validated when doing _get_first_packet() and _move_to_next() */

  guint          item_offset;  /* current item offset for navigating SDES */
  guint          item_count;   /* current item count */
  guint          entry_offset; /* current entry offset for navigating SDES items */
};

/* creating buffers */

GST_RTP_API
GstBuffer*      gst_rtcp_buffer_new_take_data     (gpointer data, guint len);

GST_RTP_API
GstBuffer*      gst_rtcp_buffer_new_copy_data     (gconstpointer data, guint len);

GST_RTP_API
gboolean        gst_rtcp_buffer_validate_data     (guint8 *data, guint len);

GST_RTP_API
gboolean        gst_rtcp_buffer_validate          (GstBuffer *buffer);

GST_RTP_API
gboolean        gst_rtcp_buffer_validate_data_reduced   (guint8 *data, guint len);

GST_RTP_API
gboolean        gst_rtcp_buffer_validate_reduced        (GstBuffer *buffer);


GST_RTP_API
GstBuffer*      gst_rtcp_buffer_new               (guint mtu);

GST_RTP_API
gboolean        gst_rtcp_buffer_map               (GstBuffer *buffer, GstMapFlags flags, GstRTCPBuffer *rtcp);

GST_RTP_API
gboolean        gst_rtcp_buffer_unmap             (GstRTCPBuffer *rtcp);

/* adding/retrieving packets */

GST_RTP_API
guint           gst_rtcp_buffer_get_packet_count  (GstRTCPBuffer *rtcp);

GST_RTP_API
gboolean        gst_rtcp_buffer_get_first_packet  (GstRTCPBuffer *rtcp, GstRTCPPacket *packet);

GST_RTP_API
gboolean        gst_rtcp_packet_move_to_next      (GstRTCPPacket *packet);

GST_RTP_API
gboolean        gst_rtcp_buffer_add_packet        (GstRTCPBuffer *rtcp, GstRTCPType type,
                                                   GstRTCPPacket *packet);

GST_RTP_API
gboolean        gst_rtcp_packet_remove            (GstRTCPPacket *packet);

/* working with packets */

GST_RTP_API
gboolean        gst_rtcp_packet_get_padding       (GstRTCPPacket *packet);

GST_RTP_API
guint8          gst_rtcp_packet_get_count         (GstRTCPPacket *packet);

GST_RTP_API
GstRTCPType     gst_rtcp_packet_get_type          (GstRTCPPacket *packet);

GST_RTP_API
guint16         gst_rtcp_packet_get_length        (GstRTCPPacket *packet);


/* sender reports */

GST_RTP_API
void            gst_rtcp_packet_sr_get_sender_info    (GstRTCPPacket *packet, guint32 *ssrc,
                                                       guint64 *ntptime, guint32 *rtptime,
                                                       guint32 *packet_count, guint32 *octet_count);

GST_RTP_API
void            gst_rtcp_packet_sr_set_sender_info    (GstRTCPPacket *packet, guint32 ssrc,
                                                       guint64 ntptime, guint32 rtptime,
                                                       guint32 packet_count, guint32 octet_count);
/* receiver reports */

GST_RTP_API
guint32         gst_rtcp_packet_rr_get_ssrc           (GstRTCPPacket *packet);

GST_RTP_API
void            gst_rtcp_packet_rr_set_ssrc           (GstRTCPPacket *packet, guint32 ssrc);


/* report blocks for SR and RR */

GST_RTP_API
guint           gst_rtcp_packet_get_rb_count          (GstRTCPPacket *packet);

GST_RTP_API
void            gst_rtcp_packet_get_rb                (GstRTCPPacket *packet, guint nth, guint32 *ssrc,
                                                       guint8 *fractionlost, gint32 *packetslost,
                                                       guint32 *exthighestseq, guint32 *jitter,
                                                       guint32 *lsr, guint32 *dlsr);

GST_RTP_API
gboolean        gst_rtcp_packet_add_rb                (GstRTCPPacket *packet, guint32 ssrc,
                                                       guint8 fractionlost, gint32 packetslost,
                                                       guint32 exthighestseq, guint32 jitter,
                                                       guint32 lsr, guint32 dlsr);

GST_RTP_API
void            gst_rtcp_packet_set_rb                (GstRTCPPacket *packet, guint nth, guint32 ssrc,
                                                       guint8 fractionlost, gint32 packetslost,
                                                       guint32 exthighestseq, guint32 jitter,
                                                       guint32 lsr, guint32 dlsr);

/* profile-specific extensions for SR and RR */

GST_RTP_API
gboolean        gst_rtcp_packet_add_profile_specific_ext        (GstRTCPPacket * packet,
                                                                 const guint8 * data, guint len);

GST_RTP_API
guint16         gst_rtcp_packet_get_profile_specific_ext_length (GstRTCPPacket * packet);

GST_RTP_API
gboolean        gst_rtcp_packet_get_profile_specific_ext        (GstRTCPPacket * packet,
                                                                 guint8 ** data, guint * len);

GST_RTP_API
gboolean        gst_rtcp_packet_copy_profile_specific_ext       (GstRTCPPacket * packet,
                                                                 guint8 ** data, guint * len);

/* source description packet */

GST_RTP_API
guint           gst_rtcp_packet_sdes_get_item_count   (GstRTCPPacket *packet);

GST_RTP_API
gboolean        gst_rtcp_packet_sdes_first_item       (GstRTCPPacket *packet);

GST_RTP_API
gboolean        gst_rtcp_packet_sdes_next_item        (GstRTCPPacket *packet);

GST_RTP_API
guint32         gst_rtcp_packet_sdes_get_ssrc         (GstRTCPPacket *packet);

GST_RTP_API
gboolean        gst_rtcp_packet_sdes_first_entry      (GstRTCPPacket *packet);

GST_RTP_API
gboolean        gst_rtcp_packet_sdes_next_entry       (GstRTCPPacket *packet);

GST_RTP_API
gboolean        gst_rtcp_packet_sdes_get_entry        (GstRTCPPacket *packet,
                                                       GstRTCPSDESType *type, guint8 *len,
                                                       guint8 **data);

GST_RTP_API
gboolean        gst_rtcp_packet_sdes_copy_entry       (GstRTCPPacket *packet,
                                                       GstRTCPSDESType *type, guint8 *len,
                                                       guint8 **data);

GST_RTP_API
gboolean        gst_rtcp_packet_sdes_add_item         (GstRTCPPacket *packet, guint32 ssrc);

GST_RTP_API
gboolean        gst_rtcp_packet_sdes_add_entry        (GstRTCPPacket *packet, GstRTCPSDESType type,
                                                       guint8 len, const guint8 *data);

/* bye packet */

GST_RTP_API
guint           gst_rtcp_packet_bye_get_ssrc_count    (GstRTCPPacket *packet);

GST_RTP_API
guint32         gst_rtcp_packet_bye_get_nth_ssrc      (GstRTCPPacket *packet, guint nth);

GST_RTP_API
gboolean        gst_rtcp_packet_bye_add_ssrc          (GstRTCPPacket *packet, guint32 ssrc);

GST_RTP_API
gboolean        gst_rtcp_packet_bye_add_ssrcs         (GstRTCPPacket *packet, guint32 *ssrc, guint len);

GST_RTP_API
guint8          gst_rtcp_packet_bye_get_reason_len    (GstRTCPPacket *packet);

GST_RTP_API
gchar*          gst_rtcp_packet_bye_get_reason        (GstRTCPPacket *packet);

GST_RTP_API
gboolean        gst_rtcp_packet_bye_set_reason        (GstRTCPPacket *packet, const gchar *reason);

/* app packets */

GST_RTP_API
void            gst_rtcp_packet_app_set_subtype       (GstRTCPPacket * packet, guint8 subtype);

GST_RTP_API
guint8          gst_rtcp_packet_app_get_subtype       (GstRTCPPacket * packet);

GST_RTP_API
void            gst_rtcp_packet_app_set_ssrc          (GstRTCPPacket * packet, guint32 ssrc);

GST_RTP_API
guint32         gst_rtcp_packet_app_get_ssrc          (GstRTCPPacket * packet);

GST_RTP_API
void            gst_rtcp_packet_app_set_name          (GstRTCPPacket * packet, const gchar *name);

GST_RTP_API
const gchar*    gst_rtcp_packet_app_get_name          (GstRTCPPacket * packet);

GST_RTP_API
guint16         gst_rtcp_packet_app_get_data_length   (GstRTCPPacket * packet);

GST_RTP_API
gboolean        gst_rtcp_packet_app_set_data_length   (GstRTCPPacket * packet, guint16 wordlen);

GST_RTP_API
guint8*         gst_rtcp_packet_app_get_data          (GstRTCPPacket * packet);

/* feedback packets */

GST_RTP_API
guint32         gst_rtcp_packet_fb_get_sender_ssrc    (GstRTCPPacket *packet);

GST_RTP_API
void            gst_rtcp_packet_fb_set_sender_ssrc    (GstRTCPPacket *packet, guint32 ssrc);

GST_RTP_API
guint32         gst_rtcp_packet_fb_get_media_ssrc     (GstRTCPPacket *packet);

GST_RTP_API
void            gst_rtcp_packet_fb_set_media_ssrc     (GstRTCPPacket *packet, guint32 ssrc);

GST_RTP_API
GstRTCPFBType   gst_rtcp_packet_fb_get_type           (GstRTCPPacket *packet);

GST_RTP_API
void            gst_rtcp_packet_fb_set_type           (GstRTCPPacket *packet, GstRTCPFBType type);

GST_RTP_API
guint16         gst_rtcp_packet_fb_get_fci_length     (GstRTCPPacket *packet);

GST_RTP_API
gboolean        gst_rtcp_packet_fb_set_fci_length     (GstRTCPPacket *packet, guint16 wordlen);

GST_RTP_API
guint8 *        gst_rtcp_packet_fb_get_fci            (GstRTCPPacket *packet);

/* helper functions */

GST_RTP_API
guint64         gst_rtcp_ntp_to_unix                  (guint64 ntptime);

GST_RTP_API
guint64         gst_rtcp_unix_to_ntp                  (guint64 unixtime);

GST_RTP_API
const gchar *   gst_rtcp_sdes_type_to_name            (GstRTCPSDESType type);

GST_RTP_API
GstRTCPSDESType gst_rtcp_sdes_name_to_type            (const gchar *name);

/* extended report */

GST_RTP_API
guint32         gst_rtcp_packet_xr_get_ssrc           (GstRTCPPacket *packet);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_first_rb           (GstRTCPPacket *packet);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_next_rb            (GstRTCPPacket * packet);

GST_RTP_API
GstRTCPXRType   gst_rtcp_packet_xr_get_block_type     (GstRTCPPacket * packet);

GST_RTP_API
guint16         gst_rtcp_packet_xr_get_block_length   (GstRTCPPacket * packet);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_get_rle_info       (GstRTCPPacket * packet,
                                                       guint32 * ssrc, guint8 * thinning,
                                                       guint16 * begin_seq, guint16 * end_seq,
                                                       guint32 * chunk_count);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_get_rle_nth_chunk  (GstRTCPPacket * packet, guint nth,
                                                       guint16 * chunk);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_get_prt_info       (GstRTCPPacket * packet,
                                                       guint32 * ssrc, guint8 * thinning,
                                                       guint16 * begin_seq, guint16 * end_seq);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_get_prt_by_seq     (GstRTCPPacket * packet, guint16 seq,
                                                       guint32 * receipt_time);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_get_rrt            (GstRTCPPacket * packet, guint64 * timestamp);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_get_dlrr_block     (GstRTCPPacket * packet,
                                                       guint nth, guint32 * ssrc,
                                                       guint32 * last_rr, guint32 * delay);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_get_summary_info   (GstRTCPPacket * packet, guint32 * ssrc,
                                                       guint16 * begin_seq, guint16 * end_seq);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_get_summary_pkt    (GstRTCPPacket * packet,
                                                       guint32 * lost_packets, guint32 * dup_packets);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_get_summary_jitter (GstRTCPPacket * packet,
                                                       guint32 * min_jitter, guint32 * max_jitter,
                                                       guint32 * mean_jitter, guint32 * dev_jitter);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_get_summary_ttl    (GstRTCPPacket * packet, gboolean * is_ipv4,
                                                       guint8 * min_ttl, guint8 * max_ttl,
                                                       guint8 * mean_ttl, guint8 * dev_ttl);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_get_voip_metrics_ssrc        (GstRTCPPacket * packet, guint32 * ssrc);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_get_voip_packet_metrics      (GstRTCPPacket * packet,
                                                                 guint8 * loss_rate, guint8 * discard_rate);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_get_voip_burst_metrics       (GstRTCPPacket * packet,
                                                                 guint8 * burst_density, guint8 * gap_density,
                                                                 guint16 * burst_duration, guint16 * gap_duration);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_get_voip_delay_metrics       (GstRTCPPacket * packet,
                                                                 guint16 * roundtrip_delay,
                                                                 guint16 * end_system_delay);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_get_voip_signal_metrics      (GstRTCPPacket * packet,
                                                                 guint8 * signal_level, guint8 * noise_level,
                                                                 guint8 * rerl, guint8 * gmin);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_get_voip_quality_metrics     (GstRTCPPacket * packet,
                                                                 guint8 * r_factor, guint8 * ext_r_factor,
                                                                 guint8 * mos_lq, guint8 * mos_cq);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_get_voip_configuration_params        (GstRTCPPacket * packet,
                                                                         guint8 * gmin, guint8 * rx_config);

GST_RTP_API
gboolean        gst_rtcp_packet_xr_get_voip_jitter_buffer_params        (GstRTCPPacket * packet,
                                                                         guint16 * jb_nominal,
                                                                         guint16 * jb_maximum,
                                                                         guint16 * jb_abs_max);

G_END_DECLS

#endif /* __GST_RTCPBUFFER_H__ */

