/* GStreamer
 * Copyright (C) <2007> Wim Taymans <wim.taymans@gmail.com>
 *
 * gstrtppayloads.h: various helper functions to deal with RTP payload
 *     types.
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

#ifndef __GST_RTPPAYLOADS_H__
#define __GST_RTPPAYLOADS_H__

#include <gst/gst.h>
#include <gst/rtp/rtp-prelude.h>

G_BEGIN_DECLS

/**
 * GstRTPPayload:
 * @GST_RTP_PAYLOAD_PCMU: ITU-T G.711. mu-law audio (RFC 3551)
 * @GST_RTP_PAYLOAD_1016: RFC 3551 says reserved
 * @GST_RTP_PAYLOAD_G721: RFC 3551 says reserved
 * @GST_RTP_PAYLOAD_GSM: GSM audio
 * @GST_RTP_PAYLOAD_G723: ITU G.723.1 audio
 * @GST_RTP_PAYLOAD_DVI4_8000: IMA ADPCM wave type (RFC 3551)
 * @GST_RTP_PAYLOAD_DVI4_16000: IMA ADPCM wave type (RFC 3551)
 * @GST_RTP_PAYLOAD_LPC: experimental linear predictive encoding
 * @GST_RTP_PAYLOAD_PCMA: ITU-T G.711 A-law audio (RFC 3551)
 * @GST_RTP_PAYLOAD_G722: ITU-T G.722 (RFC 3551)
 * @GST_RTP_PAYLOAD_L16_STEREO: stereo PCM
 * @GST_RTP_PAYLOAD_L16_MONO: mono PCM
 * @GST_RTP_PAYLOAD_QCELP: EIA & TIA standard IS-733
 * @GST_RTP_PAYLOAD_CN: Comfort Noise (RFC 3389)
 * @GST_RTP_PAYLOAD_MPA: Audio MPEG 1-3.
 * @GST_RTP_PAYLOAD_G728: ITU-T G.728 Speech coder (RFC 3551)
 * @GST_RTP_PAYLOAD_DVI4_11025: IMA ADPCM wave type (RFC 3551)
 * @GST_RTP_PAYLOAD_DVI4_22050: IMA ADPCM wave type (RFC 3551)
 * @GST_RTP_PAYLOAD_G729: ITU-T G.729 Speech coder (RFC 3551)
 * @GST_RTP_PAYLOAD_CELLB: See RFC 2029
 * @GST_RTP_PAYLOAD_JPEG: ISO Standards 10918-1 and 10918-2 (RFC 2435)
 * @GST_RTP_PAYLOAD_NV: nv encoding by Ron Frederick
 * @GST_RTP_PAYLOAD_H261: ITU-T Recommendation H.261 (RFC 2032)
 * @GST_RTP_PAYLOAD_MPV: Video MPEG 1 & 2 (RFC 2250)
 * @GST_RTP_PAYLOAD_MP2T: MPEG-2 transport stream (RFC 2250)
 * @GST_RTP_PAYLOAD_H263: Video H263 (RFC 2190)
 *
 * Standard predefined fixed payload types.
 *
 * The official list is at:
 * http://www.iana.org/assignments/rtp-parameters
 *
 * Audio:
 * reserved: 19
 * unassigned: 20-23,
 *
 * Video:
 * unassigned: 24, 27, 29, 30, 35-71, 77-95
 * Reserved for RTCP conflict avoidance: 72-76
 */
typedef enum
{
  /* Audio: */
  GST_RTP_PAYLOAD_PCMU = 0,
  GST_RTP_PAYLOAD_1016 = 1, /* RFC 3551 says reserved */
  GST_RTP_PAYLOAD_G721 = 2, /* RFC 3551 says reserved */
  GST_RTP_PAYLOAD_GSM = 3,
  GST_RTP_PAYLOAD_G723 = 4,
  GST_RTP_PAYLOAD_DVI4_8000 = 5,
  GST_RTP_PAYLOAD_DVI4_16000 = 6,
  GST_RTP_PAYLOAD_LPC = 7,
  GST_RTP_PAYLOAD_PCMA = 8,
  GST_RTP_PAYLOAD_G722 = 9,
  GST_RTP_PAYLOAD_L16_STEREO = 10,
  GST_RTP_PAYLOAD_L16_MONO = 11,
  GST_RTP_PAYLOAD_QCELP = 12,
  GST_RTP_PAYLOAD_CN = 13,
  GST_RTP_PAYLOAD_MPA = 14,
  GST_RTP_PAYLOAD_G728 = 15,
  GST_RTP_PAYLOAD_DVI4_11025 = 16,
  GST_RTP_PAYLOAD_DVI4_22050 = 17,
  GST_RTP_PAYLOAD_G729 = 18,

  /* Video: */

  GST_RTP_PAYLOAD_CELLB = 25,
  GST_RTP_PAYLOAD_JPEG = 26,
  GST_RTP_PAYLOAD_NV = 28,
  GST_RTP_PAYLOAD_H261 = 31,
  GST_RTP_PAYLOAD_MPV = 32,
  GST_RTP_PAYLOAD_MP2T = 33,
  GST_RTP_PAYLOAD_H263 = 34,

  /* BOTH */
} GstRTPPayload;

/* backward compatibility */
#define GST_RTP_PAYLOAD_G723_63 16
#define GST_RTP_PAYLOAD_G723_53 17
#define GST_RTP_PAYLOAD_TS48 18
#define GST_RTP_PAYLOAD_TS41 19

#define GST_RTP_PAYLOAD_G723_63_STRING "16"
#define GST_RTP_PAYLOAD_G723_53_STRING "17"
#define GST_RTP_PAYLOAD_TS48_STRING "18"
#define GST_RTP_PAYLOAD_TS41_STRING "19"

/* Defining the above as strings, to make the declaration of pad_templates
 * easier. So if please keep these synchronized with the above.
 */
#define GST_RTP_PAYLOAD_PCMU_STRING "0"
#define GST_RTP_PAYLOAD_1016_STRING "1"
#define GST_RTP_PAYLOAD_G721_STRING "2"
#define GST_RTP_PAYLOAD_GSM_STRING "3"
#define GST_RTP_PAYLOAD_G723_STRING "4"
#define GST_RTP_PAYLOAD_DVI4_8000_STRING "5"
#define GST_RTP_PAYLOAD_DVI4_16000_STRING "6"
#define GST_RTP_PAYLOAD_LPC_STRING "7"
#define GST_RTP_PAYLOAD_PCMA_STRING "8"
#define GST_RTP_PAYLOAD_G722_STRING "9"
#define GST_RTP_PAYLOAD_L16_STEREO_STRING "10"
#define GST_RTP_PAYLOAD_L16_MONO_STRING "11"
#define GST_RTP_PAYLOAD_QCELP_STRING "12"
#define GST_RTP_PAYLOAD_CN_STRING "13"
#define GST_RTP_PAYLOAD_MPA_STRING "14"
#define GST_RTP_PAYLOAD_G728_STRING "15"
#define GST_RTP_PAYLOAD_DVI4_11025_STRING "16"
#define GST_RTP_PAYLOAD_DVI4_22050_STRING "17"
#define GST_RTP_PAYLOAD_G729_STRING "18"

#define GST_RTP_PAYLOAD_CELLB_STRING "25"
#define GST_RTP_PAYLOAD_JPEG_STRING "26"
#define GST_RTP_PAYLOAD_NV_STRING "28"

#define GST_RTP_PAYLOAD_H261_STRING "31"
#define GST_RTP_PAYLOAD_MPV_STRING "32"
#define GST_RTP_PAYLOAD_MP2T_STRING "33"
#define GST_RTP_PAYLOAD_H263_STRING "34"

#define GST_RTP_PAYLOAD_DYNAMIC_STRING "[96, 127]"

/**
 * GST_RTP_PAYLOAD_IS_DYNAMIC:
 * @pt: a payload type
 *
 * Check if @pt is a dynamic payload type.
 */
#define GST_RTP_PAYLOAD_IS_DYNAMIC(pt) ((pt) >= 96 && (pt) <= 127)

typedef struct _GstRTPPayloadInfo GstRTPPayloadInfo;

/**
 * GstRTPPayloadInfo:
 * @payload_type: payload type, -1 means dynamic
 * @media: the media type(s), usually "audio", "video", "application", "text",
 * "message".
 * @encoding_name: the encoding name of @pt
 * @clock_rate: default clock rate, 0 = unknown/variable
 * @encoding_parameters: encoding parameters. For audio this is the number of
 * channels. NULL = not applicable.
 * @bitrate: the bitrate of the media. 0 = unknown/variable.
 *
 * Structure holding default payload type information.
 */
struct _GstRTPPayloadInfo
{
  guint8       payload_type;
  const gchar *media;
  const gchar *encoding_name;
  guint        clock_rate;
  const gchar *encoding_parameters;
  guint        bitrate;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_RTP_API
const GstRTPPayloadInfo *  gst_rtp_payload_info_for_pt     (guint8 payload_type);

GST_RTP_API
const GstRTPPayloadInfo *  gst_rtp_payload_info_for_name   (const gchar *media, const gchar *encoding_name);

G_END_DECLS

#endif /* __GST_RTPPAYLOADS_H__ */

