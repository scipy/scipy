/* GStreamer
 * Copyright (C) <2005,2006> Wim Taymans <wim@fluendo.com>
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
/*
 * Unless otherwise indicated, Source Code is licensed under MIT license.
 * See further explanation attached in License Statement (distributed in the file
 * LICENSE).
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to do
 * so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef __GST_RTSP_TRANSPORT_H__
#define __GST_RTSP_TRANSPORT_H__

#include <gst/gstconfig.h>
#include <gst/rtsp/gstrtspdefs.h>
#include <gst/rtsp/gstrtsp-enumtypes.h>

G_BEGIN_DECLS

/**
 * GstRTSPTransMode:
 * @GST_RTSP_TRANS_UNKNOWN: invalid tansport mode
 * @GST_RTSP_TRANS_RTP: transfer RTP data
 * @GST_RTSP_TRANS_RDT: transfer RDT (RealMedia) data
 *
 * The transfer mode to use.
 */
typedef enum {
  GST_RTSP_TRANS_UNKNOWN =  0,
  GST_RTSP_TRANS_RTP     = (1 << 0),
  GST_RTSP_TRANS_RDT     = (1 << 1)
} GstRTSPTransMode;

/**
 * GstRTSPProfile:
 * @GST_RTSP_PROFILE_UNKNOWN: invalid profile
 * @GST_RTSP_PROFILE_AVP: the Audio/Visual profile (RFC 3551)
 * @GST_RTSP_PROFILE_SAVP: the secure Audio/Visual profile (RFC 3711)
 * @GST_RTSP_PROFILE_AVPF: the Audio/Visual profile with feedback (RFC 4585)
 * @GST_RTSP_PROFILE_SAVPF: the secure Audio/Visual profile with feedback (RFC 5124)
 *
 * The transfer profile to use.
 */
/* FIXME 2.0: This should probably be an enum, not flags and maybe be replaced
 * by GstRTPTransport */
typedef enum {
  GST_RTSP_PROFILE_UNKNOWN =  0,
  GST_RTSP_PROFILE_AVP     = (1 << 0),
  GST_RTSP_PROFILE_SAVP    = (1 << 1),
  GST_RTSP_PROFILE_AVPF    = (1 << 2),
  GST_RTSP_PROFILE_SAVPF   = (1 << 3),
} GstRTSPProfile;

/**
 * GstRTSPLowerTrans:
 * @GST_RTSP_LOWER_TRANS_UNKNOWN: invalid transport flag
 * @GST_RTSP_LOWER_TRANS_UDP: stream data over UDP
 * @GST_RTSP_LOWER_TRANS_UDP_MCAST: stream data over UDP multicast
 * @GST_RTSP_LOWER_TRANS_TCP: stream data over TCP
 * @GST_RTSP_LOWER_TRANS_HTTP: stream data tunneled over HTTP.
 * @GST_RTSP_LOWER_TRANS_TLS: encrypt TCP and HTTP with TLS
 *
 * The different transport methods.
 */
typedef enum {
  GST_RTSP_LOWER_TRANS_UNKNOWN   = 0,
  GST_RTSP_LOWER_TRANS_UDP       = (1 << 0),
  GST_RTSP_LOWER_TRANS_UDP_MCAST = (1 << 1),
  GST_RTSP_LOWER_TRANS_TCP       = (1 << 2),
  GST_RTSP_LOWER_TRANS_HTTP      = (1 << 4),
  GST_RTSP_LOWER_TRANS_TLS       = (1 << 5)
} GstRTSPLowerTrans;

typedef struct _GstRTSPRange GstRTSPRange;
typedef struct _GstRTSPTransport GstRTSPTransport;

/**
 * GstRTSPRange:
 * @min: minimum value of the range
 * @max: maximum value of the range
 *
 * A type to specify a range.
 */

struct _GstRTSPRange {
  gint min;
  gint max;
};

/**
 * GstRTSPTransport:
 * @trans: the transport mode
 * @profile: the tansport profile
 * @lower_transport: the lower transport
 * @destination: the destination ip/hostname
 * @source: the source ip/hostname
 * @layers: the number of layers
 * @mode_play: if play mode was selected
 * @mode_record: if record mode was selected
 * @append: is append mode was selected
 * @interleaved: the interleave range
 * @ttl: the time to live for multicast UDP
 * @port: the port pair for multicast sessions
 * @client_port: the client port pair for receiving data. For TCP
 *   based transports, applications can use this field to store the
 *   sender and receiver ports of the client.
 * @server_port: the server port pair for receiving data. For TCP
 *   based transports, applications can use this field to store the
 *   sender and receiver ports of the server.
 * @ssrc: the ssrc that the sender/receiver will use
 *
 * A structure holding the RTSP transport values.
 */

struct _GstRTSPTransport {
  GstRTSPTransMode  trans;
  GstRTSPProfile    profile;
  GstRTSPLowerTrans lower_transport;

  gchar         *destination;
  gchar         *source;
  guint          layers;
  gboolean       mode_play;
  gboolean       mode_record;
  gboolean       append;
  GstRTSPRange   interleaved;

  /* multicast specific */
  guint  ttl;
  GstRTSPRange   port;

  /* UDP/TCP specific */
  GstRTSPRange   client_port;
  GstRTSPRange   server_port;
  /* RTP specific */
  guint          ssrc;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_RTSP_API
GstRTSPResult      gst_rtsp_transport_new          (GstRTSPTransport **transport);

GST_RTSP_API
GstRTSPResult      gst_rtsp_transport_init         (GstRTSPTransport *transport);

GST_RTSP_API
GstRTSPResult      gst_rtsp_transport_parse        (const gchar *str, GstRTSPTransport *transport);

GST_RTSP_API
gchar*             gst_rtsp_transport_as_text      (GstRTSPTransport *transport);

GST_RTSP_DEPRECATED_FOR(gst_rtsp_transport_get_media_type)
GstRTSPResult      gst_rtsp_transport_get_mime     (GstRTSPTransMode trans, const gchar **mime);

GST_RTSP_API
GstRTSPResult      gst_rtsp_transport_get_manager  (GstRTSPTransMode trans, const gchar **manager, guint option);

GST_RTSP_API
GstRTSPResult      gst_rtsp_transport_get_media_type (GstRTSPTransport *transport,
                                                      const gchar **media_type);

GST_RTSP_API
GstRTSPResult      gst_rtsp_transport_free         (GstRTSPTransport *transport);

G_END_DECLS

#endif /* __GST_RTSP_TRANSPORT_H__ */
