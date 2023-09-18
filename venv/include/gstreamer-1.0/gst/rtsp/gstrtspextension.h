/* GStreamer RTSP Extension
 * Copyright (C) 2007 Wim Taymans <wim@fluendo.com>
 *
 * gstrtspextension.h: RTSP Extension interface.
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

#ifndef __GST_RTSP_EXTENSION_H__
#define __GST_RTSP_EXTENSION_H__

#include <gst/gst.h>

#include <gst/sdp/gstsdpmessage.h>
#include <gst/rtsp/gstrtsptransport.h>
#include <gst/rtsp/gstrtspmessage.h>
#include <gst/rtsp/gstrtspurl.h>

G_BEGIN_DECLS

#define GST_TYPE_RTSP_EXTENSION \
  (gst_rtsp_extension_get_type ())
#define GST_RTSP_EXTENSION(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_RTSP_EXTENSION, GstRTSPExtension))
#define GST_IS_RTSP_EXTENSION(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_RTSP_EXTENSION))
#define GST_RTSP_EXTENSION_GET_INTERFACE(inst) \
  (G_TYPE_INSTANCE_GET_INTERFACE ((inst), GST_TYPE_RTSP_EXTENSION, GstRTSPExtensionInterface))

/**
 * GstRTSPExtensionInterface:
 *
 * An interface representing RTSP extensions.
 */
typedef struct _GstRTSPExtension GstRTSPExtension;
typedef struct _GstRTSPExtensionInterface GstRTSPExtensionInterface;

struct _GstRTSPExtensionInterface {
  GTypeInterface parent;

  /* vfunctions */
  gboolean      (*detect_server)    (GstRTSPExtension *ext, GstRTSPMessage *resp);

  GstRTSPResult (*before_send)      (GstRTSPExtension *ext, GstRTSPMessage *req);
  GstRTSPResult (*after_send)       (GstRTSPExtension *ext, GstRTSPMessage *req, GstRTSPMessage *resp);

  GstRTSPResult (*parse_sdp)        (GstRTSPExtension *ext, GstSDPMessage *sdp, GstStructure *s);
  GstRTSPResult (*setup_media)      (GstRTSPExtension *ext, GstSDPMedia *media);

  gboolean      (*configure_stream) (GstRTSPExtension *ext, GstCaps *caps);

  GstRTSPResult (*get_transports)   (GstRTSPExtension *ext, GstRTSPLowerTrans protocols, gchar **transport);

  GstRTSPResult (*stream_select)    (GstRTSPExtension *ext, GstRTSPUrl *url);

  /* signals */
  GstRTSPResult (*send)             (GstRTSPExtension *ext, GstRTSPMessage *req, GstRTSPMessage *resp);

  /* ABI: more vfunctions added later */
  GstRTSPResult (*receive_request)  (GstRTSPExtension *ext, GstRTSPMessage *req);

  /*< private >*/
  gpointer                 _gst_reserved[GST_PADDING];
};

GST_RTSP_API
GType           gst_rtsp_extension_get_type          (void);

/* invoke vfunction on interface */

GST_RTSP_API
gboolean        gst_rtsp_extension_detect_server     (GstRTSPExtension *ext, GstRTSPMessage *resp);

GST_RTSP_API
GstRTSPResult   gst_rtsp_extension_before_send       (GstRTSPExtension *ext, GstRTSPMessage *req);

GST_RTSP_API
GstRTSPResult   gst_rtsp_extension_after_send        (GstRTSPExtension *ext, GstRTSPMessage *req,
                                                      GstRTSPMessage *resp);

GST_RTSP_API
GstRTSPResult   gst_rtsp_extension_parse_sdp         (GstRTSPExtension *ext, GstSDPMessage *sdp,
                                                      GstStructure *s);

GST_RTSP_API
GstRTSPResult   gst_rtsp_extension_setup_media       (GstRTSPExtension *ext, GstSDPMedia *media);

GST_RTSP_API
gboolean        gst_rtsp_extension_configure_stream  (GstRTSPExtension *ext, GstCaps *caps);

GST_RTSP_API
GstRTSPResult   gst_rtsp_extension_get_transports    (GstRTSPExtension *ext, GstRTSPLowerTrans protocols,
                                                      gchar **transport);

GST_RTSP_API
GstRTSPResult   gst_rtsp_extension_stream_select     (GstRTSPExtension *ext, GstRTSPUrl *url);

GST_RTSP_API
GstRTSPResult   gst_rtsp_extension_receive_request   (GstRTSPExtension *ext, GstRTSPMessage *req);

/* signal emission */

GST_RTSP_API
GstRTSPResult   gst_rtsp_extension_send              (GstRTSPExtension *ext, GstRTSPMessage *req,
                                                      GstRTSPMessage *resp);

G_END_DECLS

#endif /* __GST_RTSP_EXTENSION_H__ */
