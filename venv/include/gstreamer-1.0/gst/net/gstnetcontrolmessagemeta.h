/* GStreamer
 * Copyright (C) <2014> William Manley <will@williammanley.net>
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

#ifndef __GST_NET_CONTROL_MESSAGE_META_H__
#define __GST_NET_CONTROL_MESSAGE_META_H__

#include <gst/gst.h>
#include <gio/gio.h>
#include <gst/net/net-prelude.h>

G_BEGIN_DECLS

typedef struct _GstNetControlMessageMeta GstNetControlMessageMeta;

/**
 * GstNetControlMessageMeta:
 * @meta: the parent type
 * @message: a #GSocketControlMessage stored as metadata
 *
 * Buffer metadata for GSocket control messages, AKA ancillary data attached to
 * data sent across a socket.
 */
struct _GstNetControlMessageMeta {
  GstMeta       meta;

  GSocketControlMessage *message;
};

GST_NET_API
GType gst_net_control_message_meta_api_get_type (void);

#define GST_NET_CONTROL_MESSAGE_META_API_TYPE \
  (gst_net_control_message_meta_api_get_type())

#define gst_buffer_get_net_control_message_meta(b) ((GstNetControlMessageMeta*)\
  gst_buffer_get_meta((b),GST_NET_CONTROL_MESSAGE_META_API_TYPE))

/* implementation */

GST_NET_API
const GstMetaInfo *gst_net_control_message_meta_get_info (void);

#define GST_NET_CONTROL_MESSAGE_META_INFO \
  (gst_net_control_message_meta_get_info())

GST_NET_API
GstNetControlMessageMeta * gst_buffer_add_net_control_message_meta (GstBuffer             * buffer,
                                                                    GSocketControlMessage * message);

G_END_DECLS

#endif /* __GST_NET_CONTROL_MESSAGE_META_H__ */

