/* GStreamer
 * Copyright (C) 1999,2000 Erik Walthinsen <omega@cse.ogi.edu>
 *                    2005 Wim Taymans <wim@fluendo.com>
 *                    2005 Andy Wingo <wingo@pobox.com>
 * Copyright (C) 2012 Collabora Ltd. <tim.muller@collabora.co.uk>
 *
 * gstnetclientclock.h: clock that synchronizes itself to a time provider over
 * the network
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


#ifndef __GST_NET_CLIENT_CLOCK_H__
#define __GST_NET_CLIENT_CLOCK_H__

#include <gst/gst.h>
#include <gst/gstsystemclock.h>
#include <gst/net/net-prelude.h>

G_BEGIN_DECLS

#define GST_TYPE_NET_CLIENT_CLOCK \
  (gst_net_client_clock_get_type())
#define GST_NET_CLIENT_CLOCK(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_NET_CLIENT_CLOCK,GstNetClientClock))
#define GST_NET_CLIENT_CLOCK_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_NET_CLIENT_CLOCK,GstNetClientClockClass))
#define GST_IS_NET_CLIENT_CLOCK(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_NET_CLIENT_CLOCK))
#define GST_IS_NET_CLIENT_CLOCK_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_NET_CLIENT_CLOCK))

typedef struct _GstNetClientClock GstNetClientClock;
typedef struct _GstNetClientClockClass GstNetClientClockClass;
typedef struct _GstNetClientClockPrivate GstNetClientClockPrivate;

/**
 * GstNetClientClock:
 *
 * Opaque #GstNetClientClock structure.
 */
struct _GstNetClientClock {
  GstSystemClock clock;

  /*< private >*/
  GstNetClientClockPrivate *priv;

  gpointer _gst_reserved[GST_PADDING];
};

struct _GstNetClientClockClass {
  GstSystemClockClass parent_class;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_NET_API
GType           gst_net_client_clock_get_type	(void);

GST_NET_API
GstClock*	gst_net_client_clock_new	(const gchar *name, const gchar *remote_address,
                                                 gint remote_port, GstClockTime base_time);

#define GST_TYPE_NTP_CLOCK \
  (gst_ntp_clock_get_type())
#define GST_NTP_CLOCK(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_NTP_CLOCK,GstNtpClock))
#define GST_NTP_CLOCK_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_NTP_CLOCK,GstNtpClockClass))
#define GST_IS_NTP_CLOCK(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_NTP_CLOCK))
#define GST_IS_NTP_CLOCK_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_NTP_CLOCK))

typedef struct _GstNetClientClock GstNtpClock;
typedef struct _GstNetClientClockClass GstNtpClockClass;

GST_NET_API
GType           gst_ntp_clock_get_type	        (void);

GST_NET_API
GstClock*	gst_ntp_clock_new	        (const gchar *name, const gchar *remote_address,
                                                 gint remote_port, GstClockTime base_time);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstNetClientClock, gst_object_unref)

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstNtpClock, gst_object_unref)

G_END_DECLS

#endif /* __GST_NET_CLIENT_CLOCK_H__ */
