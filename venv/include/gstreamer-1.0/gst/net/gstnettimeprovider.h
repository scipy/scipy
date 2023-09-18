/* GStreamer
 * Copyright (C) 2005 Andy Wingo <wingo@pobox.com>
 *               2006 Joni Valtanen <joni.valtanen@movial.fi>
 * Copyright (C) 2012 Collabora Ltd. <tim.muller@collabora.co.uk>
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


#ifndef __GST_NET_TIME_PROVIDER_H__
#define __GST_NET_TIME_PROVIDER_H__

#include <gst/gst.h>
#include <gst/net/net-prelude.h>

G_BEGIN_DECLS

#define GST_TYPE_NET_TIME_PROVIDER \
  (gst_net_time_provider_get_type())
#define GST_NET_TIME_PROVIDER(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_NET_TIME_PROVIDER,GstNetTimeProvider))
#define GST_NET_TIME_PROVIDER_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_NET_TIME_PROVIDER,GstNetTimeProviderClass))
#define GST_IS_NET_TIME_PROVIDER(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_NET_TIME_PROVIDER))
#define GST_IS_NET_TIME_PROVIDER_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_NET_TIME_PROVIDER))

typedef struct _GstNetTimeProvider GstNetTimeProvider;
typedef struct _GstNetTimeProviderClass GstNetTimeProviderClass;
typedef struct _GstNetTimeProviderPrivate GstNetTimeProviderPrivate;

/**
 * GstNetTimeProvider:
 *
 * Opaque #GstNetTimeProvider structure.
 */
struct _GstNetTimeProvider {
  GstObject parent;

  /*< private >*/
  GstNetTimeProviderPrivate *priv;

  gpointer _gst_reserved[GST_PADDING];
};

struct _GstNetTimeProviderClass {
  GstObjectClass parent_class;

  gpointer _gst_reserved[GST_PADDING];
};

GST_NET_API
GType                   gst_net_time_provider_get_type  (void);

GST_NET_API
GstNetTimeProvider*     gst_net_time_provider_new       (GstClock *clock,
                                                         const gchar *address,
                                                         gint port);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstNetTimeProvider, gst_object_unref)

G_END_DECLS


#endif /* __GST_NET_TIME_PROVIDER_H__ */
