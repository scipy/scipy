/* GStreamer Color Balance
 * Copyright (C) 2003 Ronald Bultje <rbultje@ronald.bitfreak.net>
 *
 * colorbalancechannel.h: individual channel object
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

#ifndef __GST_COLOR_BALANCE_CHANNEL_H__
#define __GST_COLOR_BALANCE_CHANNEL_H__

#include <gst/gst.h>
#include <gst/video/video-prelude.h>

G_BEGIN_DECLS

#define GST_TYPE_COLOR_BALANCE_CHANNEL \
  (gst_color_balance_channel_get_type ())
#define GST_COLOR_BALANCE_CHANNEL(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_COLOR_BALANCE_CHANNEL, \
                               GstColorBalanceChannel))
#define GST_COLOR_BALANCE_CHANNEL_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST ((klass), GST_TYPE_COLOR_BALANCE_CHANNEL, \
                            GstColorBalanceChannelClass))
#define GST_IS_COLOR_BALANCE_CHANNEL(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_COLOR_BALANCE_CHANNEL))
#define GST_IS_COLOR_BALANCE_CHANNEL_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_TYPE ((klass), GST_TYPE_COLOR_BALANCE_CHANNEL))

typedef struct _GstColorBalanceChannel GstColorBalanceChannel;
typedef struct _GstColorBalanceChannelClass GstColorBalanceChannelClass;

/**
 * GstColorBalanceChannel:
 * @label: A string containing a descriptive name for this channel
 * @min_value: The minimum valid value for this channel.
 * @max_value: The maximum valid value for this channel.
 */
struct _GstColorBalanceChannel {
  GObject parent;

  /*< public >*/
  gchar  *label;
  gint    min_value;
  gint    max_value;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstColorBalanceChannelClass:
 * @parent: the parent class
 * @value_changed: default handler for value changed notification
 *
 * Color-balance channel class.
 */
struct _GstColorBalanceChannelClass {
  GObjectClass parent;

  /* signals */
  void (* value_changed) (GstColorBalanceChannel *channel,
                          gint                    value);

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_VIDEO_API
GType   gst_color_balance_channel_get_type (void);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstColorBalanceChannel, gst_object_unref)

G_END_DECLS

#endif /* __GST_COLOR_BALANCE_CHANNEL_H__ */
