/* GStreamer Color Balance
 * Copyright (C) 2003 Ronald Bultje <rbultje@ronald.bitfreak.net>
 *
 * color-balance.h: image color balance interface design
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

#ifndef __GST_COLOR_BALANCE_H__
#define __GST_COLOR_BALANCE_H__

#include <gst/gst.h>
#include <gst/video/colorbalancechannel.h>

G_BEGIN_DECLS

#define GST_TYPE_COLOR_BALANCE \
  (gst_color_balance_get_type ())
#define GST_COLOR_BALANCE(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_COLOR_BALANCE, GstColorBalance))
#define GST_IS_COLOR_BALANCE(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_COLOR_BALANCE))
#define GST_COLOR_BALANCE_GET_INTERFACE(inst) \
  (G_TYPE_INSTANCE_GET_INTERFACE ((inst), GST_TYPE_COLOR_BALANCE, GstColorBalanceInterface))

typedef struct _GstColorBalance GstColorBalance;
typedef struct _GstColorBalanceInterface GstColorBalanceInterface;

/**
 * GstColorBalanceType:
 * @GST_COLOR_BALANCE_HARDWARE: Color balance is implemented with dedicated
 *         hardware.
 * @GST_COLOR_BALANCE_SOFTWARE: Color balance is implemented via software
 *         processing.
 *
 * An enumeration indicating whether an element implements color balancing
 * operations in software or in dedicated hardware. In general, dedicated
 * hardware implementations (such as those provided by xvimagesink) are
 * preferred.
 */
typedef enum
{
  GST_COLOR_BALANCE_HARDWARE,
  GST_COLOR_BALANCE_SOFTWARE
} GstColorBalanceType;

/**
 * GstColorBalanceInterface:
 * @iface: the parent interface
 * @get_balance_type: implementation type
 * @list_channels: list handled channels
 * @set_value: set a channel value
 * @get_value: get a channel value
 * @value_changed: default handler for value changed notification
 *
 * Color-balance interface.
 */
struct _GstColorBalanceInterface {
  GTypeInterface iface;

  /* virtual functions */
  const GList * (* list_channels) (GstColorBalance        *balance);

  void          (* set_value)     (GstColorBalance        *balance,
                                   GstColorBalanceChannel *channel,
                                   gint                    value);
  gint          (* get_value)     (GstColorBalance        *balance,
                                   GstColorBalanceChannel *channel);
  GstColorBalanceType (*get_balance_type)  (GstColorBalance *balance);

  /* signals */
  void (* value_changed) (GstColorBalance        *balance,
                          GstColorBalanceChannel *channel,
                          gint                    value);

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_VIDEO_API
GType   gst_color_balance_get_type      (void);

/* virtual class function wrappers */

GST_VIDEO_API
const GList *
        gst_color_balance_list_channels (GstColorBalance        *balance);

GST_VIDEO_API
void    gst_color_balance_set_value     (GstColorBalance        *balance,
                                         GstColorBalanceChannel *channel,
                                         gint                    value);

GST_VIDEO_API
gint    gst_color_balance_get_value     (GstColorBalance        *balance,
                                         GstColorBalanceChannel *channel);

GST_VIDEO_API
GstColorBalanceType
        gst_color_balance_get_balance_type (GstColorBalance        *balance);

/* trigger signal */

GST_VIDEO_API
void    gst_color_balance_value_changed (GstColorBalance        *balance,
                                         GstColorBalanceChannel *channel,
                                         gint                    value);

G_END_DECLS

#endif /* __GST_COLOR_BALANCE_H__ */
