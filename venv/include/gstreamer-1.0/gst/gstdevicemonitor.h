/* GStreamer
 * Copyright (C) 2013 Olivier Crete <olivier.crete@collabora.com>
 *
 * gstdevicemonitor.c: Device monitor
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
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#ifndef __GST_DEVICE_MONITOR_H__
#define __GST_DEVICE_MONITOR_H__

#include <gst/gstobject.h>
#include <gst/gstdevice.h>
#include <gst/gstdeviceprovider.h>
#include <gst/gstdeviceproviderfactory.h>

G_BEGIN_DECLS

typedef struct _GstDeviceMonitor GstDeviceMonitor;
typedef struct _GstDeviceMonitorPrivate GstDeviceMonitorPrivate;
typedef struct _GstDeviceMonitorClass GstDeviceMonitorClass;

#define GST_TYPE_DEVICE_MONITOR                 (gst_device_monitor_get_type())
#define GST_IS_DEVICE_MONITOR(obj)              (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_DEVICE_MONITOR))
#define GST_IS_DEVICE_MONITOR_CLASS(klass)      (G_TYPE_CHECK_CLASS_TYPE ((klass), GST_TYPE_DEVICE_MONITOR))
#define GST_DEVICE_MONITOR_GET_CLASS(obj)       (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_DEVICE_MONITOR, GstDeviceMonitorClass))
#define GST_DEVICE_MONITOR(obj)                 (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_DEVICE_MONITOR, GstDeviceMonitor))
#define GST_DEVICE_MONITOR_CLASS(klass)         (G_TYPE_CHECK_CLASS_CAST ((klass), GST_TYPE_DEVICE_MONITOR, GstDeviceMonitorClass))
#define GST_DEVICE_MONITOR_CAST(obj)            ((GstDeviceMonitor *)(obj))

/**
 * GstDeviceMonitor:
 * @parent: the parent #GstObject structure
 *
 * Opaque device monitor object structure.
 *
 * Since: 1.4
 */
struct _GstDeviceMonitor {
  GstObject                parent;

  /*< private >*/

  GstDeviceMonitorPrivate *priv;

  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstDeviceMonitorClass:
 * @parent_class: the parent #GstObjectClass structure
 *
 * Opaque device monitor class structure.
 *
 * Since: 1.4
 */
struct _GstDeviceMonitorClass {
  GstObjectClass           parent_class;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_API
GType     gst_device_monitor_get_type (void);

GST_API
GstDeviceMonitor * gst_device_monitor_new  (void);

GST_API
GstBus *  gst_device_monitor_get_bus (GstDeviceMonitor * monitor);

GST_API
GList *   gst_device_monitor_get_devices (GstDeviceMonitor * monitor);


GST_API
gboolean  gst_device_monitor_start (GstDeviceMonitor * monitor);

GST_API
void      gst_device_monitor_stop  (GstDeviceMonitor * monitor);


GST_API
guint     gst_device_monitor_add_filter (GstDeviceMonitor * monitor,
                                         const gchar      * classes,
                                         GstCaps          * caps);
GST_API
gboolean  gst_device_monitor_remove_filter (GstDeviceMonitor * monitor,
                                            guint filter_id);
GST_API
gchar **  gst_device_monitor_get_providers (GstDeviceMonitor * monitor);

GST_API
void      gst_device_monitor_set_show_all_devices (GstDeviceMonitor * monitor, gboolean show_all);

GST_API
gboolean  gst_device_monitor_get_show_all_devices (GstDeviceMonitor * monitor);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstDeviceMonitor, gst_object_unref)

G_END_DECLS

#endif /* __GST_DEVICE_MONITOR_H__ */
