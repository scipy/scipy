/* GStreamer
 * Copyright (C) 2004 Wim Taymans <wim@fluendo.com>
 *
 * gstbus.h: Header for GstBus subsystem
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

#ifndef __GST_BUS_H__
#define __GST_BUS_H__

typedef struct _GstBus GstBus;
typedef struct _GstBusPrivate GstBusPrivate;
typedef struct _GstBusClass GstBusClass;

#include <gst/gstmessage.h>
#include <gst/gstclock.h>

G_BEGIN_DECLS

/* --- standard type macros --- */
#define GST_TYPE_BUS              (gst_bus_get_type ())
#define GST_BUS(bus)              (G_TYPE_CHECK_INSTANCE_CAST ((bus), GST_TYPE_BUS, GstBus))
#define GST_IS_BUS(bus)           (G_TYPE_CHECK_INSTANCE_TYPE ((bus), GST_TYPE_BUS))
#define GST_BUS_CLASS(bclass)     (G_TYPE_CHECK_CLASS_CAST ((bclass), GST_TYPE_BUS, GstBusClass))
#define GST_IS_BUS_CLASS(bclass)  (G_TYPE_CHECK_CLASS_TYPE ((bclass), GST_TYPE_BUS))
#define GST_BUS_GET_CLASS(bus)    (G_TYPE_INSTANCE_GET_CLASS ((bus), GST_TYPE_BUS, GstBusClass))
#define GST_BUS_CAST(bus)         ((GstBus*)(bus))

/**
 * GstBusFlags:
 * @GST_BUS_FLUSHING: The bus is currently dropping all messages
 * @GST_BUS_FLAG_LAST: offset to define more flags
 *
 * The standard flags that a bus may have.
 */
typedef enum {
  GST_BUS_FLUSHING      = (GST_OBJECT_FLAG_LAST << 0),
  /* padding */
  GST_BUS_FLAG_LAST     = (GST_OBJECT_FLAG_LAST << 1)
} GstBusFlags;

/**
 * GstBusSyncReply:
 * @GST_BUS_DROP: drop the message
 * @GST_BUS_PASS: pass the message to the async queue
 * @GST_BUS_ASYNC: pass message to async queue, continue if message is handled
 *
 * The result values for a GstBusSyncHandler.
 */
typedef enum
{
  GST_BUS_DROP = 0,
  GST_BUS_PASS = 1,
  GST_BUS_ASYNC = 2
} GstBusSyncReply;

/**
 * GstBusSyncHandler:
 * @bus: the #GstBus that sent the message
 * @message: the #GstMessage
 * @user_data: user data that has been given, when registering the handler
 *
 * Handler will be invoked synchronously, when a new message has been injected
 * into the bus. This function is mostly used internally. Only one sync handler
 * can be attached to a given bus.
 *
 * If the handler returns %GST_BUS_DROP, it should unref the message, else the
 * message should not be unreffed by the sync handler.
 *
 * Returns: #GstBusSyncReply stating what to do with the message
 */
typedef GstBusSyncReply (*GstBusSyncHandler)    (GstBus * bus, GstMessage * message, gpointer user_data);

/**
 * GstBusFunc:
 * @bus: the #GstBus that sent the message
 * @message: the #GstMessage
 * @user_data: user data that has been given, when registering the handler
 *
 * Specifies the type of function passed to gst_bus_add_watch() or
 * gst_bus_add_watch_full(), which is called from the mainloop when a message
 * is available on the bus.
 *
 * The message passed to the function will be unreffed after execution of this
 * function so it should not be freed in the function.
 *
 * Note that this function is used as a #GSourceFunc which means that returning
 * %FALSE will remove the #GSource from the mainloop.
 *
 * Returns: %FALSE if the event source should be removed.
 */
typedef gboolean        (*GstBusFunc)           (GstBus * bus, GstMessage * message, gpointer user_data);

/**
 * GstBus:
 * @object: the parent structure
 *
 * The opaque #GstBus data structure.
 */
struct _GstBus
{
  GstObject         object;

  /*< private >*/
  GstBusPrivate    *priv;

  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstBusClass:
 * @parent_class: the parent class structure
 *
 * GStreamer bus class.
 */
struct _GstBusClass
{
  GstObjectClass parent_class;

  /* signals */

  /**
   * GstBusClass::message:
   * @bus: the #GstBus
   * @message: the message that has been posted asynchronously
   *
   * A message has been posted on the bus.
   */
  void (*message)       (GstBus *bus, GstMessage *message);

  /**
   * GstBusClass::sync_message:
   * @bus: the #GstBus
   * @message: the message that has been posted synchronously
   *
   * A message has been posted on the bus.
   */
  void (*sync_message)  (GstBus *bus, GstMessage *message);

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_API
GType                   gst_bus_get_type                (void);

GST_API
GstBus*                 gst_bus_new                     (void);

GST_API
gboolean                gst_bus_post                    (GstBus * bus, GstMessage * message);

GST_API
gboolean                gst_bus_have_pending            (GstBus * bus);

GST_API
GstMessage *            gst_bus_peek                    (GstBus * bus);

GST_API
GstMessage *            gst_bus_pop                     (GstBus * bus);

GST_API
GstMessage *            gst_bus_pop_filtered            (GstBus * bus, GstMessageType types);

GST_API
GstMessage *            gst_bus_timed_pop               (GstBus * bus, GstClockTime timeout);

GST_API
GstMessage *            gst_bus_timed_pop_filtered      (GstBus * bus, GstClockTime timeout, GstMessageType types);

GST_API
void                    gst_bus_set_flushing            (GstBus * bus, gboolean flushing);

/* synchronous dispatching */

GST_API
void                    gst_bus_set_sync_handler        (GstBus * bus, GstBusSyncHandler func,
                                                         gpointer user_data, GDestroyNotify notify);

/* asynchronous message notifications */

GST_API
void                    gst_bus_get_pollfd              (GstBus * bus, GPollFD *fd);

/* GSource based dispatching */

GST_API
GSource *               gst_bus_create_watch            (GstBus * bus);

GST_API
guint                   gst_bus_add_watch_full          (GstBus * bus,
                                                         gint priority,
                                                         GstBusFunc func,
                                                         gpointer user_data,
                                                         GDestroyNotify notify);
GST_API
guint                   gst_bus_add_watch               (GstBus * bus,
                                                         GstBusFunc func,
                                                         gpointer user_data);
GST_API
gboolean                gst_bus_remove_watch            (GstBus * bus);

/* polling the bus */

GST_API
GstMessage*             gst_bus_poll                    (GstBus *bus, GstMessageType events,
                                                         GstClockTime timeout);

/* signal based dispatching helper functions. */

GST_API
gboolean                gst_bus_async_signal_func       (GstBus *bus, GstMessage *message,
                                                         gpointer data);
GST_API
GstBusSyncReply         gst_bus_sync_signal_handler     (GstBus *bus, GstMessage *message,
                                                         gpointer data);

/* convenience api to add/remove a gsource that emits the async signals */

GST_API
void                    gst_bus_add_signal_watch        (GstBus * bus);

GST_API
void                    gst_bus_add_signal_watch_full   (GstBus * bus, gint priority);

GST_API
void                    gst_bus_remove_signal_watch     (GstBus * bus);

GST_API
void                    gst_bus_enable_sync_message_emission (GstBus * bus);

GST_API
void                    gst_bus_disable_sync_message_emission (GstBus * bus);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstBus, gst_object_unref)

G_END_DECLS

#endif /* __GST_BUS_H__ */
