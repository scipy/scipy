/* GStreamer
 * Copyright (C) 2012 Olivier Crete <olivier.crete@collabora.com>
 *
 * gstdeviceprovider.h: Device probing and monitoring
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

#ifndef __GST_DEVICE_PROVIDER_H__
#define __GST_DEVICE_PROVIDER_H__

#include <gst/gstelement.h>

/**
 * GST_DEVICE_PROVIDER_REGISTER_DEFINE_CUSTOM:
 * @d_p: The device provider name in lower case, with words separated by '_'.
 * Used to generate `gst_device_provider_register_*(GstPlugin* plugin)`.
 * @register_func: pointer to a method with the format: `gboolean register_func (GstPlugin* plugin);`
 *
 * A convenience macro to define the entry point of a
 * device provider `gst_device_provider_register_*(GstPlugin* plugin)` which uses
 * register_func as the main registration method for the device provider.
 * As an example, you may define the device provider named "device-provider"
 * with the namespace `my` as following using `device_provider_register_custom`:
 *
 * ```
 *
 * gboolean my_device_provider_register_custom (GstPlugin * plugin)
 * {
 *    gboolean ret = FALSE;
 *    ret |= gst_device_provider_register (plugin, "my-device-provider",
         GST_RANK_PRIMARY, GST_TYPE_MY_DEVICE_PROVIDER);
 *    return TRUE;
 * }
 *
 * GST_DEVICE_PROVIDER_REGISTER_DEFINE_CUSTOM (my_device_provider, my_device_provider_register_custom)
 * ```
 *
 * Since: 1.20
 */
#define GST_DEVICE_PROVIDER_REGISTER_DEFINE_CUSTOM(d_p, register_func) \
G_BEGIN_DECLS \
gboolean G_PASTE (gst_device_provider_register_, d_p) (GstPlugin * plugin) \
{ \
  return register_func (plugin); \
} \
G_END_DECLS

/**
 * GST_DEVICE_PROVIDER_REGISTER_DEFINE:
 * @d_p: The device provider name in lower case, with words separated by '_'.
 * Used to generate `gst_device_provider_register_*(GstPlugin* plugin)`.
 * @d_p_n: The public name of the device provider
 * @r: The #GstRank of the device provider (higher rank means more importance when autoplugging, see #GstRank)
 * @t: The #GType of the device provider.
 *
 * A convenience macro to define the entry point of a
 * device provider `gst_device_provider_register_*(GstPlugin* plugin)`.
 *
 * Since: 1.20
 */
#define GST_DEVICE_PROVIDER_REGISTER_DEFINE(d_p, d_p_n, r, t) \
G_BEGIN_DECLS \
gboolean G_PASTE (gst_device_provider_register_, d_p) (GstPlugin * plugin) \
{ \
  return gst_device_provider_register (plugin, d_p_n, r, t); \
} \
G_END_DECLS

/**
 * GST_DEVICE_PROVIDER_REGISTER_DECLARE:
 * @d_p: The device provider name in lower case, with words separated by '_'.
 *
 * This macro can be used to declare a new device provider.
 * It has to be used in combination with #GST_DEVICE_PROVIDER_REGISTER_DEFINE macro
 * and must be placed outside any block to declare the device provider registration
 * function.
 *
 * Since: 1.20
 */
#define GST_DEVICE_PROVIDER_REGISTER_DECLARE(d_p) \
G_BEGIN_DECLS \
gboolean G_PASTE(gst_device_provider_register_, d_p) (GstPlugin * plugin); \
G_END_DECLS

/**
 * GST_DEVICE_PROVIDER_REGISTER:
 * @d_p: The device provider name in lower case, with words separated by '_'.
 * @plugin: The #GstPlugin where to register the device provider.
 *
 * This macro can be used to register a device provider into a #GstPlugin.
 * This method will be usually called in the plugin init function
 * but can also be called with a NULL plugin.
 *
 * Since: 1.20
 */
#define GST_DEVICE_PROVIDER_REGISTER(d_p, plugin) G_PASTE(gst_device_provider_register_, d_p) (plugin)

G_BEGIN_DECLS

typedef struct _GstDeviceProvider GstDeviceProvider;
typedef struct _GstDeviceProviderClass GstDeviceProviderClass;
typedef struct _GstDeviceProviderPrivate GstDeviceProviderPrivate;

#include <gst/gstdeviceproviderfactory.h>

#define GST_TYPE_DEVICE_PROVIDER                 (gst_device_provider_get_type())
#define GST_IS_DEVICE_PROVIDER(obj)              (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GST_TYPE_DEVICE_PROVIDER))
#define GST_IS_DEVICE_PROVIDER_CLASS(klass)      (G_TYPE_CHECK_CLASS_TYPE ((klass), GST_TYPE_DEVICE_PROVIDER))
#define GST_DEVICE_PROVIDER_GET_CLASS(obj)       (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_DEVICE_PROVIDER, GstDeviceProviderClass))
#define GST_DEVICE_PROVIDER(obj)                 (G_TYPE_CHECK_INSTANCE_CAST ((obj), GST_TYPE_DEVICE_PROVIDER, GstDeviceProvider))
#define GST_DEVICE_PROVIDER_CLASS(klass)         (G_TYPE_CHECK_CLASS_CAST ((klass), GST_TYPE_DEVICE_PROVIDER, GstDeviceProviderClass))
#define GST_DEVICE_PROVIDER_CAST(obj)            ((GstDeviceProvider *)(obj))


/**
 * GstDeviceProvider:
 * @parent: The parent #GstObject
 * @devices: a #GList of the #GstDevice objects
 *
 * The structure of the base #GstDeviceProvider
 *
 * Since: 1.4
 */
struct _GstDeviceProvider {
  GstObject         parent;

  /* Protected by the Object lock */
  GList *devices;

  /*< private >*/

  GstDeviceProviderPrivate *priv;

  gpointer _gst_reserved[GST_PADDING];
};

/**
 * GstDeviceProviderClass:
 * @parent_class: the parent #GstObjectClass structure
 * @factory: a pointer to the #GstDeviceProviderFactory that creates this
 *  provider
 * @probe: Returns a list of devices that are currently available.
 *  This should never block. The devices should not have a parent and should
 *  be floating.
 * @start: Starts monitoring for new devices. Only subclasses that can know
 *  that devices have been added or remove need to implement this method.
 * @stop: Stops monitoring for new devices. Only subclasses that implement
 *  the start() method need to implement this method.
 *
 * The structure of the base #GstDeviceProviderClass
 *
 * Since: 1.4
 */

struct _GstDeviceProviderClass {
  GstObjectClass    parent_class;

  GstDeviceProviderFactory     *factory;

  GList*      (*probe) (GstDeviceProvider * provider);

  gboolean    (*start) (GstDeviceProvider * provider);
  void        (*stop)  (GstDeviceProvider * provider);

  /*< private >*/
  gpointer metadata;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_API
GType       gst_device_provider_get_type (void);


GST_API
GList *     gst_device_provider_get_devices    (GstDeviceProvider * provider);

GST_API
gboolean    gst_device_provider_start          (GstDeviceProvider * provider);

GST_API
void        gst_device_provider_stop           (GstDeviceProvider * provider);

GST_API
gboolean    gst_device_provider_can_monitor    (GstDeviceProvider * provider);

GST_API
GstBus *    gst_device_provider_get_bus        (GstDeviceProvider * provider);

GST_API
void        gst_device_provider_device_add     (GstDeviceProvider * provider,
                                                GstDevice * device);
GST_API
void        gst_device_provider_device_remove  (GstDeviceProvider * provider,
                                                GstDevice * device);
GST_API
gchar **    gst_device_provider_get_hidden_providers (GstDeviceProvider * provider);

GST_API
void        gst_device_provider_hide_provider        (GstDeviceProvider * provider,
                                                      const gchar       * name);
GST_API
void        gst_device_provider_unhide_provider      (GstDeviceProvider * provider,
                                                      const gchar       * name);

GST_API
const gchar * gst_device_provider_get_metadata       (GstDeviceProvider * provider,
                                                      const gchar * key);

GST_API
gboolean    gst_device_provider_is_started     (GstDeviceProvider * provider);

/* device provider class meta data */

GST_API
void        gst_device_provider_class_set_metadata         (GstDeviceProviderClass *klass,
                                                            const gchar     *longname,
                                                            const gchar     *classification,
                                                            const gchar     *description,
                                                            const gchar     *author);
GST_API
void        gst_device_provider_class_set_static_metadata  (GstDeviceProviderClass *klass,
                                                            const gchar     *longname,
                                                            const gchar     *classification,
                                                            const gchar     *description,
                                                            const gchar     *author);
GST_API
void        gst_device_provider_class_add_metadata         (GstDeviceProviderClass * klass,
                                                            const gchar * key, const gchar * value);
GST_API
void        gst_device_provider_class_add_static_metadata  (GstDeviceProviderClass * klass,
                                                            const gchar * key, const gchar * value);
GST_API
const gchar * gst_device_provider_class_get_metadata       (GstDeviceProviderClass * klass,
                                                            const gchar * key);

GST_API
void gst_device_provider_device_changed                    (GstDeviceProvider * provider,
                                                            GstDevice *device,
                                                            GstDevice *changed_device);

/* factory management */

GST_API
GstDeviceProviderFactory * gst_device_provider_get_factory (GstDeviceProvider * provider);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstDeviceProvider, gst_object_unref)

G_END_DECLS

#endif /* __GST_DEVICE_PROVIDER_H__ */
