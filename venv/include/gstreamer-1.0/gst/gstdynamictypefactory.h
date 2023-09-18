/* GStreamer
 * Copyright (C) 2015 Jan Schmidt <jan@centricular.com>
 *
 * gstdynamictypefactory.h: Header for GstDynamicTypeFactory
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

#ifndef __GST_DYNAMIC_TYPE_FACTORY_H__
#define __GST_DYNAMIC_TYPE_FACTORY_H__

/**
 * GST_DYNAMIC_TYPE_REGISTER_DEFINE:
 * @t_n: The dynamic type name in lower case, with words separated by '_'.
 * Used to generate `gst_dynamic_type_register_*(GstPlugin* plugin)`.
 * @t: The #GType of the dynamic type

 * A convenience macro to define the entry point of a
 * dynamic type `gst_dynamic_type_register_*(GstPlugin* plugin)`.
 *
 * Since: 1.20
 */
#define GST_DYNAMIC_TYPE_REGISTER_DEFINE(t_n, t) \
G_BEGIN_DECLS \
gboolean G_PASTE (gst_dynamic_type_register_, t_n) (GstPlugin * plugin) \
{ \
  return gst_dynamic_type_register (plugin, t); \
} \
G_END_DECLS

/**
 * GST_DYNAMIC_TYPE_REGISTER_DECLARE:
 * @t_f: The dynamic type name in lower case, with words separated by '_'.
 *
 * This macro can be used to declare a new dynamic type.
 * It has to be used in combination with #GST_DYNAMIC_TYPE_REGISTER_DEFINE macro
 * and must be placed outside any block to declare the type find registration
 * function.
 *
 * Since: 1.20
 */
#define GST_DYNAMIC_TYPE_REGISTER_DECLARE(t_n) \
G_BEGIN_DECLS \
gboolean G_PASTE(gst_dynamic_type_register_, t_n) (GstPlugin * plugin); \
G_END_DECLS

/**
 * GST_DYNAMIC_TYPE_REGISTER:
 * @t_n: The dynamic type name to register
 * @plugin: The #GstPlugin where to register the dynamic type.
 *
 * This macro can be used to register a dynamic type into a #GstPlugin.
 * This method will be usually called in the plugin init function
 * but can also be called with a NULL plugin.
 *
 * Since: 1.20
 */
#define GST_DYNAMIC_TYPE_REGISTER(t_n, plugin) G_PASTE(gst_dynamic_type_register_, t_n) (plugin)

/**
 * GstDynamicTypeFactory:
 *
 * The opaque #GstDynamicTypeFactory data structure.
 *
 * Since: 1.12
 */
typedef struct _GstDynamicTypeFactory GstDynamicTypeFactory;
typedef struct _GstDynamicTypeFactoryClass GstDynamicTypeFactoryClass;

#include <gst/gstconfig.h>
#include <gst/gstplugin.h>
#include <gst/gstpluginfeature.h>

G_BEGIN_DECLS

#define GST_TYPE_DYNAMIC_TYPE_FACTORY           (gst_dynamic_type_factory_get_type())
#define GST_DYNAMIC_TYPE_FACTORY(obj)           (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_DYNAMIC_TYPE_FACTORY,\
                                                 GstDynamicTypeFactory))
#define GST_DYNAMIC_TYPE_CLASS(klass)           (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_DYNAMIC_TYPE_FACTORY,\
                                                 GstDynamicTypeFactoryClass))
#define GST_IS_DYNAMIC_TYPE_FACTORY(obj)        (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_DYNAMIC_TYPE_FACTORY))
#define GST_IS_DYNAMIC_TYPE_FACTORY_CLASS(klass)   (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_DYNAMIC_TYPE_FACTORY))
#define GST_DYNAMIC_TYPE_FACTORY_CAST(obj)      ((GstDynamicTypeFactory *)(obj))

GST_API
GType  gst_dynamic_type_factory_get_type        (void);

GST_API
GType  gst_dynamic_type_factory_load            (const gchar *factoryname);

GST_API
gboolean gst_dynamic_type_register              (GstPlugin *plugin, GType type);

G_END_DECLS

#endif
