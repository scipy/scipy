/* GStreamer
 *
 * Copyright (C) 2011 Stefan Sauer <ensonic@users.sf.net>
 *
 * gstdirectcontrolbinding.h: Direct attachment for control sources
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

#ifndef __GST_DIRECT_CONTROL_BINDING_H__
#define __GST_DIRECT_CONTROL_BINDING_H__

#include <gst/gstconfig.h>

#include <glib-object.h>

#include <gst/gstcontrolsource.h>
#include <gst/controller/controller-prelude.h>

G_BEGIN_DECLS

#define GST_TYPE_DIRECT_CONTROL_BINDING \
  (gst_direct_control_binding_get_type())
#define GST_DIRECT_CONTROL_BINDING(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_DIRECT_CONTROL_BINDING,GstDirectControlBinding))
#define GST_DIRECT_CONTROL_BINDING_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_CAST((klass),GST_TYPE_DIRECT_CONTROL_BINDING,GstDirectControlBindingClass))
#define GST_IS_DIRECT_CONTROL_BINDING(obj) \
  (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_DIRECT_CONTROL_BINDING))
#define GST_IS_DIRECT_CONTROL_BINDING_CLASS(klass) \
  (G_TYPE_CHECK_CLASS_TYPE((klass),GST_TYPE_DIRECT_CONTROL_BINDING))
#define GST_DIRECT_CONTROL_BINDING_GET_CLASS(obj) \
  (G_TYPE_INSTANCE_GET_CLASS ((obj), GST_TYPE_CONTOL_SOURCE, GstDirectControlBindingClass))

typedef struct _GstDirectControlBinding GstDirectControlBinding;
typedef struct _GstDirectControlBindingClass GstDirectControlBindingClass;

/**
 * GstDirectControlBindingConvertValue:
 * @self: the #GstDirectControlBinding instance
 * @src_value: the value returned by the cotnrol source
 * @dest_value: the target location
 *
 * Function to map a control-value to the target plain data type.
 */
typedef void (* GstDirectControlBindingConvertValue) (GstDirectControlBinding *self, gdouble src_value, gpointer dest_value);

/**
 * GstDirectControlBindingConvertGValue:
 * @self: the #GstDirectControlBinding instance
 * @src_value: the value returned by the cotnrol source
 * @dest_value: the target GValue
 *
 * Function to map a control-value to the target GValue.
 */
typedef void (* GstDirectControlBindingConvertGValue) (GstDirectControlBinding *self, gdouble src_value, GValue *dest_value);

/**
 * GstDirectControlBinding:
 * @name: name of the property of this binding
 *
 * The instance structure of #GstDirectControlBinding.
 */
struct _GstDirectControlBinding {
  GstControlBinding parent;
  
  /*< private >*/
  GstControlSource *cs;    /* GstControlSource for this property */
  GValue cur_value;
  gdouble last_value;
  gint byte_size;

  GstDirectControlBindingConvertValue convert_value;
  GstDirectControlBindingConvertGValue convert_g_value;

  union {
    gpointer _gst_reserved[GST_PADDING];
    struct {
      gboolean want_absolute;
    } abi;
  } ABI;
};

/**
 * GstDirectControlBindingClass:
 * @parent_class: Parent class
 * @convert: Class method to convert control-values
 *
 * The class structure of #GstDirectControlBinding.
 */

struct _GstDirectControlBindingClass
{
  GstControlBindingClass parent_class;

  /*< private >*/
  gpointer _gst_reserved[GST_PADDING];
};

GST_CONTROLLER_API
GType gst_direct_control_binding_get_type (void);

/* Functions */

GST_CONTROLLER_API
GstControlBinding * gst_direct_control_binding_new (GstObject * object, const gchar * property_name,
                                                    GstControlSource * cs);
GST_CONTROLLER_API
GstControlBinding * gst_direct_control_binding_new_absolute (GstObject * object, const gchar * property_name, 
                                                    GstControlSource * cs);

G_DEFINE_AUTOPTR_CLEANUP_FUNC(GstDirectControlBinding, gst_object_unref)

G_END_DECLS

#endif /* __GST_DIRECT_CONTROL_BINDING_H__ */
